module rm @private @feat(RM_TAA);

import std::math;

// ================================
// Private: TAA Definitions
// ================================

const int   TAA_SAMPLES     = 16;      // jitter sequence length
const float TAA_AGE_MOVING  = 16.0f;   // history length while moving (at motion = 1)
const float TAA_AGE_STILL   = 128.0f;  // history length at rest
const float TAA_CLAMP_PAD   = 0.02f;   // padding around the neighborhood box
const float TAA_REJECT      = 8.0f;    // distance to the box * this factor = history decay
const float TAA_REST_PAD    = 1.0f;    // same once settled at rest (1 = clamp disabled)
const float TAA_REST_REJECT = 0.0f;    // same once settled at rest (0 = never reject)

// ================================
// Private: TAA Functions
// ================================

fn void taa_init(int w, int h)
{
	sz n = (sz)w * h * float[<4>]::size;
	ctx.curr      = ((float[<4>]*)mem::calloc_aligned(n, CACHE_LINE))[:w*h];
	ctx.hist      = ((float[<4>]*)mem::calloc_aligned(n, CACHE_LINE))[:w*h];
	ctx.hist_prev = ((float[<4>]*)mem::calloc_aligned(n, CACHE_LINE))[:w*h];

	ctx.w = w;
	ctx.h = h;

	ctx.motion      = 0.5f;
	ctx.stability   = 1.0f;
	ctx.force_reset = true; // no previous view yet

	for (int i = 0; i < TAA_SAMPLES; i++)
	{
		ctx.jitter_seq[i] = {
			halton(i + 1, 2) - 0.5f,
			halton(i + 1, 3) - 0.5f,
		};
	}
}

fn void taa_quit()
{
	mem::free_aligned(ctx.hist_prev);
	mem::free_aligned(ctx.hist);
	mem::free_aligned(ctx.curr);
	ctx = {};
}

fn void taa_reset() @inline
{
	ctx.force_reset = true;
}

fn void taa_configure(float motion, float stability)
{
	ctx.motion    = saturate(motion);
	ctx.stability = saturate(stability);
}

fn Vec2 taa_jitter() @inline
{
	return ctx.jitter;
}

fn void taa_prepare(RmCamera camera)
{
	// The history written last frame is read this frame
	float[<4>][] tmp = ctx.hist;
	ctx.hist         = ctx.hist_prev;
	ctx.hist_prev    = tmp;

	float tan_half = camera.tan_half();

	ctx.prev = ctx.cur;
	ctx.cur  = {
		.pos   = camera.position,
		.fwd   = camera.forward(),
		.right = camera.right(),
		.up    = camera.up(),
		.ax    = camera.aspect * tan_half,
		.ay    = tan_half,
	};
	ctx.dpos = ctx.cur.pos - ctx.prev.pos;

	bool moved = ctx.force_reset
		|| ctx.dpos.length_sq() > 1e-10f
		|| (ctx.cur.fwd - ctx.prev.fwd).length_sq() > 1e-10f
		|| (ctx.cur.up  - ctx.prev.up).length_sq()  > 1e-10f
		|| ctx.cur.ax != ctx.prev.ax
		|| ctx.cur.ay != ctx.prev.ay;

	ctx.still = moved ? 0 : math::min(ctx.still + 1, TAA_SAMPLES);

	// Maximum history length for this frame
	ctx.age_cap   = ctx.force_reset ? 0.0f : moved ? ctx.motion * TAA_AGE_MOVING : TAA_AGE_STILL;
	ctx.reproject = moved && ctx.age_cap > 0.0f; // no history, nothing to reproject

	// Trust in the history (0 right after a motion, 1 once a full jitter cycle has been accumulated at rest)
	ctx.rest = (float)ctx.still / (float)TAA_SAMPLES * ctx.stability;

	// The jitter fades out with the history length (no history, no jitter)
	ctx.jitter = ctx.jitter_seq[ctx.frame] * math::min(ctx.age_cap, 1.0f);
	ctx.frame  = (ctx.frame + 1) % TAA_SAMPLES;

	ctx.force_reset = false;
}

// z = view depth of the sample (distance along the camera forward axis)
fn void taa_store(int x, int y, RmV3 c, RmF z) @inline
{
	RmV3 s = c.saturate();
	int  i = y * ctx.w + x;

	for (int l = 0; l < LANES; l++)
	{
		ctx.curr[i + l] = { s.x[l], s.y[l], s.z[l], z[l] };
	}
}

fn void taa_resolve_item(void* arg, int y)
{
	// Once settled at rest the history is exact (no box, no rejection)
	bool  boxed  = ctx.rest < 1.0f;
	float pad    = lerp(TAA_CLAMP_PAD, TAA_REST_PAD,    ctx.rest);
	float reject = lerp(TAA_REJECT,    TAA_REST_REJECT, ctx.rest);

	for (int x = 0; x < ctx.w; x += LANES)
	{
		RmV3 out;

		for (int l = 0; l < LANES; l++)
		{
			int px = x + l;
			int i  = y * ctx.w + px;

			float[<4>] c4  = ctx.curr[i];
			float[<3>] c   = c4.xyz;
			float[<4>] h   = ctx.hist_prev[i];
			float[<3>] hc  = h.xyz;
			float      rej = 0.0f;

			if (boxed)
			{
				float[<3>] lo = c;
				float[<3>] hi = c;
				float      zn = c4[3];

				// Neighborhood box + closest depth
				for (int dy = -1; dy <= 1; dy++)
				{
					int row = math::clamp(y + dy, 0, ctx.h - 1) * ctx.w;

					for (int dx = -1; dx <= 1; dx++)
					{
						float[<4>] n = ctx.curr[row + math::clamp(px + dx, 0, ctx.w - 1)];
						lo = math::min(lo, n.xyz);
						hi = math::max(hi, n.xyz);
						zn = math::min(zn, n[3]);
					}
				}

				// Reprojected while the camera moves, same pixel otherwise
				if (ctx.reproject) h = taa_history(px, y, zn);

				hc = math::clamp(h.xyz, lo - pad, hi + pad);

				float[<3>] d = h.xyz - hc;
				rej = math::min(math::sqrt((d * d).sum()) * reject, 1.0f);
			}

			float age = math::min(h[3], ctx.age_cap) * (1.0f - rej);

			float[<3>] res = hc + (c - hc) * (1.0f / (age + 1.0f));

			ctx.hist[i] = { res.x, res.y, res.z, math::min(age + 1.0f, ctx.age_cap) };

			out.set_lane(l, res);
		}

		$if RM_SRGB_OUTPUT:
			out = linear_to_srgb(out);
		$endif

		$if RM_DEBANDING:
			out = debanding((RmF)x + lane_index(), (RmF)y, out);
		$endif

		store_packet(x, y, out);
	}
}

// ================================
// Local: TAA Functions
// ================================

fn float[<4>] taa_history(int px, int py, float z) @local
{
	const float[<4>] NONE = { 0.0f, 0.0f, 0.0f, 0.0f };

	float ndc_x = ((float)px + 0.5f) * (2.0f / (float)ctx.w) - 1.0f;
	float ndc_y = 1.0f - ((float)py + 0.5f) * (2.0f / (float)ctx.h);

	// World point = current camera + ray * z
	Vec3 dir = ctx.cur.fwd + ctx.cur.right * (ctx.cur.ax * ndc_x) + ctx.cur.up * (ctx.cur.ay * ndc_y);
	Vec3 v   = ctx.dpos + dir * z;

	// Into the previous camera
	float zp = v.dot(ctx.prev.fwd);
	if (zp < 1e-4f) return NONE;

	float inv = 1.0f / zp;
	float nx  = v.dot(ctx.prev.right) * inv / ctx.prev.ax;
	float ny  = v.dot(ctx.prev.up)    * inv / ctx.prev.ay;

	// Continuous coordinates in pixel (center) space
	float fx = math::clamp((nx + 1.0f) * 0.5f * (float)ctx.w - 0.5f, -2.0f, (float)ctx.w + 1.0f);
	float fy = math::clamp((1.0f - ny) * 0.5f * (float)ctx.h - 0.5f, -2.0f, (float)ctx.h + 1.0f);

	float flx = math::floor(fx);
	float fly = math::floor(fy);
	int   x0  = (int)flx;
	int   y0  = (int)fly;

	float[2] wx = { 1.0f - (fx - flx), fx - flx };
	float[2] wy = { 1.0f - (fy - fly), fy - fly };

	float[<4>] acc  = NONE;
	float      wsum = 0.0f;

	for (int j = 0; j < 2; j++)
	{
		for (int i = 0; i < 2; i++)
		{
			int sx = x0 + i;
			int sy = y0 + j;
			if ((uint)sx >= (uint)ctx.w || (uint)sy >= (uint)ctx.h) continue;

			float wt = wx[i] * wy[j];

			acc  += ctx.hist_prev[sy * ctx.w + sx] * wt;
			wsum += wt;
		}
	}

	if (wsum < 1e-3f) return NONE;

	// Color renormalized over the valid taps
	float k = 1.0f / wsum;
	return acc * (float[<4>]){ k, k, k, 1.0f };
}

fn float halton(int index, int base) @local
{
	float f = 1.0f;
	float r = 0.0f;
	while (index > 0)
	{
		f /= (float)base;
		r += f * (float)(index % base);
		index /= base;
	}
	return r;
}

// ================================
// Local: TAA Context
// ================================

struct TaaView @local
{
	Vec3  pos, fwd, right, up;
	float ax, ay; // image plane half size at depth 1
}

struct Context @local
{
	float[<4>][] curr;      // rgb + depth
	float[<4>][] hist;      // rgb + age (written)
	float[<4>][] hist_prev; // rgb + age (read)

	int w, h;

	Vec2[TAA_SAMPLES] jitter_seq;
	Vec2  jitter;
	int   frame;       // jitter index
	int   still;       // frames without motion

	float motion;      // user: AA while moving [0, 1]
	float stability;   // user: trust at rest [0, 1]
	float age_cap;     // max history length
	float rest;        // 0 = just moved, 1 = settled
	bool  reproject;
	bool  force_reset;

	TaaView cur, prev;
	Vec3    dpos;      // cur.pos - prev.pos
}

Context ctx @local;
