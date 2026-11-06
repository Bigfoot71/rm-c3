module rm @private @feat(RM_TAA);

import std::math;

// ================================
// Private: TAA Definitions
// ================================

const int   TAA_SAMPLES   = 16;     // jitter sequence length
const float TAA_MAX_AGE   = 14.0f;  // maximum age of a frame
const float TAA_CLAMP_PAD = 0.01f;  // padding around the neighborhood box
const float TAA_REJECT    = 12.0f;  // distance to clamp * this factor = history decay

// ================================
// Private: TAA Functions
// ================================

fn void taa_init(int w, int h)
{
	sz n = (sz)w * h * float[<4>]::size;
	ctx.curr = ((float[<4>]*)mem::calloc_aligned(n, CACHE_LINE))[:w*h];
	ctx.hist = ((float[<4>]*)mem::calloc_aligned(n, CACHE_LINE))[:w*h];

	ctx.w = w;
	ctx.h = h;

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
	mem::free_aligned(ctx.hist);
	mem::free_aligned(ctx.curr);
}

fn void taa_reset() @inline
{
	ctx.force_reset = true;
}

fn Vec2 taa_jitter() @inline
{
	return ctx.jitter;
}

fn void taa_prepare(RmCamera camera)
{
	Vec3 fwd = camera.forward();
	Vec3 up  = camera.up();

	bool moved = ctx.force_reset
		|| (camera.position - ctx.prev_pos).length_sq() > 1e-10f
		|| (fwd - ctx.prev_fwd).length_sq()             > 1e-10f
		|| (up  - ctx.prev_up).length_sq()              > 1e-10f
		|| camera.fovy != ctx.prev_fovy
		|| camera.aspect != ctx.prev_aspect;

	ctx.still = moved ? 0 : math::min(ctx.still + 1, 1 << 20);
	ctx.force_reset = false;

	ctx.prev_pos = camera.position;  ctx.prev_fwd = fwd;  ctx.prev_up = up;
	ctx.prev_fovy = camera.fovy;     ctx.prev_aspect = camera.aspect;

	ctx.jitter = ctx.still == 0 ? (Vec2){ 0, 0 } : ctx.jitter_seq[(ctx.still - 1) % TAA_SAMPLES];
}

fn void taa_store(int x, int y, RmV3 c) @inline
{
	RmV3 s = c.saturate();
	int  i = y * ctx.w + x;

	for (int l = 0; l < LANES; l++)
	{
		ctx.curr[i + l] = { s.x[l], s.y[l], s.z[l], 0.0f };
	}
}

fn void taa_resolve_item(void* arg, int y)
{
	bool valid = ctx.still > 0;

	int ym = math::max(y - 1, 0);
	int yp = math::min(y + 1, ctx.h - 1);

	for (int x = 0; x < ctx.w; x += LANES)
	{
		RmV3 out;

		for (int l = 0; l < LANES; l++)
		{
			int px = x + l;
			int xm = math::max(px - 1, 0);
			int xp = math::min(px + 1, ctx.w - 1);
			int i  = y * ctx.w + px;

			float[<3>] c  = ctx.curr[i].xyz;
			float[<3>] lo = c;
			float[<3>] hi = c;

			int[3] rows = { ym * ctx.w, y * ctx.w, yp * ctx.w };
			int[3] cols = { xm, px, xp };

			foreach (r : rows)
			{
				foreach (cc : cols)
				{
					float[<3>] n = ctx.curr[r + cc].xyz;
					lo = math::min(lo, n);
					hi = math::max(hi, n);
				}
			}

			float[<4>] h   = ctx.hist[i];
			float[<3>] hc  = math::clamp(h.xyz, lo - TAA_CLAMP_PAD, hi + TAA_CLAMP_PAD);
			float[<3>] d   = h.xyz - hc;

			float rej = math::min(math::sqrt((d * d).sum()) * TAA_REJECT, 1.0f);
			float age = valid ? h[3] * (1.0f - rej) : 0.0f;

			float[<3>] res = hc + (c - hc) * (1.0f / (age + 1.0f));

			ctx.hist[i] = { res.x, res.y, res.z, math::min(age + 1.0f, TAA_MAX_AGE) };

			out.x[l] = res.x;
			out.y[l] = res.y;
			out.z[l] = res.z;
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

struct Context @local
{
	float[<4>][] curr;
	float[<4>][] hist;

	int w, h;

	Vec2[TAA_SAMPLES] jitter_seq;
	Vec2  jitter;
	int   still;
	Vec3  prev_pos, prev_fwd, prev_up;
	float prev_fovy, prev_aspect;
	bool  force_reset;
}

Context ctx @local;
