module rm @private @feat(RM_TAA);

import std::math;

// ================================
// Private: TAA Definitions
// ================================

const int   TAA_SAMPLES     = 16;      // jitter sequence length
const float TAA_AGE_MOVING  = 8.0f;    // history length while moving
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
	ctx.curr = ((float[<4>]*)mem::calloc_aligned(n, CACHE_LINE))[:w*h];
	ctx.hist = ((float[<4>]*)mem::calloc_aligned(n, CACHE_LINE))[:w*h];

	ctx.w = w;
	ctx.h = h;

	ctx.motion    = 0.0f;
	ctx.stability = 1.0f;
	ctx.age_cap   = 0.0f;
	ctx.pad       = TAA_CLAMP_PAD;
	ctx.reject    = TAA_REJECT;
	ctx.use_box   = true;

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
	Vec3 fwd = camera.forward();
	Vec3 up  = camera.up();

	bool moved = ctx.force_reset
		|| (camera.position - ctx.prev_pos).length_sq() > 1e-10f
		|| (fwd - ctx.prev_fwd).length_sq()             > 1e-10f
		|| (up  - ctx.prev_up).length_sq()              > 1e-10f
		|| camera.fovy != ctx.prev_fovy
		|| camera.aspect != ctx.prev_aspect;

	ctx.still = moved ? 0 : math::min(ctx.still + 1, 1 << 20);

	// Maximum history length for this frame
	float cap = TAA_AGE_STILL;
	if (ctx.force_reset) { cap = 0.0f; }
	else if (moved)      { cap = ctx.motion * TAA_AGE_MOVING; }
	ctx.age_cap = cap;

	// Trust in the history (0 right after a motion, 1 once a full jitter cycle has been accumulated at rest)
	float rest = math::min((float)ctx.still / (float)TAA_SAMPLES, 1.0f) * ctx.stability;
	ctx.pad     = lerp(TAA_CLAMP_PAD, TAA_REST_PAD, rest);
	ctx.reject  = lerp(TAA_REJECT, TAA_REST_REJECT, rest);
	ctx.use_box = rest < 1.0f;

	ctx.force_reset = false;

	ctx.prev_pos = camera.position;  ctx.prev_fwd = fwd;  ctx.prev_up = up;
	ctx.prev_fovy = camera.fovy;     ctx.prev_aspect = camera.aspect;

	// The jitter fades out with the history length (no history, no jitter)
	ctx.jitter = ctx.jitter_seq[ctx.frame] * math::min(cap, 1.0f);
	ctx.frame  = (ctx.frame + 1) % TAA_SAMPLES;
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
			float[<4>] h  = ctx.hist[i];
			float[<3>] hc = h.xyz;
			float      rej = 0.0f;

			// Once settled at rest the history is exact (no box, no rejection)
			if (ctx.use_box)
			{
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

				hc = math::clamp(h.xyz, lo - ctx.pad, hi + ctx.pad);

				float[<3>] d = h.xyz - hc;
				rej = math::min(math::sqrt((d * d).sum()) * ctx.reject, 1.0f);
			}

			float age = math::min(h[3], ctx.age_cap) * (1.0f - rej);

			float[<3>] res = hc + (c - hc) * (1.0f / (age + 1.0f));

			ctx.hist[i] = { res.x, res.y, res.z, math::min(age + 1.0f, ctx.age_cap) };

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
	int   frame;       // index in the jitter sequence
	int   still;       // consecutive frames without camera motion

	float motion;      // user: AA while moving, [0, 1]
	float stability;   // user: trust in the history at rest [0, 1]
	float age_cap;     // history length allowed this frame
	float pad;         // neighborhood box padding this frame
	float reject;      // rejection gain this frame
	bool  use_box;     // false = clamp / rejection skipped this frame

	Vec3  prev_pos, prev_fwd, prev_up;
	float prev_fovy, prev_aspect;
	bool  force_reset;
}

Context ctx @local;
