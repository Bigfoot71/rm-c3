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
	ctx.w      = w;
	ctx.h      = h;
	ctx.stride = w + LANES;

	// One block of 3 buffers (curr, hist, hist_prev) x 4 planes (r, g, b, a)
	sz     n = (sz)ctx.stride * h + 2 * LANES;
	float* m = mem::calloc_aligned(n * 12 * float::size, CACHE_LINE);
	m += LANES;

	ctx.curr      = { m,         m + n,      m + 2 * n,  m + 3 * n  };
	ctx.hist      = { m + 4 * n, m + 5 * n,  m + 6 * n,  m + 7 * n  };
	ctx.hist_prev = { m + 8 * n, m + 9 * n,  m + 10 * n, m + 11 * n };

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
	mem::free_aligned(ctx.curr.r - LANES);
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
	// Last frame output is this frame input
	@swap(ctx.hist, ctx.hist_prev);

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

	// No history after a reset, a short one while moving, a long one at rest
	ctx.age_cap   = ctx.force_reset ? 0.0f : moved ? ctx.motion * TAA_AGE_MOVING : TAA_AGE_STILL;
	ctx.reproject = moved && ctx.age_cap > 0.0f;

	// Goes from 0 right after a motion to 1 after a full jitter cycle at rest
	ctx.rest = (float)ctx.still / (float)TAA_SAMPLES * ctx.stability;

	// Fades out with the history length
	ctx.jitter = ctx.jitter_seq[ctx.frame] * math::min(ctx.age_cap, 1.0f);
	ctx.frame  = (ctx.frame + 1) % TAA_SAMPLES;

	// Previous camera position is d + z * (f + gx * ndc_x + gy * ndc_y), with x and y already divided by the previous ax and ay
	Vec3 pr = ctx.prev.right / ctx.prev.ax;
	Vec3 pu = ctx.prev.up    / ctx.prev.ay;
	Vec3 pf = ctx.prev.fwd;
	Vec3 gx = ctx.cur.right * ctx.cur.ax;
	Vec3 gy = ctx.cur.up    * ctx.cur.ay;

	ctx.rp.d  = { ctx.dpos.dot(pr),    ctx.dpos.dot(pu),    ctx.dpos.dot(pf) };
	ctx.rp.f  = { ctx.cur.fwd.dot(pr), ctx.cur.fwd.dot(pu), ctx.cur.fwd.dot(pf) };
	ctx.rp.gx = { gx.dot(pr),          gx.dot(pu),          gx.dot(pf) };
	ctx.rp.gy = { gy.dot(pr),          gy.dot(pu),          gy.dot(pf) };

	ctx.force_reset = false;
}

// z = view depth of the sample (distance along the camera forward axis)
fn void taa_store(int x, int y, RmV3 c, RmF z) @inline
{
	RmV3 s = c.saturate();
	int  i = y * ctx.stride + x;

	*(RmF*)&ctx.curr.r[i] = s.x;
	*(RmF*)&ctx.curr.g[i] = s.y;
	*(RmF*)&ctx.curr.b[i] = s.z;
	*(RmF*)&ctx.curr.a[i] = z;

	if (x == 0)
	{
		ctx.curr.r[i - 1] = s.x[0];
		ctx.curr.g[i - 1] = s.y[0];
		ctx.curr.b[i - 1] = s.z[0];
		ctx.curr.a[i - 1] = z[0];
	}

	if (x + LANES == ctx.w)
	{
		ctx.curr.r[i + LANES] = s.x[LANES - 1];
		ctx.curr.g[i + LANES] = s.y[LANES - 1];
		ctx.curr.b[i + LANES] = s.z[LANES - 1];
		ctx.curr.a[i + LANES] = z[LANES - 1];
	}
}

fn void taa_resolve_item(void* arg, int y)
{
	// Once settled at rest the history is exact (no box, no rejection)
	bool boxed  = ctx.rest < 1.0f;
	RmF  pad    = (RmF)lerp(TAA_CLAMP_PAD, TAA_REST_PAD,    ctx.rest);
	RmF  reject = (RmF)lerp(TAA_REJECT,    TAA_REST_REJECT, ctx.rest);
	RmF  cap    = (RmF)ctx.age_cap;
	RmF  one    = (RmF)1.0f;
	RmF  zero   = (RmF)0.0f;

	int[3] rows = {
		math::max(y - 1, 0) * ctx.stride,
		y * ctx.stride,
		math::min(y + 1, ctx.h - 1) * ctx.stride,
	};

	// Row constants of the reprojection
	float ndc_y = 1.0f - ((float)y + 0.5f) * (2.0f / (float)ctx.h);
	Vec3  ky    = ctx.rp.f + ctx.rp.gy * ndc_y;
	float hw    = 0.5f * (float)ctx.w;
	float hh    = 0.5f * (float)ctx.h;

	for (int x = 0; x < ctx.w; x += LANES)
	{
		int i = y * ctx.stride + x;

		RmV3 c = {
			*(RmF*)&ctx.curr.r[i],
			*(RmF*)&ctx.curr.g[i],
			*(RmF*)&ctx.curr.b[i],
		};

		RmV3 lo = c;
		RmV3 hi = c;
		RmF  zn = *(RmF*)&ctx.curr.a[i];

		// Neighborhood box + closest depth
		if (boxed)
		{
			for (int r = 0; r < 3; r++)
			{
				for (int dx = -1; dx <= 1; dx++)
				{
					int j = rows[r] + x + dx;

					RmV3 n = {
						mem::load((RmF*)&ctx.curr.r[j], 4),
						mem::load((RmF*)&ctx.curr.g[j], 4),
						mem::load((RmF*)&ctx.curr.b[j], 4),
					};

					lo = lo.vmin(n);
					hi = hi.vmax(n);
					zn = math::min(zn, mem::load((RmF*)&ctx.curr.a[j], 4));
				}
			}
		}

		RmV3 h;
		RmF  hage;

		// Reprojected while the camera moves, same pixel otherwise
		if (ctx.reproject)
		{
			RmF ndc_x = ((RmF)x + lane_index() + 0.5f) * (2.0f / (float)ctx.w) - 1.0f;

			// Previous camera space position
			RmF px = (ndc_x * ctx.rp.gx[0] + ky[0]) * zn + ctx.rp.d[0];
			RmF py = (ndc_x * ctx.rp.gx[1] + ky[1]) * zn + ctx.rp.d[1];
			RmF pz = (ndc_x * ctx.rp.gx[2] + ky[2]) * zn + ctx.rp.d[2];

			RmM front = comp_ge(pz, 1e-4f);
			RmF inv   = one / math::max(pz, (RmF)1e-4f);

			// Continuous coordinates in pixel (center) space
			RmF fx = math::clamp(px * inv * hw + (hw - 0.5f),         (RmF)-2.0f, (RmF)((float)ctx.w + 1.0f));
			RmF fy = math::clamp((hh - 0.5f) - py * inv * hh,         (RmF)-2.0f, (RmF)((float)ctx.h + 1.0f));

			RmF flx = fx.floor();
			RmF fly = fy.floor();
			RmF tx  = fx - flx;
			RmF ty  = fy - fly;
			RmI x0  = (RmI)flx;
			RmI y0  = (RmI)fly;

			RmV3 acc;
			RmF  aacc;
			RmF  wsum;

			for (int jy = 0; jy < 2; jy++)
			{
				for (int jx = 0; jx < 2; jx++)
				{
					RmI sx = x0 + jx;
					RmI sy = y0 + jy;

					RmM ok = front
						& comp_ge(sx, 0) & comp_lt(sx, ctx.w)
						& comp_ge(sy, 0) & comp_lt(sy, ctx.h);

					RmF wt = select(ok, (jx == 0 ? one - tx : tx) * (jy == 0 ? one - ty : ty), zero);

					RmI idx = math::clamp(sy, (RmI)0, (RmI)(ctx.h - 1)) * ctx.stride
					        + math::clamp(sx, (RmI)0, (RmI)(ctx.w - 1));

					acc.x += gather(ctx.hist_prev.r, idx) * wt;
					acc.y += gather(ctx.hist_prev.g, idx) * wt;
					acc.z += gather(ctx.hist_prev.b, idx) * wt;
					aacc  += gather(ctx.hist_prev.a, idx) * wt;
					wsum  += wt;
				}
			}

			// Color renormalized over the valid taps
			RmM good = comp_ge(wsum, 1e-3f);
			RmF k    = select(good, one / math::max(wsum, (RmF)1e-3f), zero);

			h    = acc * k;
			hage = select(good, aacc, zero);
		}
		else
		{
			h = {
				*(RmF*)&ctx.hist_prev.r[i],
				*(RmF*)&ctx.hist_prev.g[i],
				*(RmF*)&ctx.hist_prev.b[i],
			};
			hage = *(RmF*)&ctx.hist_prev.a[i];
		}

		// At rest the box is padded out (TAA_REST_PAD) and the rejection is null
		RmV3 hc  = h.clamp(lo - pad, hi + pad);
		RmF  rej = math::min((h - hc).length() * reject, one);

		RmF  age = math::min(hage, cap) * (one - rej);
		RmV3 res = hc + (c - hc) * (one / (age + one));

		*(RmF*)&ctx.hist.r[i] = res.x;
		*(RmF*)&ctx.hist.g[i] = res.y;
		*(RmF*)&ctx.hist.b[i] = res.z;
		*(RmF*)&ctx.hist.a[i] = math::min(age + one, cap);

		RmV3 out = res;

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

struct Planes @local
{
	float* r, g, b, a; // a = depth (curr) / age (hist)
}

struct TaaView @local
{
	Vec3  pos, fwd, right, up;
	float ax, ay; // image plane half size at depth 1
}

struct TaaReproj @local
{
	Vec3 d, f, gx, gy; // (x, y, z) = previous right / up / forward axes
}

struct Context @local
{
	Planes curr;      // rgb + depth
	Planes hist;      // rgb + age (written)
	Planes hist_prev; // rgb + age (read)

	int w, h, stride;

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

	TaaView   cur, prev;
	TaaReproj rp;
	Vec3      dpos;    // cur.pos - prev.pos
}

Context ctx @local;
