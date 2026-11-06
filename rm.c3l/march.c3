module rm;

import std::math;

// ================================
// Private: March Constants
// ================================

const float ANAA_HIT_EFF @local @feat( RM_ANAA) = ANAA_HIT;
const float ANAA_HIT_EFF @local @feat(!RM_ANAA) = 1.0f;

// ================================
// Public: March Definitions
// ================================

struct RmRays
{
	Vec3 origin;
	RmV3 dir;
}

struct RmHits
{
	RmM  mask;
	RmF  dist;  // distance along the rays
	RmF  d;     // SDF value at the hit (to project onto the surface)
	RmI  mat;   // material id at the hit
	RmI  steps; // iterations used (debug)
	RmV3 pos;

	RmF gz_cov  @feat(RM_ANAA); // coverage of the grazed surface, 0 = none
	RmF gz_dist @feat(RM_ANAA); // distance along the ray at the closest approach
	RmF gz_d    @feat(RM_ANAA); // SDF value there
	RmI gz_mat  @feat(RM_ANAA);
}

// ================================
// Public: Rays Methods
// ================================

fn RmV3 RmRays.at(self, RmF t) @inline
{
	return self.dir * t + v3_from(self.origin);
}

// ================================
// Private: March Functions
// ================================

fn RmHits march(RmSceneFn scene, RmRays rays, float cone, float t0, float max_dist, int max_steps, float eps, float omega0) @private
{
	const float BIG = 1.0e9f;

	RmF t      = (RmF)math::min(t0, max_dist);
	RmF step   = (RmF)0.0f;
	RmF prev_r = (RmF)0.0f;
	RmF omega  = (RmF)omega0;
	RmM active = t.comp_lt((RmF)max_dist);

	RmHits out = { .dist = (RmF)max_dist };

$if RM_ANAA:
	RmF cur_r = (RmF)BIG;  RmF cur_t = (RmF)1.0f;  RmI cur_m = (RmI)0;
	RmF gr    = (RmF)BIG;  RmF gt    = (RmF)1.0f;  RmI gm    = (RmI)0;
$endif

	for (int i = 0; i < max_steps && mask_any(active); i++)
	{
		RmSd s    = scene(rays.at(t));
		RmF  r    = s.d;
		RmF  foot = t * cone;

		RmM back = omega.comp_gt((RmF)1.0f) & (r + prev_r).comp_lt(step);
		RmM got  = active & mask_not(back) & r.comp_lt(math::max((RmF)eps, foot * ANAA_HIT_EFF));

	$if RM_ANAA:
		RmM live = active & mask_not(back) & mask_not(got);

		RmM closer = live & (r * cur_t).comp_lt(cur_r * t);
		RmM leave  = live & mask_not(closer) & cur_r.comp_lt(cur_t * (ANAA_RAMP * cone));
		RmM better = leave & (cur_r * gt).comp_lt(gr * cur_t);

		gr = select(better, cur_r, gr);
		gt = select(better, cur_t, gt);
		gm = select(better, cur_m, gm);

		RmM restart = closer | leave;
		cur_r = select(restart, r,     cur_r);
		cur_t = select(restart, t,     cur_t);
		cur_m = select(restart, s.mat, cur_m);
	$endif

		out.mask  = out.mask | got;
		out.dist  = select(got, t, out.dist);
		out.d     = select(got, r, out.d);
		out.mat   = select(got, s.mat, out.mat);
		out.steps = select(active, (RmI)(i + 1), out.steps);
		active    = active & mask_not(got);

		RmF newstep = r * omega;

		t      = select(back, t + prev_r - step, t + newstep);
		step   = select(back, prev_r, newstep);
		prev_r = select(back, prev_r, r);
		omega  = select(back, (RmF)1.0f, omega);

		active = active & t.comp_lt((RmF)max_dist);
	}

$if RM_ANAA:
	RmM tail = mask_not(out.mask) & cur_r.comp_lt(cur_t * (ANAA_RAMP * cone)) & (cur_r * gt).comp_lt(gr * cur_t);

	gr = select(tail, cur_r, gr);
	gt = select(tail, cur_t, gt);
	gm = select(tail, cur_m, gm);

	RmF gq = gr / math::max(gt * cone, (RmF)1e-6f);

	out.gz_cov  = saturate(((RmF)ANAA_RAMP - gq) * (1.0f / (ANAA_RAMP - ANAA_HIT)));
	out.gz_dist = gt;
	out.gz_d    = gr;
	out.gz_mat  = gm;
$endif

	out.pos = rays.at(out.dist);
	return out;
}

// Marches a cone covering the tile and returns safe start distances per ray
fn RmF march_beam(RmSceneFn scene, RmRays rays, float beam_cone, float max_dist, int max_steps, float eps, float safety) @private
{
	RmF  t      = (RmF)0.0f;
	RmF  result = (RmF)max_dist;
	RmM  active = (RmM)true;
	RmF  inv    = (RmF)(1.0f / (1.0f + beam_cone));

	for (int i = 0; i < max_steps && mask_any(active); i++)
	{
		RmF r      = scene(rays.at(t)).d;
		RmF margin = t * beam_cone + eps;

		RmM stop = active & r.comp_le(margin);
		result   = select(stop, t * safety, result);
		active   = active & mask_not(stop);

		t += (r - margin) * inv;
		active = active & t.comp_lt((RmF)max_dist);
	}

	return select(active, t * safety, result);
}
