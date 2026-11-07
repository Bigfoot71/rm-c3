module rm;

import std::math;

// ================================
// Public: Surface Definitions
// ================================

struct RmSurface
{
	RmM  hit;
	RmV3 pos;  // projected onto the surface
	RmV3 n;
	RmI  mat;
	RmF  dist;
	RmF  foot; // pixel footprint at the hit (cone * dist)
}

// ================================
// Public: Surface Functions
// ================================

// Normal + point projected onto the surface
fn RmSurface get_surface(RmSceneFn scene, RmHits h, float cone) @inline
{
	RmF  foot = h.dist * cone;
	RmF  neps = math::clamp(foot, (RmF)0.001f, (RmF)0.01f);
	RmV3 n    = get_normal(scene, h.pos, neps);
	RmV3 ps   = h.pos - n * h.d;
	//   n    = get_normal(scene, ps, neps);

	return { .hit = h.mask, .pos = ps, .n = n, .mat = h.mat, .dist = h.dist, .foot = foot };
}

fn RmV3 get_normal(RmSceneFn scene, RmV3 p, RmF eps) @inline
{
	RmSd s0 = scene({ p.x + eps, p.y - eps, p.z - eps });
	RmSd s1 = scene({ p.x - eps, p.y - eps, p.z + eps });
	RmSd s2 = scene({ p.x - eps, p.y + eps, p.z - eps });
	RmSd s3 = scene({ p.x + eps, p.y + eps, p.z + eps });

	RmV3 n = {
		 s0.d - s1.d - s2.d + s3.d,
		-s0.d - s1.d + s2.d + s3.d,
		-s0.d + s1.d - s2.d + s3.d,
	};
	return n.normalize();
}

fn RmV3 get_heat(RmF t) @inline
{
	t = saturate(t);
	return {
		saturate(t * 2.0f - 0.5f),
		(RmF)1.0f - math::abs(t * 2.0f - 1.0f),
		saturate((RmF)1.0f - t * 2.0f + 0.5f),
	};
}
