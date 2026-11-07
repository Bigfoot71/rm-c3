module rm @private @feat(RM_ANAA);

// ================================
// Private: ANAA Definitions
// ================================

const float ANAA_MIN_COV = 1.0f / 255; // below this the edge is not worth a second shading
const float ANAA_HIT     = 0.25f;      // hit when within this fraction of the pixel footprint (t * cone)
const float ANAA_RAMP    = 1.5f;       // coverage falls to 0 at this many footprints from the surface

// ================================
// Local: ANAA Functions
// ================================

fn RmV3 anaa_resolve(RmV3 color, RmRays rays, RmHits hits, TileArgs* args) @inline
{
	RmM edge = hits.gz_cov.comp_gt((RmF)ANAA_MIN_COV);

	if (mask_any(edge))
	{
		RmHits gh = {
			.mask = edge,
			.dist = hits.gz_dist,
			.d    = hits.gz_d,
			.mat  = hits.gz_mat,
			.pos  = rays.at(hits.gz_dist),
		};

		RmV3 gc    = args.shade(rays, gh, args.cone);
		RmF  c     = hits.gz_cov;
		RmV3 mixed = color.lerp(gc, c);

		color = select3(edge, mixed, color);
	}

	return color;
}
