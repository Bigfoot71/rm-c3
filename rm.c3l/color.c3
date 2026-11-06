module rm @private;

import std::math;

// ================================
// Private: Color Functions
// ================================

fn RmV3 linear_to_srgb(RmV3 color) @inline @feat(RM_SRGB_OUTPUT)
{
	return {
		rmf_linear_to_srgb(color.x),
		rmf_linear_to_srgb(color.y),
		rmf_linear_to_srgb(color.z),
	};
}

fn RmV3 debanding(RmF x, RmF y, RmV3 color) @inline @feat(RM_DEBANDING)
{
	const RmF DITHER_STRENGTH = 1.0f / 255.0f;
	RmF d = (ign(x, y) - 0.5f) * DITHER_STRENGTH;
	return { color.x + d, color.y + d, color.z + d };
}

// ================================
// Local: Helper Functions
// ================================

fn RmF rmf_linear_to_srgb(RmF color) @local @inline @feat(RM_SRGB_OUTPUT)
{
    // Approximation from http://chilliant.blogspot.com/2012/08/srgb-approximations-for-hlsl.html
	return math::max(1.055f * math::pow(color, 0.416666667f) - 0.055f, (RmF)0.0f);
}

fn RmF rmf_srgb_to_linear(RmF color) @local @inline @feat(RM_SRGB_OUTPUT)
{
    // Approximation from http://chilliant.blogspot.com/2012/08/srgb-approximations-for-hlsl.html
    return color * (color * (color * 0.305306011 + 0.682171111) + 0.012522878);
}
