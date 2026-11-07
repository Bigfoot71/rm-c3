module rm;

import std::math;

// ================================
// Public: Fresnel
// ================================

fn RmF fresnel_schlick(RmF cos_theta, RmF f0) @inline
{
	return f0 + ((RmF)1.0f - f0) * pow5((RmF)1.0f - saturate(cos_theta));
}

fn RmF fresnel_schlick_rough(RmF cos_theta, RmF f0, RmF rough) @inline
{
	RmF top = math::max((RmF)1.0f - rough, f0);
	return f0 + (top - f0) * pow5((RmF)1.0f - saturate(cos_theta));
}

// ================================
// Public: BRDF Inputs
// ================================

struct RmBrdf
{
	RmV3 n;
	RmV3 v;     // towards the viewer
	RmF  ndv;
	RmF  rough; // perceptual roughness (0..1)
	RmF  a;     // alpha = rough^2 (clamped for highlights)
	RmF  f0;    // reflectance at normal incidence (0.04 for dielectrics)
}

fn RmBrdf get_brdf(RmSurface s, RmRays ray, RmF rough, RmF f0) @inline
{
	RmV3 v = -ray.dir;
	RmF  r = saturate(rough);

	return {
		.n     = s.n,
		.v     = v,
		.ndv   = math::max(s.n.dot(v), (RmF)1e-4f),
		.rough = r,
		.a     = math::max(r * r, (RmF)0.02f),
		.f0    = f0,
	};
}

fn RmV3 get_reflect_dir(RmBrdf b) @inline
{
	return (-b.v).reflect(b.n);
}

// ================================
// Public: Diffuse Terms
// ================================

fn RmF diffuse_lambert(RmBrdf b, RmV3 l) @inline
{
	return saturate(b.n.dot(l));
}

fn RmF diffuse_burley(RmBrdf b, RmV3 l) @inline
{
	RmV3 h    = (b.v + l).normalize();
	RmF  ndl  = saturate(b.n.dot(l));
	RmF  ldh  = saturate(l.dot(h));
	RmF  fd90 = (RmF)0.5f + ldh * ldh * b.rough * 2.0f;
	RmF  fl   = pow5((RmF)1.0f - ndl);
	RmF  fv   = pow5((RmF)1.0f - b.ndv);

	return ((RmF)1.0f + (fd90 - 1.0f) * fl) * ((RmF)1.0f + (fd90 - 1.0f) * fv) * ndl;
}

// ================================
// Public: Specular Terms
// ================================

fn RmF specular_blinn_phong(RmBrdf b, RmV3 l) @inline
{
	RmV3 h   = (b.v + l).normalize();
	RmF  ndl = saturate(b.n.dot(l));
	RmF  ndh = saturate(b.n.dot(h));
	RmF  vdh = saturate(b.v.dot(h));

	RmF a2   = b.a * b.a;
	RmF p    = (RmF)2.0f / a2 - 2.0f;
	RmF norm = (p + 8.0f) * 0.125f;

	return pow(ndh, p) * norm * fresnel_schlick(vdh, b.f0) * ndl;
}

// GGX (Trowbridge-Reitz) + height-correlated Smith visibility (approx) + Schlick
fn RmF specular_ggx(RmBrdf b, RmV3 l) @inline
{
	RmV3 h   = (b.v + l).normalize();
	RmF  ndl = saturate(b.n.dot(l));
	RmF  ndh = saturate(b.n.dot(h));
	RmF  vdh = saturate(b.v.dot(h));

	RmF a2 = b.a * b.a;
	RmF dd = ndh * ndh * (a2 - 1.0f) + 1.0f;
	RmF dist = a2 / (dd * dd);

	RmF vis_den = lerp(ndl * b.ndv * 2.0f, ndl + b.ndv, b.a);
	RmF vis     = (RmF)0.5f / math::max(vis_den, (RmF)1e-5f);

	return dist * vis * fresnel_schlick(vdh, b.f0) * ndl;
}

// ================================
// Public: Lights
// ================================

struct RmLit
{
	RmV3 diff;
	RmV3 spec;
}

fn RmLit RmLit.add(self, RmLit o) @operator(+) @inline
{
	return { self.diff + o.diff, self.spec + o.spec };
}

fn RmV3 RmLit.apply(self, RmV3 albedo) @inline
{
	return self.diff * albedo + self.spec;
}

fn RmF point_atten(RmF d2, float range) @inline
{
	RmF r2 = (RmF)(range * range);
	RmF x  = d2 / r2;
	RmF q  = saturate((RmF)1.0f - x * x);
	return q * q / (d2 + 1.0f);
}

// Usage: light_dir(diffuse_burley, specular_ggx, b, dir, color, shadow)
macro RmLit @light_dir(#diff, #spec, RmBrdf b, Vec3 dir, RmV3 color, RmF shadow)
{
	RmV3 l = v3_from(dir);

	return { color * (#diff(b, l) * shadow), color * (#spec(b, l) * shadow) };
}

macro RmLit @light_point(#diff, #spec, RmBrdf b, RmV3 pos, Vec3 light_pos, RmV3 color, float range, RmF shadow)
{
	RmV3 to = v3_from(light_pos) - pos;
	RmF  d2 = math::max(to.length_sq(), (RmF)1e-6f);
	RmV3 l  = to * ((RmF)1.0f / math::sqrt(d2));
	RmF  k  = point_atten(d2, range) * shadow;

	return { color * (#diff(b, l) * k), color * (#spec(b, l) * k) };
}

macro RmLit @light_spot(#diff, #spec, RmBrdf b, RmV3 pos, Vec3 light_pos, Vec3 axis, float cos_inner, float cos_outer, RmV3 color, float range, RmF shadow)
{
	RmV3 to = v3_from(light_pos) - pos;
	RmF  d2 = math::max(to.length_sq(), (RmF)1e-6f);
	RmV3 l  = to * ((RmF)1.0f / math::sqrt(d2));

	RmF cd   = -l.dot(v3_from(axis));
	RmF cone = smoothstep(cos_outer, cos_inner, cd);
	RmF k    = point_atten(d2, range) * cone * shadow;

	return { color * (#diff(b, l) * k), color * (#spec(b, l) * k) };
}

// ================================
// Public: Shadows
// ================================

// Soft shadow toward a (uniform) light direction
fn RmF get_shadow(RmSceneFn scene, RmSurface s, Vec3 light, RmM live, float tmin = 0.005f, float tmax = 20.0f, float k = 8.0f, int max_steps = 96) @inline
{
	if (!mask_any(live)) return (RmF)1.0f;

	RmV3 origin = s.pos + s.n * math::max((RmF)0.004f, s.foot * 2.0f);
	RmV3 dir    = v3_from(light);

	RmF res    = (RmF)1.0f;
	RmF t      = (RmF)tmin;
	RmM active = live;

	for (int i = 0; i < max_steps; i++)
	{
		RmF h = scene(origin + dir * t).d;

		RmM blocked = active & h.comp_lt((RmF)0.001f);
		res    = select(blocked, (RmF)0.0f, res);
		active = active & mask_not(blocked);

		res    = select(active, math::min(res, h * k / t), res);
		t      = select(active, t + h, t);
		active = active & t.comp_lt((RmF)tmax);

		if (!mask_any(active)) break;
	}
	return res;
}

// ================================
// Public: Ambient
// ================================

fn RmV3 get_ambient(RmV3 n, RmV3 sky, RmV3 ground, RmF ao) @inline
{
	RmF t = n.y * 0.5f + 0.5f;
	return ground.lerp(sky, t) * ao;
}

// 1 = open, 0 = fully occluded
fn RmF get_ambient_occlusion(RmSceneFn scene, RmSurface s, float strength = 3.0f, int max_steps = 3) @inline
{
	RmF occ = (RmF)0.0f;
	float sca = 1.0f;

	for (int i = 0; i < max_steps; i++)
	{
		float h = 0.01f + 0.12f * (float)i / (float)(max_steps - 1);
		RmF   d = scene(s.pos + s.n * h).d;
		occ += ((RmF)h - d) * sca;
		sca *= 0.95f;
	}
	return saturate((RmF)1.0f - occ * strength);
}

fn RmF get_specular_occlusion(RmBrdf b, RmF ao) @inline
{
	RmF e = math::exp((b.rough * (-16.0f) - 1.0f) * 0.6931472f); // 2^(-16 * rough - 1)
	return saturate(pow(b.ndv + ao, e) - 1.0f + ao);
}

fn RmF get_reflectivity(RmBrdf b, RmF ao) @inline
{
	return fresnel_schlick_rough(b.ndv, b.f0, b.rough) * get_specular_occlusion(b, ao);
}

// ================================
// Public: Sky
// ================================

fn RmV3 sky_gradient(RmV3 dir, RmV3 horizon, RmV3 top, float sharpness = 1.0f) @inline
{
	RmF t = saturate(dir.y);
	if (sharpness != 1.0f) t = pow(t, (RmF)sharpness);
	return horizon.lerp(top, t);
}

fn RmV3 sky_sun(RmV3 dir, Vec3 sun_dir, RmV3 color, float disc_cos = 0.9995f, float halo_power = 48.0f, float disc_gain = 4.0f) @inline
{
	RmF c    = saturate(dir.dot(v3_from(sun_dir)));
	RmF disc = smoothstep(disc_cos - (1.0f - disc_cos) * 0.5f, disc_cos, c);
	RmF halo = pow(c, (RmF)halo_power) * 0.25f;

	return color * (disc * disc_gain + halo);
}

// ================================
// Public: Fog
// ================================

fn RmV3 fog_simple(RmV3 col, RmV3 fog_col, RmF dist, float density) @inline
{
	RmF f = (RmF)1.0f - math::exp(dist * (-density));
	return col.lerp(fog_col, f);
}

fn RmV3 fog_full(RmV3 col, RmV3 fog_col, RmV3 sun_col, RmRays ray, RmF dist, Vec3 sun_dir, float density, float falloff, float sun_power = 8.0f) @inline
{
	RmF dy    = ray.dir.y;
	RmF x     = math::clamp(dist * dy * falloff, (RmF)(-30.0f), (RmF)30.0f);
	RmM small = comp_lt(math::abs(dy), 1e-4f);

	// (1 - exp(-x)) / dy; with its limit dist * falloff when dy -> 0
	RmF safe_dy  = select(small, (RmF)1.0f, dy);
	RmF integral = select(small, dist * falloff, ((RmF)1.0f - math::exp(-x)) / safe_dy);

	RmF amount = integral * (density * math::exp(-ray.origin.y * falloff));
	RmF f      = (RmF)1.0f - math::exp(-amount);

	RmF  sd   = saturate(ray.dir.dot(v3_from(sun_dir)));
	RmV3 tint = fog_col.lerp(sun_col, pow(sd, (RmF)sun_power));

	return col.lerp(tint, f);
}

// ================================
// Public: Emissive & Glow
// ================================

fn RmV3 add_emissive(RmV3 lit, RmV3 emission, RmF strength) @inline
{
	return lit + emission * strength;
}

fn RmF get_glow(RmSceneFn emitters, RmRays ray, RmF t_end, float falloff = 6.0f, int max_steps = 24, float min_step = 0.03f, float max_step = 0.5f) @inline
{
	RmF t      = (RmF)0.0f;
	RmF acc    = (RmF)0.0f;
	RmM active = comp_lt(t, t_end);

	for (int i = 0; i < max_steps; i++)
	{
		RmF d  = math::max(emitters(ray.at(t)).d, (RmF)0.0f);
		RmF dt = math::clamp(d, (RmF)min_step, (RmF)max_step);
		RmF w  = math::exp(d * (-falloff));

		acc    = select(active, acc + w * dt, acc);
		t      = select(active, t + dt, t);
		active = active & t.comp_lt(t_end);

		if (!mask_any(active)) break;
	}
	return acc;
}

// ================================
// Public: Translucency
// ================================

// 0 = solid, 1 = thin; samples the SDF inward along -n
fn RmF get_thinness(RmSceneFn scene, RmSurface s, float depth = 0.4f, int steps = 4) @inline
{
	RmF sum = (RmF)0.0f;

	for (int i = 0; i < steps; i++)
	{
		float h = depth * (float)(i + 1) / (float)steps;
		RmF   d = scene(s.pos - s.n * h).d;

		sum += saturate((d + h) * (1.0f / h));
	}
	return sum * (1.0f / (float)steps);
}

// Light bleeding through a thin object (Barre-Brisebois)
fn RmV3 get_translucency(RmBrdf b, Vec3 light, RmV3 color, RmF thin, float distortion = 0.3f, float power = 4.0f, float scale = 1.0f) @inline
{
	RmV3 l   = v3_from(light);
	RmV3 h   = (l + b.n * distortion).normalize();
	RmF  vdh = saturate(-b.v.dot(h));

	return color * (pow(vdh, (RmF)power) * thin * scale);
}

// ================================
// Public: Reflections
// ================================

struct RmReflect
{
	RmM  hit;
	RmV3 pos;
	RmV3 n;
	RmV3 dir;
	RmI  mat;
	RmF  dist;
}

// Use a bigger `cone` for blurrier / cheaper reflections
fn RmReflect get_reflection(RmSceneFn scene, RmSurface s, RmBrdf b, RmM live, float cone, float max_dist = 20.0f, int max_steps = 48, float eps = 0.001f) @inline
{
	RmV3 dir    = get_reflect_dir(b);
	RmV3 origin = s.pos + s.n * math::max((RmF)0.004f, s.foot * 2.0f);

	RmReflect out = { .dir = dir, .dist = (RmF)max_dist };
	if (!mask_any(live)) return out;

	RmF t      = (RmF)0.02f;
	RmM active = live;

	for (int i = 0; i < max_steps; i++)
	{
		RmSd sd = scene(origin + dir * t);
		RmF  r  = sd.d;

		RmM got = active & r.comp_lt(math::max((RmF)eps, t * cone));

		out.hit  = out.hit | got;
		out.dist = select(got, t, out.dist);
		out.mat  = select(got, sd.mat, out.mat);
		active   = active & mask_not(got);

		t      = select(active, t + r, t);
		active = active & t.comp_lt((RmF)max_dist);

		if (!mask_any(active)) break;
	}

	out.pos = origin + dir * out.dist;

	if (mask_any(out.hit))
	{
		RmF neps = math::clamp(out.dist * cone, (RmF)0.001f, (RmF)0.01f);
		out.n    = get_normal(scene, out.pos, neps);
	}
	return out;
}

// Resolves the reflected color
macro RmV3 @reflection(#on_hit, #on_miss, RmReflect r)
{
	RmV3 c = v3(0.0f, 0.0f, 0.0f);

	if (!mask_all(r.hit)) c = #on_miss(r.dir);
	if (mask_any(r.hit))  c = select3(r.hit, #on_hit(r), c);

	return c;
}

// ================================
// Local: Math Helpers
// ================================

fn RmF pow5(RmF x) @local @inline
{
	RmF x2 = x * x;
	return x2 * x2 * x;
}

fn RmF smoothstep(float e0, float e1, RmF x) @local @inline
{
	RmF t = saturate((x - (RmF)e0) * (RmF)(1.0f / (e1 - e0)));
	return t * t * ((RmF)3.0f - t * 2.0f);
}
