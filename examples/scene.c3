import std::math;
import raylib6;
import rm;

const WIN_W = 1280;
const WIN_H = 720;

// ================================
// Constants
// ================================

const int MAT_FLOOR  = 1;
const int MAT_MIRROR = 2;
const int MAT_RED    = 3;
const int MAT_GOLD   = 4;
const int MAT_BLUE   = 5;

const Vec3 SUN = { 0.4f, 0.8f, 0.45f };

// ================================
// Scene and Shading
// ================================

fn RmSd scene(RmV3 p)
{
	RmSd floor  = rm::sd_plane_y(p, 0.0f).material(MAT_FLOOR);
	RmSd mirror = rm::sd_sphere(p - rm::v3(0.0f, 1.0f, 0.0f), 1.0f).material(MAT_MIRROR);
	RmSd box    = rm::sd_round_box(p - rm::v3(-2.4f, 0.5f, 0.8f), { 0.5f, 0.5f, 0.5f }, 0.1f).material(MAT_RED);
	RmSd torus  = rm::sd_torus(p - rm::v3(2.4f, 0.25f, 0.4f), { 0.7f, 0.25f }).material(MAT_GOLD);
	RmSd caps   = rm::sd_vertical_capsule(p - rm::v3(0.6f, 0.0f, -2.2f), 1.2f, 0.3f).material(MAT_BLUE);

	return floor | mirror | box | torus | caps;
}

fn RmV3 shade(RmRays ray, RmHits hits, float cone)
{
	RmV3 sky = sky_color(ray.dir);
	if (!rm::mask_any(hits.mask)) return sky;

	RmSurface s = rm::get_surface(&scene, hits, cone);

	RmV3 albedo = albedo_of(s.mat, s.pos);
	RmF  rough  = rm::select(s.mat.comp_eq((RmI)MAT_MIRROR), (RmF)0.05f, (RmF)0.5f);
	RmF  f0     = rm::select(s.mat.comp_eq((RmI)MAT_MIRROR), (RmF)0.9f,  (RmF)0.04f);
	f0          = rm::select(s.mat.comp_eq((RmI)MAT_GOLD),   (RmF)0.5f,  f0);

	RmBrdf b = rm::get_brdf(s, ray, rough, f0);

	// Sun + soft shadows + AO + ambient
	RmF  shadow = rm::get_shadow(&scene, s, SUN, s.hit);
	RmF  ao     = rm::get_ambient_occlusion(&scene, s);
	RmLit lit   = rm::@light_dir(rm::diffuse_burley, rm::specular_ggx, b, SUN, rm::v3(1.0f, 0.9f, 0.75f) * (RmF)1.6f, shadow);
	RmV3 ambient = rm::get_ambient(s.n, rm::v3(0.25f, 0.35f, 0.55f), rm::v3(0.15f, 0.12f, 0.1f), ao);

	RmV3 col = lit.apply(albedo) + albedo * ambient;

	// Reflections
	RmReflect refl = rm::get_reflection(&scene, s, b, s.hit, cone * 4.0f);
	RmV3 rc = rm::@reflection(reflect_hit, sky_color, refl);
	col = col.lerp(rc, rm::get_reflectivity(b, ao));

	// Fog
	col = rm::fog_simple(col, rm::v3(0.75f, 0.85f, 0.95f), hits.dist, 0.02f);

	return rm::select3(hits.mask, col, sky);
}

// ================================
// Sky and Materials
// ================================

fn RmV3 sky_color(RmV3 dir)
{
	RmV3 c = rm::sky_gradient(dir, rm::v3(0.75f, 0.85f, 0.95f), rm::v3(0.15f, 0.35f, 0.75f), 0.6f);
	return c + rm::sky_sun(dir, SUN, rm::v3(1.0f, 0.85f, 0.6f));
}

fn RmV3 albedo_of(RmI mat, RmV3 p)
{
	// checkerboard floor
	RmF s   = math::floor(p.x) + math::floor(p.z);
	RmF chk = s - math::floor(s * 0.5f) * 2.0f; // 0 or 1
	RmV3 c  = rm::v3(0.1f, 0.1f, 0.1f).lerp(rm::v3(0.8f, 0.8f, 0.8f), chk);

	c = rm::select3(mat.comp_eq((RmI)MAT_MIRROR), rm::v3(0.02f, 0.02f, 0.02f), c);
	c = rm::select3(mat.comp_eq((RmI)MAT_RED),    rm::v3(0.8f, 0.1f, 0.1f),    c);
	c = rm::select3(mat.comp_eq((RmI)MAT_GOLD),   rm::v3(1.0f, 0.7f, 0.2f),    c);
	c = rm::select3(mat.comp_eq((RmI)MAT_BLUE),   rm::v3(0.1f, 0.3f, 0.9f),    c);
	return c;
}

// Minimal lighting for what is seen in the reflections
fn RmV3 reflect_hit(RmReflect r)
{
	RmF ndl = rm::saturate(r.n.dot(rm::v3_from(SUN)));
	return albedo_of(r.mat, r.pos) * (ndl * 0.8f + 0.2f);
}

// ================================
// Main
// ================================

fn int main(String[] args)
{
	rl::init_window(WIN_W, WIN_H, "RM Scene");
	defer rl::close_window();
	rl::set_target_fps(60);

	rm::init(WIN_W, WIN_H);
	defer rm::quit();

	RLTexture rt = {
		.id      = rl::gl_load_texture(null, WIN_W, WIN_H, (int)RLPixelFormat.UNCOMPRESSED_R8G8B8A8, 1),
		.width   = WIN_W,
		.height  = WIN_H,
		.format  = UNCOMPRESSED_R8G8B8A8,
		.mipmaps = 1,
	};
	defer rl::unload_texture(rt);

	RmCamera camera = rm::camera_new({ 0.0f, 2.5f, 7.0f }, { 0.0f, 0.8f, 0.0f }, (float)WIN_W / (float)WIN_H, 60.0f, 40.0f);

	while (!rl::window_should_close())
	{
		if (!rl::is_mouse_button_down(LEFT))
		{
			camera.orbit({ 0.0f, 0.8f, 0.0f }, rl::get_frame_time() * 0.4f, 0.0f);
		}

		char[<4>][] colors = rm::render(camera, &scene, &shade);
		rt.update(colors.ptr);

		rl::@drawing()
		{
			rl::clear_background({ 255, 0, 255, 255 });
			rl::draw_texture(rt, 0, 0, rl::WHITE);
			rl::draw_fps(10, 10);
		};
	}

	return 0;
}
