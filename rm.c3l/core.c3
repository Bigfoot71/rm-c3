module rm;

import std::math;

// ================================
// Public: Core Constants
// ================================

const bool RM_SRGB_OUTPUT = $feat(RM_SRGB_OUTPUT);
const bool RM_DEBANDING   = $feat(RM_DEBANDING);
const bool RM_THREADED    = $feat(RM_THREADED);
const bool RM_ANAA        = $feat(RM_ANAA);
const bool RM_TAA         = $feat(RM_TAA);

// ================================
// Private: Core Constants
// ================================

const int LANES       @private = 8;
const int TILE        @private = 16;
const int THREADS     @private = env::WASM ? 5 : 0; // 0 = one per logical CPU
const int CACHE_LINE  @private = 64;                // Color buffer alignment take a cache line so two threads never write the same one

const float BEAM_SAFETY @private = 0.98f;

$assert TILE % LANES == 0 : "TILE must be a multiple of LANES";

// ================================
// Public: Core Definitions
// ================================

alias RmShadeFn = fn RmV3(RmRays ray, RmHits hits, float cone);
alias RmSceneFn = fn RmSd(RmV3 p);

// ================================
// Public: Core Functions
// ================================

fn void init(int w, int h)
{
	w = (w + LANES - 1) & ~(LANES - 1);

	ctx.w       = w;
	ctx.h       = h;
	ctx.tiles_x = (w + TILE - 1) / TILE;
	ctx.tiles_y = (h + TILE - 1) / TILE;

	ctx.colors = ((char[<4>]*)mem::calloc_aligned((sz)w*h * char[<4>]::size, CACHE_LINE))[:w*h];
	ctx.beam   = mem::new_array(float, (sz)ctx.tiles_x * ctx.tiles_y);

	$if RM_THREADED:
		threads_init(THREADS);
	$endif

	$if RM_TAA:
		taa_init(w, h);
	$endif
}

fn void quit()
{
	$if RM_TAA:
		taa_quit();
	$endif

	$if RM_THREADED:
		threads_quit();
	$endif

	mem::free_aligned(ctx.colors);
	mem::free(ctx.beam);
}

fn char[<4>][] render(RmCamera camera, RmSceneFn scene, RmShadeFn shade, int max_steps = 128, float eps = 0.0005f, float omega0 = 1.6f)
{
	float tan_half = camera.tan_half();

	ctx.cam_pos = camera.position;
	ctx.fwd     = v3_from(camera.forward());
	ctx.right_s = v3_from(camera.right() * (camera.aspect * tan_half));
	ctx.up_s    = v3_from(camera.up()    * tan_half);

	$if RM_TAA:
		taa_prepare(camera);
	$endif

	float cone_x    = camera.aspect * tan_half / (float)ctx.w;
	float cone_y    = tan_half / (float)ctx.h;
	float cone      = math::max(cone_x, cone_y);
	float beam_cone = (float)TILE * math::sqrt(cone_x * cone_x + cone_y * cone_y) + cone;

	BeamArgs beam = {
		.scene     = scene,
		.beam_cone = beam_cone,
		.max_dist  = camera.far,
		.max_steps = max_steps >> 1,
		.eps       = eps,
	};

	TileArgs tile = {
		.shade     = shade,
		.scene     = scene,
		.cone      = cone,
		.max_dist  = camera.far,
		.max_steps = max_steps,
		.eps       = eps,
		.omega0    = omega0,
	};

	for_each(&beam_item, &beam, beam_item_count());
	for_each(&tile_item, &tile, tile_item_count()); // needs the beam pass to be complete

	$if RM_TAA:
		for_each(&taa_resolve_item, null, ctx.h);
	$endif

	return ctx.colors;
}

fn int[<2>] resolution() @inline
{
	return { ctx.w, ctx.h };
}

fn void reset_history() @inline
{
	$if RM_TAA:
		taa_reset();
	$endif
}

// motion:    0 = no AA while the camera moves, 1 = full accumulation (ghosting)
// stability: 1 = converges perfectly at rest (static scene), 0 = keeps rejecting changes (animated scene)
fn void set_taa(float motion, float stability) @inline
{
	$if RM_TAA:
		taa_configure(motion, stability);
	$endif
}

// ================================
// Local: Frame Functions
// ================================

// Generates rays from continuous pixel coordinates (centers at +0.5)
fn RmRays pixel_rays(RmF cx, float cy) @local @inline
{
	$if RM_TAA:
		cx += taa_jitter().x;
		cy += taa_jitter().y;
	$endif

	RmF   px = cx * (2.0f / (float)ctx.w) - 1.0f;
	float py = 1.0f - cy * (2.0f / (float)ctx.h);

	RmV3 dir = ctx.fwd + ctx.right_s * px + ctx.up_s * py;
	return { .origin = ctx.cam_pos, .dir = dir.normalize() };
}

fn void store_packet(int x, int y, RmV3 color) @private @inline
{
	RmI r = (RmI)(saturate(color.x) * 255.0f + 0.5f);
	RmI g = (RmI)(saturate(color.y) * 255.0f + 0.5f);
	RmI b = (RmI)(saturate(color.z) * 255.0f + 0.5f);

	RmI packed = r | (g << 8) | (b << 16) | (255 << 24);

	*(RmI*)&ctx.colors[y * ctx.w + x] = packed;
}

// ================================
// Private: Work Definitions
// ================================

alias ItemFn @private = fn void(void* arg, int item);

struct BeamArgs @private
{
	RmSceneFn scene;
	float beam_cone;
	float max_dist;
	int   max_steps;
	float eps;
}

struct TileArgs @private
{
	RmShadeFn shade;
	RmSceneFn scene;
	float cone;
	float max_dist;
	int   max_steps;
	float eps;
	float omega0;
}

// ================================
// Local: Work Functions
// ================================

// Runs item_fn(arg, i) for every i in 0 .. count - 1
fn void for_each(ItemFn item_fn, void* arg, int count) @local @inline
{
	$if RM_THREADED:
		threads_run(item_fn, arg, count);
	$else
		for (int i = 0; i < count; i++)
		{
			item_fn(arg, i);
		}
	$endif
}

// LANES neighbouring tiles of one tile row (one beam packet)
fn int beam_item_count() @local
{
	return ctx.tiles_y * ((ctx.tiles_x + LANES - 1) / LANES);
}

fn void beam_item(void* arg, int item) @local
{
	BeamArgs* args = arg;

	int groups = (ctx.tiles_x + LANES - 1) / LANES;
	int ty     = item / groups;
	int tx     = (item % groups) * LANES;

	float cy = (float)(ty * TILE) + (float)TILE * 0.5f;
	RmF   cx = (RmF)(float)(tx * TILE) + (float)TILE * 0.5f + lane_index() * (float)TILE;
	RmF   t0 = march_beam(args.scene, pixel_rays(cx, cy), args.beam_cone, args.max_dist, args.max_steps, args.eps, BEAM_SAFETY);

	for (int l = 0; l < LANES; l++)
	{
		if (tx + l >= ctx.tiles_x) break;
		ctx.beam[ty * ctx.tiles_x + tx + l] = t0[l];
	}
}

fn int tile_item_count() @local
{
	return ctx.tiles_x * ctx.tiles_y;
}

fn void tile_item(void* arg, int item) @local
{
	TileArgs* args = arg;

	int tx = item % ctx.tiles_x;
	int ty = item / ctx.tiles_x;

	int x0 = tx * TILE;
	int x1 = math::min(x0 + TILE, ctx.w);
	int y0 = ty * TILE;
	int y1 = math::min(y0 + TILE, ctx.h);

	float t0 = ctx.beam[item];

	for (int y = y0; y < y1; y++)
	{
		for (int x = x0; x < x1; x += LANES)
		{
			RmF xl = (RmF)x + lane_index();

			RmRays rays = pixel_rays(xl + 0.5f, (float)y + 0.5f);
			RmHits hits = march(args.scene, rays, args.cone, t0, args.max_dist, args.max_steps, args.eps, args.omega0);
			RmV3  color = args.shade(rays, hits, args.cone);

			$if RM_ANAA:
				color = anaa_resolve(color, rays, hits, args);
			$endif

			$if RM_TAA:
				taa_store(x, y, color, hits.dist * rays.dir.dot(ctx.fwd));
			$else
				$if RM_SRGB_OUTPUT:
					color = linear_to_srgb(color);
				$endif
				$if RM_DEBANDING:
					color = debanding(xl, (RmF)y, color);
				$endif
				store_packet(x, y, color);
			$endif
		}
	}
}

// ================================
// Local: Core Context
// ================================

struct Context @local
{
	char[<4>][] colors;
	float[]     beam;     // safe start distance per tile

	int w, h;
	int tiles_x, tiles_y;

	Vec3  cam_pos;
	RmV3  fwd, right_s, up_s;
	float tan_half;
	float aspect;
}

Context ctx @local;
