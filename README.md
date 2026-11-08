# rm-c3

Software raymarcher running on the CPU, written in [C3](https://c3-lang.org/).

It can even run in a browser, I made a mini game compiled to WASM which you can try here:
- Link => https://bigfoot71.github.io/demos/rm-racer

## Example

By default `c3c run` launches the example scene. It contains a few primitives with reflections, soft shadows, AO, sky and fog.

You can also hold the left mouse button to pause the camera and see the TAA in action.

Everything is rendered at 1280x720 with all "rm's features" enabled, so it can be a bit slow for real time.

The example uses raylib to display the rendered image, but raylib is only required by the example, not by the lib itself.

## Usage
```c3
rm::init(w, h);
char[<4>][] colors = rm::render(camera, &scene, &shade);
rm::quit();
```

- `scene(RmV3 p) -> RmSd`: the signed distance and material at a given point.
- `shade(RmRays, RmHits, float cone) -> RmV3`: the color of each ray packet.

## Features

They are configured in `project.json` and are all enabled by default.

| Flag | Effect |
|---|---|
| `RM_THREADED` | Multithreading |
| `RM_ANAA` | Analytic anti-aliasing |
| `RM_TAA` | Temporal anti-aliasing |
| `RM_SRGB_OUTPUT` | Linear to sRGB output conversion |
| `RM_DEBANDING` | Debanding (dithering) |

## Optimizations

- SIMD rendering in ray packets, with masks for active or terminated rays
- 16x16 tile-based rendering and cache-line aligned buffers to avoid false sharing between threads
- Beam pass: a coarse cone per tile providing a safe starting distance
- Over-relaxed sphere tracing with backtracking
- Adaptive intersection threshold based on distance (pixel cone)

## Contents

- `camera`: quaternion-based camera (movement, orbit, look-at)
- `sdf`: primitives, operators (`|`, `&`, `-`, smooth variants) and modifiers (`material`, `grow`, `shell`)
- `surface`: normals and surface projection
- `shading`: PBR BRDF (Burley / GGX), directional, point and spot lights, soft shadows, AO, reflections, sky, fog, emissive and glow, translucency
- Exposed SIMD math functions to help with shading
