module rm;

import std::math;

// ================================
// Public: Math Constants
// ================================

const PI      = math::PI;
const TAU     = math::TAU;
const INV_TAU = 1.0 / TAU;
const HALF_PI = math::HALF_PI;

// ================================
// Public: Lane Types
// ================================

alias RmF = float [<LANES>]; // one float per lane
alias RmI = int   [<LANES>]; // one int per lane
alias RmM = bool  [<LANES>]; // one bool per lane (comparison result)

// ================================
// Public: Lane Helpers
// ================================

// Per-lane choice: m ? a : b  (see select2 and select3 for vectors)
macro select(RmM m, a, b) => math::select(m, a, b);

macro RmM comp_lt(a, b) => a.comp_lt(($Typeof(a))b);
macro RmM comp_le(a, b) => a.comp_le(($Typeof(a))b);
macro RmM comp_eq(a, b) => a.comp_eq(($Typeof(a))b);
macro RmM comp_gt(a, b) => a.comp_gt(($Typeof(a))b);
macro RmM comp_ge(a, b) => a.comp_ge(($Typeof(a))b);
macro RmM comp_ne(a, b) => a.comp_ne(($Typeof(a))b);

fn bool mask_any(RmM m) @inline => m.or();
fn bool mask_all(RmM m) @inline => m.and();
fn RmM  mask_not(RmM m) @inline => m ^ (RmM)true;

macro RmF lane_index()
{
	RmF v;
	$for int $i = 0; $i < LANES; $i++:
		v[$i] = (float)$i;
	$endfor

	return v;
}

// ================================
// Public: Math Helpers
// ================================

macro pow(x, p)
{
	return math::exp(math::ln(math::max(x, ($Typeof(x))1e-6f)) * p);
}

macro sin(x)
{
	var t = x * ($Typeof(x))INV_TAU;
	t = t - math::floor(t + 0.5f);
	var y = t * ($Typeof(x))TAU; // [-pi, pi)
	var s = y * (1.2732395f - 0.4052847f * math::abs(y));
	return s * (0.225f * (math::abs(s) - 1.0f) + 1.0f);
}

macro cos(x)
{
	return sin(x + ($Typeof(x))HALF_PI);
}

macro saturate(x)
{
	return math::clamp(x, ($Typeof(x))0, ($Typeof(x))1);
}

macro fract(x)
{
	return x - math::floor(x);
}

macro lerp(a, b, t)
{
	return a + (b - a) * t;
}

macro ign(x, y)
{
	const MX = ($Typeof(x))0.06711056;
	const MY = ($Typeof(x))0.00583715;
	const MZ = ($Typeof(x))52.9829189;

	return fract(MZ * fract(x * MX + y * MY));
}

// ================================
// Public: RmV2 Definitions
// ================================

struct RmV2
{
	RmF x, y;
}

fn RmV2 v2(float x, float y) @inline
{
	return { (RmF)x, (RmF)y };
}

fn RmV2 v2_splat(RmF s) @inline
{
	return { s, s };
}

fn RmV2 v2_from(Vec2 v) @inline
{
	return { (RmF)v.x, (RmF)v.y };
}

// ================================
// Public: RmV2 Operators
// ================================

fn RmV2 RmV2.add(self, RmV2 o) @operator(+) @inline
{
	return { self.x + o.x, self.y + o.y };
}

fn RmV2 RmV2.sub(self, RmV2 o) @operator(-) @inline
{
	return { self.x - o.x, self.y - o.y };
}

fn RmV2 RmV2.mul(self, RmV2 o) @operator(*) @inline
{
	return { self.x * o.x, self.y * o.y };
}

fn RmV2 RmV2.div(self, RmV2 o) @operator(/) @inline
{
	return { self.x / o.x, self.y / o.y };
}

fn RmV2 RmV2.neg(self) @operator(-) @inline
{
	return { -self.x, -self.y };
}

fn RmV2 RmV2.add_each(self, RmF s) @operator_s(+) @inline
{
	return { self.x + s, self.y + s };
}

fn RmV2 RmV2.sub_each(self, RmF s) @operator(-) @inline
{
	return { self.x - s, self.y - s };
}

fn RmV2 RmV2.sub_each_r(self, RmF s) @operator_r(-) @inline
{
	return { s - self.x, s - self.y };
}

fn RmV2 RmV2.scale(self, RmF s) @operator_s(*) @inline
{
	return { self.x * s, self.y * s };
}

fn RmV2 RmV2.div_each(self, RmF s) @operator(/) @inline
{
	return { self.x / s, self.y / s };
}

fn RmV2 RmV2.div_each_r(self, RmF s) @operator_r(/) @inline
{
	return { s / self.x, s / self.y };
}

fn RmV2 RmV2.add_eachf(self, float s) @operator_s(+) @inline
{
	return self.add_each((RmF)s);
}

fn RmV2 RmV2.sub_eachf(self, float s) @operator(-) @inline
{
	return self.sub_each((RmF)s);
}

fn RmV2 RmV2.sub_each_rf(self, float s) @operator_r(-) @inline
{
	return self.sub_each_r((RmF)s);
}

fn RmV2 RmV2.scalef(self, float s) @operator_s(*) @inline
{
	return self.scale((RmF)s);
}

fn RmV2 RmV2.div_eachf(self, float s) @operator(/) @inline
{
	return self.div_each((RmF)s);
}

fn RmV2 RmV2.div_each_rf(self, float s) @operator_r(/) @inline
{
	return self.div_each_r((RmF)s);
}

// ================================
// Public: RmV2 Comparison
// ================================

fn RmM RmV2.comp_eq(self, RmV2 o) @inline
{
	return self.x.comp_eq(o.x) & self.y.comp_eq(o.y);
}

fn RmM RmV2.comp_ne(self, RmV2 o) @inline
{
	return self.x.comp_ne(o.x) | self.y.comp_ne(o.y);
}

fn RmM RmV2.is_normalized(self, float eps = 1e-4f) @inline
{
	return comp_lt(math::abs(self.length_sq() - 1.0f), eps);
}

// ================================
// Public: RmV2 Reductions
// ================================

fn RmF RmV2.sum(self) @inline
{
	return self.x + self.y;
}

fn RmF RmV2.product(self) @inline
{
	return self.x * self.y;
}

fn RmF RmV2.max(self) @inline
{
	return math::max(self.x, self.y);
}

fn RmF RmV2.min(self) @inline
{
	return math::min(self.x, self.y);
}

// ================================
// Public: RmV2 Geometry
// ================================

fn RmF RmV2.dot(self, RmV2 o) @inline
{
	return self.x * o.x + self.y * o.y;
}

fn RmF RmV2.cross(self, RmV2 o) @inline
{
	return self.x * o.y - self.y * o.x;
}

fn RmF RmV2.length_sq(self) @inline
{
	return self.dot(self);
}

fn RmF RmV2.length(self) @inline
{
	return math::sqrt(self.length_sq());
}

fn RmF RmV2.distance_sq(self, RmV2 o) @inline
{
	return (self - o).length_sq();
}

fn RmF RmV2.distance(self, RmV2 o) @inline
{
	return (self - o).length();
}

fn RmV2 RmV2.normalize(self) @inline
{
	RmF inv = (RmF)1.0f / math::sqrt(math::max(self.length_sq(), (RmF)1e-20f));
	return self * inv;
}

fn RmV2 RmV2.perpendicular(self) @inline
{
	return { -self.y, self.x };
}

fn RmV2 RmV2.reflect(self, RmV2 n) @inline
{
	return self - n * (self.dot(n) * 2.0f);
}

macro RmV2 RmV2.refract(self, RmV2 n, eta)
{
	RmF e   = (RmF)eta;
	RmF ndi = n.dot(self);
	RmF k   = (RmF)1.0f - e * e * ((RmF)1.0f - ndi * ndi);
	RmV2 r  = self * e - n * (e * ndi + math::sqrt(math::max(k, (RmF)0.0f)));

	return select2(comp_ge(k, 0.0f), r, v2(0.0f, 0.0f));
}

fn RmV2 RmV2.rejection(self, RmV2 o) @inline
{
	return self - o * (self.dot(o) / o.dot(o));
}

macro RmV2 RmV2.rotate(self, angle)
{
	RmF a = (RmF)angle;
	RmF c = cos(a);
	RmF s = sin(a);

	return { self.x * c - self.y * s, self.x * s + self.y * c };
}

macro RmV2 RmV2.lerp(self, RmV2 b, t)
{
	return self + (b - self) * t;
}

macro RmV2 RmV2.clamp_length(self, lo, hi)
{
	RmF len = self.length();
	RmF tgt = math::clamp(len, (RmF)lo, (RmF)hi);

	return select2(comp_gt(len, 0.0f), self * (tgt / len), self);
}

macro RmV2 RmV2.move_towards(self, RmV2 target, max_distance)
{
	RmV2 d = target - self;
	RmF  l = d.length();
	RmF  m = (RmF)max_distance;

	return select2(comp_le(l, m) | comp_lt(l, 1e-8f), target, self + d * (m / l));
}

// ================================
// Public: RmV2 Per-Component
// ================================

fn RmV2 RmV2.abs(self) @inline
{
	return { math::abs(self.x), math::abs(self.y) };
}

fn RmV2 RmV2.ceil(self) @inline
{
	return { self.x.ceil(), self.y.ceil() };
}

fn RmV2 RmV2.floor(self) @inline
{
	return { self.x.floor(), self.y.floor() };
}

fn RmV2 RmV2.round(self) @inline
{
	return { self.x.round(), self.y.round() };
}

fn RmV2 RmV2.trunc(self) @inline
{
	return { self.x.trunc(), self.y.trunc() };
}

fn RmV2 RmV2.fract(self) @inline
{
	return { fract(self.x), fract(self.y) };
}

fn RmV2 RmV2.saturate(self) @inline
{
	return { saturate(self.x), saturate(self.y) };
}

fn RmV2 RmV2.copysign(self, RmV2 sgn) @inline
{
	return { self.x.copysign(sgn.x), self.y.copysign(sgn.y) };
}

fn RmV2 RmV2.fma(self, RmV2 b, RmV2 c) @inline
{
	return { self.x.fma(b.x, c.x), self.y.fma(b.y, c.y) };
}

macro RmV2 RmV2.pow(self, e)
{
	return { math::pow(self.x, (RmF)e), math::pow(self.y, (RmF)e) };
}

fn RmV2 RmV2.vmin(self, RmV2 o) @inline
{
	return { math::min(self.x, o.x), math::min(self.y, o.y) };
}

fn RmV2 RmV2.vmax(self, RmV2 o) @inline
{
	return { math::max(self.x, o.x), math::max(self.y, o.y) };
}

macro RmV2 RmV2.min_each(self, s)
{
	return { math::min(self.x, (RmF)s), math::min(self.y, (RmF)s) };
}

macro RmV2 RmV2.max_each(self, s)
{
	return { math::max(self.x, (RmF)s), math::max(self.y, (RmF)s) };
}

fn RmV2 RmV2.clamp(self, RmV2 lo, RmV2 hi) @inline
{
	return self.vmax(lo).vmin(hi);
}

macro RmV2 RmV2.clamp_each(self, lo, hi)
{
	return { math::clamp(self.x, (RmF)lo, (RmF)hi), math::clamp(self.y, (RmF)lo, (RmF)hi) };
}

// ================================
// Public: RmV2 Conversion
// ================================

fn RmV3 RmV2.extend(self, RmF z) @inline
{
	return { self.x, self.y, z };
}

fn Vec2 RmV2.lane(self, int i) @inline
{
	return { self.x[i], self.y[i] };
}

fn void RmV2.set_lane(&self, int i, Vec2 v) @inline
{
	self.x[i] = v.x;
	self.y[i] = v.y;
}

fn RmV2 select2(RmM m, RmV2 a, RmV2 b) @inline
{
	return { select(m, a.x, b.x), select(m, a.y, b.y) };
}

fn RmV2 gather2(Vec2[] table, RmI id) @inline
{
	RmV2 r;
	for (int i = 0; i < LANES; i++)
	{
		Vec2 c = table[id[i]];
		r.x[i] = c.x;
		r.y[i] = c.y;
	}
	return r;
}

// ================================
// Public: RmV3 Definitions
// ================================

struct RmV3
{
	RmF x, y, z;
}

fn RmV3 v3(float x, float y, float z) @inline
{
	return { (RmF)x, (RmF)y, (RmF)z };
}

fn RmV3 v3_splat(RmF s) @inline
{
	return { s, s, s };
}

fn RmV3 v3_from(Vec3 v) @inline
{
	return { (RmF)v.x, (RmF)v.y, (RmF)v.z };
}

// ================================
// Public: RmV3 Operators
// ================================

fn RmV3 RmV3.add(self, RmV3 o) @operator(+) @inline
{
	return { self.x + o.x, self.y + o.y, self.z + o.z };
}

fn RmV3 RmV3.sub(self, RmV3 o) @operator(-) @inline
{
	return { self.x - o.x, self.y - o.y, self.z - o.z };
}

fn RmV3 RmV3.mul(self, RmV3 o) @operator(*) @inline
{
	return { self.x * o.x, self.y * o.y, self.z * o.z };
}

fn RmV3 RmV3.div(self, RmV3 o) @operator(/) @inline
{
	return { self.x / o.x, self.y / o.y, self.z / o.z };
}

fn RmV3 RmV3.neg(self) @operator(-) @inline
{
	return { -self.x, -self.y, -self.z };
}

fn RmV3 RmV3.add_each(self, RmF s) @operator_s(+) @inline
{
	return { self.x + s, self.y + s, self.z + s };
}

fn RmV3 RmV3.sub_each(self, RmF s) @operator(-) @inline
{
	return { self.x - s, self.y - s, self.z - s };
}

fn RmV3 RmV3.sub_each_r(self, RmF s) @operator_r(-) @inline
{
	return { s - self.x, s - self.y, s - self.z };
}

fn RmV3 RmV3.scale(self, RmF s) @operator_s(*) @inline
{
	return { self.x * s, self.y * s, self.z * s };
}

fn RmV3 RmV3.div_each(self, RmF s) @operator(/) @inline
{
	return { self.x / s, self.y / s, self.z / s };
}

fn RmV3 RmV3.div_each_r(self, RmF s) @operator_r(/) @inline
{
	return { s / self.x, s / self.y, s / self.z };
}

fn RmV3 RmV3.add_eachf(self, float s) @operator_s(+) @inline
{
	return self.add_each((RmF)s);
}

fn RmV3 RmV3.sub_eachf(self, float s) @operator(-) @inline
{
	return self.sub_each((RmF)s);
}

fn RmV3 RmV3.sub_each_rf(self, float s) @operator_r(-) @inline
{
	return self.sub_each_r((RmF)s);
}

fn RmV3 RmV3.scalef(self, float s) @operator_s(*) @inline
{
	return self.scale((RmF)s);
}

fn RmV3 RmV3.div_eachf(self, float s) @operator(/) @inline
{
	return self.div_each((RmF)s);
}

fn RmV3 RmV3.div_each_rf(self, float s) @operator_r(/) @inline
{
	return self.div_each_r((RmF)s);
}

// ================================
// Public: RmV3 Comparison
// ================================

fn RmM RmV3.comp_eq(self, RmV3 o) @inline
{
	return self.x.comp_eq(o.x) & self.y.comp_eq(o.y) & self.z.comp_eq(o.z);
}

fn RmM RmV3.comp_ne(self, RmV3 o) @inline
{
	return self.x.comp_ne(o.x) | self.y.comp_ne(o.y) | self.z.comp_ne(o.z);
}

fn RmM RmV3.is_normalized(self, float eps = 1e-4f) @inline
{
	return comp_lt(math::abs(self.length_sq() - 1.0f), eps);
}

// ================================
// Public: RmV3 Reductions
// ================================

fn RmF RmV3.sum(self) @inline
{
	return self.x + self.y + self.z;
}

fn RmF RmV3.product(self) @inline
{
	return self.x * self.y * self.z;
}

fn RmF RmV3.max(self) @inline
{
	return math::max(self.x, math::max(self.y, self.z));
}

fn RmF RmV3.min(self) @inline
{
	return math::min(self.x, math::min(self.y, self.z));
}

// ================================
// Public: RmV3 Geometry
// ================================

fn RmF RmV3.dot(self, RmV3 o) @inline
{
	return self.x * o.x + self.y * o.y + self.z * o.z;
}

fn RmV3 RmV3.cross(self, RmV3 o) @inline
{
	return {
		self.y * o.z - self.z * o.y,
		self.z * o.x - self.x * o.z,
		self.x * o.y - self.y * o.x,
	};
}

fn RmF RmV3.length_sq(self) @inline
{
	return self.dot(self);
}

fn RmF RmV3.length(self) @inline
{
	return math::sqrt(self.length_sq());
}

fn RmF RmV3.distance_sq(self, RmV3 o) @inline
{
	return (self - o).length_sq();
}

fn RmF RmV3.distance(self, RmV3 o) @inline
{
	return (self - o).length();
}

fn RmV3 RmV3.normalize(self) @inline
{
	RmF inv = (RmF)1.0f / math::sqrt(math::max(self.length_sq(), (RmF)1e-20f));
	return self * inv;
}

fn RmV3 RmV3.perpendicular(self) @inline
{
	RmV3 a = self.abs();

	RmM mx = comp_le(a.x, a.y) & comp_le(a.x, a.z);
	RmM my = mask_not(mx) & comp_le(a.y, a.z);
	RmM mz = mask_not(mx) & mask_not(my);

	RmF  one  = (RmF)1.0f;
	RmF  zero = (RmF)0.0f;
	RmV3 axis = { select(mx, one, zero), select(my, one, zero), select(mz, one, zero) };

	return self.cross(axis);
}

fn RmV3 RmV3.reflect(self, RmV3 n) @inline
{
	return self - n * (self.dot(n) * 2.0f);
}

macro RmV3 RmV3.refract(self, RmV3 n, eta)
{
	RmF  e   = (RmF)eta;
	RmF  ndi = n.dot(self);
	RmF  k   = (RmF)1.0f - e * e * ((RmF)1.0f - ndi * ndi);
	RmV3 r   = self * e - n * (e * ndi + math::sqrt(math::max(k, (RmF)0.0f)));

	return select3(comp_ge(k, 0.0f), r, v3(0.0f, 0.0f, 0.0f));
}

fn RmV3 RmV3.rejection(self, RmV3 o) @inline
{
	return self - o * (self.dot(o) / o.dot(o));
}

fn RmV3 RmV3.barycentric(self, RmV3 a, RmV3 b, RmV3 c) @inline
{
	RmV3 v0 = b - a;
	RmV3 v1 = c - a;
	RmV3 v2 = self - a;

	RmF d00 = v0.dot(v0);
	RmF d01 = v0.dot(v1);
	RmF d11 = v1.dot(v1);
	RmF d20 = v2.dot(v0);
	RmF d21 = v2.dot(v1);

	RmF inv = (RmF)1.0f / (d00 * d11 - d01 * d01);
	RmF v   = (d11 * d20 - d01 * d21) * inv;
	RmF w   = (d00 * d21 - d01 * d20) * inv;

	return { (RmF)1.0f - v - w, v, w };
}

macro RmV3 RmV3.lerp(self, RmV3 b, t)
{
	return self + (b - self) * t;
}

macro RmV3 RmV3.clamp_length(self, lo, hi)
{
	RmF len = self.length();
	RmF tgt = math::clamp(len, (RmF)lo, (RmF)hi);

	return select3(comp_gt(len, 0.0f), self * (tgt / len), self);
}

macro RmV3 RmV3.move_towards(self, RmV3 target, max_distance)
{
	RmV3 d = target - self;
	RmF  l = d.length();
	RmF  m = (RmF)max_distance;

	return select3(comp_le(l, m) | comp_lt(l, 1e-8f), target, self + d * (m / l));
}

// ================================
// Public: RmV3 Per-Component
// ================================

fn RmV3 RmV3.abs(self) @inline
{
	return { math::abs(self.x), math::abs(self.y), math::abs(self.z) };
}

fn RmV3 RmV3.ceil(self) @inline
{
	return { self.x.ceil(), self.y.ceil(), self.z.ceil() };
}

fn RmV3 RmV3.floor(self) @inline
{
	return { self.x.floor(), self.y.floor(), self.z.floor() };
}

fn RmV3 RmV3.round(self) @inline
{
	return { self.x.round(), self.y.round(), self.z.round() };
}

fn RmV3 RmV3.trunc(self) @inline
{
	return { self.x.trunc(), self.y.trunc(), self.z.trunc() };
}

fn RmV3 RmV3.fract(self) @inline
{
	return { fract(self.x), fract(self.y), fract(self.z) };
}

fn RmV3 RmV3.saturate(self) @inline
{
	return { saturate(self.x), saturate(self.y), saturate(self.z) };
}

fn RmV3 RmV3.copysign(self, RmV3 sgn) @inline
{
	return { self.x.copysign(sgn.x), self.y.copysign(sgn.y), self.z.copysign(sgn.z) };
}

fn RmV3 RmV3.fma(self, RmV3 b, RmV3 c) @inline
{
	return { self.x.fma(b.x, c.x), self.y.fma(b.y, c.y), self.z.fma(b.z, c.z) };
}

macro RmV3 RmV3.pow(self, e)
{
	return { math::pow(self.x, (RmF)e), math::pow(self.y, (RmF)e), math::pow(self.z, (RmF)e) };
}

fn RmV3 RmV3.vmin(self, RmV3 o) @inline
{
	return { math::min(self.x, o.x), math::min(self.y, o.y), math::min(self.z, o.z) };
}

fn RmV3 RmV3.vmax(self, RmV3 o) @inline
{
	return { math::max(self.x, o.x), math::max(self.y, o.y), math::max(self.z, o.z) };
}

macro RmV3 RmV3.min_each(self, s)
{
	return { math::min(self.x, (RmF)s), math::min(self.y, (RmF)s), math::min(self.z, (RmF)s) };
}

macro RmV3 RmV3.max_each(self, s)
{
	return { math::max(self.x, (RmF)s), math::max(self.y, (RmF)s), math::max(self.z, (RmF)s) };
}

fn RmV3 RmV3.clamp(self, RmV3 lo, RmV3 hi) @inline
{
	return self.vmax(lo).vmin(hi);
}

macro RmV3 RmV3.clamp_each(self, lo, hi)
{
	return {
		math::clamp(self.x, (RmF)lo, (RmF)hi),
		math::clamp(self.y, (RmF)lo, (RmF)hi),
		math::clamp(self.z, (RmF)lo, (RmF)hi),
	};
}

// ================================
// Public: RmV3 Conversion
// ================================

fn RmV2 RmV3.xy(self) @inline => { self.x, self.y };
fn RmV2 RmV3.xz(self) @inline => { self.x, self.z };
fn RmV2 RmV3.yz(self) @inline => { self.y, self.z };

fn Vec3 RmV3.lane(self, int i) @inline
{
	return { self.x[i], self.y[i], self.z[i] };
}

fn void RmV3.set_lane(&self, int i, Vec3 v) @inline
{
	self.x[i] = v.x;
	self.y[i] = v.y;
	self.z[i] = v.z;
}

fn RmV3 RmV3.at(self, Vec3 offset) @inline
{
	return self - v3_from(offset);
}

fn RmV3 select3(RmM m, RmV3 a, RmV3 b) @inline
{
	return { select(m, a.x, b.x), select(m, a.y, b.y), select(m, a.z, b.z) };
}

fn RmV3 gather3(Vec3[] table, RmI id) @inline
{
	RmV3 r;
	for (int i = 0; i < LANES; i++)
	{
		Vec3 c = table[id[i]];
		r.x[i] = c.x;
		r.y[i] = c.y;
		r.z[i] = c.z;
	}
	return r;
}
