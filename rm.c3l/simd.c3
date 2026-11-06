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

// Per-lane choice: m ? a : b (works for RmF, RmI, RmM)
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

macro fpow(x, p)
{
	return math::exp(math::ln(math::max(x, ($Typeof(x))1e-6f)) * p);
}

macro fsin(x)
{
	var t = x * ($Typeof(x))INV_TAU;
	t = t - math::floor(t + 0.5f);
	var y = t * ($Typeof(x))TAU; // [-pi, pi)
	var s = y * (1.2732395f - 0.4052847f * math::abs(y));
	return s * (0.225f * (math::abs(s) - 1.0f) + 1.0f);
}

macro fcos(x) => fsin(x + ($Typeof(x))HALF_PI);

macro saturate(x)
{
	return math::clamp(x, ($Typeof(x))0, ($Typeof(x))1);
}

macro fract(x)
{
	return x - math::floor(x);
}

macro mix(a, b, t)
{
	return a + (b - a) * t;
}

fn RmF ign(RmF x, RmF y) @inline
{
	const RmF MX = 0.06711056f;
	const RmF MY = 0.00583715f;
	const RmF MZ = 52.9829189f;

	return fract(MZ * fract(x * MX + y * MY));
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

fn RmV3 v3_from(Vec3 v) @inline
{
	return { (RmF)v.x, (RmF)v.y, (RmF)v.z };
}

fn RmV3 RmV3.add(self, RmV3 o) @operator(+) @inline
{
	return { self.x + o.x, self.y + o.y, self.z + o.z };
}

fn RmV3 RmV3.sub(self, RmV3 o) @operator(-) @inline
{
	return { self.x - o.x, self.y - o.y, self.z - o.z };
}

fn RmV3 RmV3.scale(self, RmF s) @operator(*) @inline
{
	return { self.x * s, self.y * s, self.z * s };
}

fn RmV3 RmV3.mul(self, RmV3 o) @inline
{
	return { self.x * o.x, self.y * o.y, self.z * o.z };
}

fn RmV3 RmV3.mulf(self, float s) @inline
{
	return self * (RmF)s;
}

fn RmF RmV3.dot(self, RmV3 o) @inline
{
	return self.x * o.x + self.y * o.y + self.z * o.z;
}

fn RmF RmV3.length(self) @inline
{
	return math::sqrt(self.dot(self));
}

fn RmV3 RmV3.abs(self) @inline
{
	return {
		math::abs(self.x),
		math::abs(self.y),
		math::abs(self.z)
	};
}

fn RmV3 RmV3.normalized(self) @inline
{
	RmF inv = (RmF)1.0f / math::sqrt(math::max(self.dot(self), (RmF)1e-20f));
	return self * inv;
}

fn RmV3 RmV3.saturate(self) @inline
{
	return {
		saturate(self.x),
		saturate(self.y),
		saturate(self.z),
	};
}

fn RmV3 RmV3.fract(self) @inline
{
	return {
		fract(self.x),
		fract(self.y),
		fract(self.z),
	};
}

fn RmV3 RmV3.mix(self, RmV3 b, RmF t) @inline
{
	return {
		mix(self.x, b.x, t),
		mix(self.y, b.y, t),
		mix(self.z, b.z, t),
	};
}

// Translate the sampling point (use before a primitive): p.at({1, 0, 0}).sphere(0.5f)
fn RmV3 RmV3.at(self, Vec3 offset) @inline
{
	return self - v3_from(offset);
}

// Per-lane choice for RmV3
fn RmV3 select3(RmM m, RmV3 a, RmV3 b) @inline
{
	return { select(m, a.x, b.x), select(m, a.y, b.y), select(m, a.z, b.z) };
}

// Gather a table entry per lane
fn RmV3 gather(Vec3[] table, RmI id) @inline
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
