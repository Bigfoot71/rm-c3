module rm;

import std::math;

// ================================
// Public: SDF Definitions
// ================================

struct RmSd
{
	RmF d;   // signed distance
	RmI mat; // material id of the closest surface
}

// ================================
// Public: SDF Primitives
// ================================

fn RmSd sd_sphere(RmV3 p, float r) @inline
{
	return to_sd(p.length() - r);
}

fn RmSd sd_box(RmV3 p, Vec3 half) @inline
{
	return to_sd(box3(p.abs() - v3_from(half)));
}

fn RmSd sd_round_box(RmV3 p, Vec3 half, float r) @inline
{
	RmV3 q = p.abs() - (v3_from(half) - r);
	return to_sd(box3(q) - r);
}

fn RmSd sd_plane(RmV3 p, Vec3 n, float h) @inline
{
	return to_sd(p.dot(v3_from(n)) + h);
}

fn RmSd sd_plane_y(RmV3 p, float height) @inline
{
	return to_sd(p.y - height);
}

fn RmSd sd_capsule(RmV3 p, Vec3 a, Vec3 b, float r) @inline
{
	Vec3  ba     = b - a;
	float inv_bb = 1.0f / ba.dot(ba);

	RmV3 vba = v3_from(ba);
	RmV3 pa  = p.at(a);
	RmF  h   = saturate(pa.dot(vba) * inv_bb);

	return to_sd((pa - vba * h).length() - r);
}

fn RmSd sd_torus(RmV3 p, Vec2 t) @inline
{
	RmV2 q = { p.xz().length() - t.x, p.y };
	return to_sd(q.length() - t.y);
}

fn RmSd sd_capped_cylinder(RmV3 p, float r, float h) @inline
{
	return to_sd(rect2({ p.xz().length() - r, math::abs(p.y) - h }));
}

fn RmSd sd_cylinder(RmV3 p, Vec3 c) @inline
{
	return to_sd((p.xz() - v2(c.x, c.y)).length() - c.z);
}

fn RmSd sd_ellipsoid(RmV3 p, Vec3 r) @inline
{
	Vec3 inv  = 1.0f / r;
	Vec3 inv2 = inv * inv;

	RmF k0 = (p * v3_from(inv)).length();
	RmF k1 = (p * v3_from(inv2)).length();
	return to_sd(k0 * (k0 - 1.0f) / k1);
}

fn RmSd sd_cone(RmV3 p, Vec2 c, float h) @inline
{
	float qx     = h * c.x / c.y;
	float qy     = -h;
	float inv_qq = 1.0f / (qx * qx + qy * qy);
	float inv_qx = 1.0f / qx;
	float k      = sgnf(qy);

	RmV2 q = v2(qx, qy);
	RmV2 w = { p.xz().length(), p.y };

	RmF  t = saturate(w.dot(q) * inv_qq);
	RmV2 a = w - q * t;

	RmF  u = saturate(w.x * inv_qx);
	RmV2 b = { w.x - u * qx, w.y - qy };

	RmF d = math::min(a.length_sq(), b.length_sq());
	RmF s = math::max(w.cross(q) * k, (w.y - qy) * k);

	return to_sd(math::sqrt(d) * sgn(s));
}

fn RmSd sd_capped_cone(RmV3 p, float h, float r1, float r2) @inline
{
	Vec2  k1    = { r2, h };
	Vec2  k2    = { r2 - r1, 2.0f * h };
	float inv_k = 1.0f / k2.dot(k2);

	RmV2 q = { p.xz().length(), p.y };

	RmV2 ca = {
		q.x - math::min(q.x, select(comp_lt(q.y, 0.0f), (RmF)r1, (RmF)r2)),
		math::abs(q.y) - h,
	};

	RmF  t  = saturate((v2_from(k1) - q).dot(v2_from(k2)) * inv_k);
	RmV2 cb = q - v2_from(k1) + v2_from(k2) * t;

	RmM inside = comp_lt(cb.x, 0.0f) & comp_lt(ca.y, 0.0f);
	RmF s      = select(inside, (RmF)(-1.0f), (RmF)1.0f);

	return to_sd(s * math::sqrt(math::min(ca.length_sq(), cb.length_sq())));
}

fn RmSd sd_vertical_capsule(RmV3 p, float h, float r) @inline
{
	RmF  y = p.y - math::clamp(p.y, (RmF)0.0f, (RmF)h);
	RmV3 q = { p.x, y, p.z };
	return to_sd(q.length() - r);
}

fn RmSd sd_rounded_cylinder(RmV3 p, float ra, float rb, float h) @inline
{
	float a = ra - rb;
	float b = h - rb;
	return to_sd(rect2({ p.xz().length() - a, math::abs(p.y) - b }) - rb);
}

fn RmSd sd_octahedron(RmV3 p, float s) @inline
{
	RmV3 a = p.abs();
	RmF  m = a.sum() - s;

	RmM c1    = comp_lt(a.x * 3.0f, m);
	RmM c2    = comp_lt(a.y * 3.0f, m);
	RmM found = c1 | c2 | comp_lt(a.z * 3.0f, m);

	if (!found.or()) return to_sd(m * 0.57735027f);

	RmF qx = select(c1, a.x, select(c2, a.y, a.z));
	RmF qy = select(c1, a.y, select(c2, a.z, a.x));
	RmF qz = select(c1, a.z, select(c2, a.x, a.y));

	RmF  k = math::clamp((qz - qy + s) * 0.5f, (RmF)0.0f, (RmF)s);
	RmV3 d = { qx, qy - s + k, qz - k };

	return to_sd(select(found, d.length(), m * 0.57735027f));
}

fn RmSd sd_octahedron_bound(RmV3 p, float s) @inline
{
	return to_sd((p.abs().sum() - s) * 0.57735027f);
}

fn RmSd sd_pyramid(RmV3 p, float h) @inline
{
	float m2     = h * h + 0.25f;
	float inv_m2 = 1.0f / m2;
	float inv_t  = 1.0f / (m2 + 0.25f);

	RmV2 ab = p.xz().abs();
	RmF  px = ab.max() - 0.5f;
	RmF  pz = ab.min() - 0.5f;

	RmV3 q = { pz, p.y * h - px * 0.5f, px * h + p.y * 0.5f };

	RmF s  = math::max(-q.x, (RmF)0.0f);
	RmF t  = saturate((q.y - pz * 0.5f) * inv_t);
	RmF a  = (q.x + s) * (q.x + s) * m2 + q.y * q.y;
	RmF bx = q.x + t * 0.5f;
	RmF by = q.y - t * m2;
	RmF b  = bx * bx * m2 + by * by;

	RmM inside = comp_gt(math::min(q.y, -q.x * m2 - q.y * 0.5f), 0.0f);
	RmF d2     = select(inside, (RmF)0.0f, math::min(a, b));

	return to_sd(math::sqrt((d2 + q.z * q.z) * inv_m2) * sgn(math::max(q.z, -p.y)));
}

// h.x = radius, h.y = half height
fn RmSd sd_hex_prism(RmV3 p, Vec2 h) @inline
{
	const float KX = -0.8660254f;
	const float KY = 0.5f;
	const float KZ = 0.57735f;

	RmV3 a  = p.abs();
	RmV2 k  = v2(KX, KY);
	RmV2 xy = a.xy();

	// p.xy -= 2 * min(dot(k.xy, p.xy), 0) * k.xy
	RmF d = math::min(xy.dot(k), (RmF)0.0f) * 2.0f;
	xy = xy - k * d;

	float lim = KZ * h.x;
	RmF   cx  = math::clamp(xy.x, (RmF)(-lim), (RmF)lim);
	RmV2  e   = { xy.x - cx, xy.y - h.x };
	RmF   dx  = e.length() * sgn(e.y);
	RmF   dy  = a.z - h.y;

	return to_sd(rect2({ dx, dy }));
}

fn RmSd sd_tri_prism(RmV3 p, Vec2 h) @inline
{
	RmV3 q = p.abs();
	RmF  a = math::max(q.x * 0.866025f + p.y * 0.5f, -p.y) - h.x * 0.5f;
	return to_sd(math::max(q.z - h.y, a));
}

fn RmSd sd_box_frame(RmV3 p, Vec3 b, float e) @inline
{
	RmV3 pp = p.abs() - v3_from(b);
	RmV3 q  = (pp + e).abs() - e;

	RmF d0 = box3({ pp.x, q.y, q.z });
	RmF d1 = box3({ q.x, pp.y, q.z });
	RmF d2 = box3({ q.x, q.y, pp.z });

	return to_sd(math::min(math::min(d0, d1), d2));
}

fn RmSd sd_round_cone(RmV3 p, float r1, float r2, float h) @inline
{
	float b  = (r1 - r2) / h;
	float a  = math::sqrt(1.0f - b * b);
	float ah = a * h;

	RmV2 q = { p.xz().length(), p.y };
	RmF  k = q.y * a - q.x * b;

	RmM lo = comp_lt(k, 0.0f);
	RmM hi = comp_gt(k, ah) & mask_not(lo);

	RmV2 lq  = { q.x, select(hi, q.y - h, q.y) };
	RmF  cap = lq.length() - select(hi, (RmF)r2, (RmF)r1);
	RmF  mid = q.x * a + q.y * b - r1;

	return to_sd(select(lo | hi, cap, mid));
}

fn RmSd sd_capped_cylinder_ab(RmV3 p, Vec3 a, Vec3 b, float r) @inline
{
	Vec3  ba      = b - a;
	float baba    = ba.dot(ba);
	float inv     = 1.0f / baba;
	float half_bb = baba * 0.5f;

	RmV3 vba  = v3_from(ba);
	RmV3 pa   = p.at(a);
	RmF  paba = pa.dot(vba);
	RmV3 v    = pa * baba - vba * paba;

	RmF x  = v.length() - r * baba;
	RmF y  = math::abs(paba - half_bb) - half_bb;
	RmF x2 = x * x;
	RmF y2 = y * y * baba;

	RmM inside = comp_lt(math::max(x, y), 0.0f);
	RmF d_out  = select(comp_gt(x, 0.0f), x2, (RmF)0.0f) + select(comp_gt(y, 0.0f), y2, (RmF)0.0f);
	RmF d      = select(inside, -math::min(x2, y2), d_out);

	return to_sd(sgn(d) * math::sqrt(math::abs(d)) * inv);
}

fn RmSd sd_capped_cone_ab(RmV3 p, Vec3 a, Vec3 b, float ra, float rb) @inline
{
	Vec3  ba       = b - a;
	float rba      = rb - ra;
	float baba     = ba.dot(ba);
	float inv_baba = 1.0f / baba;
	float inv_k    = 1.0f / (rba * rba + baba);

	RmV3 pa   = p.at(a);
	RmF  papa = pa.length_sq();
	RmF  paba = pa.dot(v3_from(ba)) * inv_baba;

	RmF x   = math::sqrt(math::max(papa - paba * paba * baba, (RmF)0.0f));
	RmF cax = math::max(x - select(comp_lt(paba, 0.5f), (RmF)ra, (RmF)rb), (RmF)0.0f);
	RmF cay = math::abs(paba - 0.5f) - 0.5f;
	RmF f   = saturate(((x - ra) * rba + paba * baba) * inv_k);
	RmF cbx = x - ra - f * rba;
	RmF cby = paba - f;

	RmM inside = comp_lt(cbx, 0.0f) & comp_lt(cay, 0.0f);
	RmF s      = select(inside, (RmF)(-1.0f), (RmF)1.0f);
	RmF d2     = math::min(cax * cax + cay * cay * baba, cbx * cbx + cby * cby * baba);

	return to_sd(s * math::sqrt(d2));
}

fn RmSd sd_round_cone_ab(RmV3 p, Vec3 a, Vec3 b, float r1, float r2) @inline
{
	Vec3  ba  = b - a;
	float l2  = ba.dot(ba);
	float rr  = r1 - r2;
	float a2  = l2 - rr * rr;
	float il2 = 1.0f / l2;
	float kk  = sgnf(rr) * rr * rr;

	RmV3 vba = v3_from(ba);
	RmV3 pa  = p.at(a);
	RmF  y   = pa.dot(vba);
	RmF  z   = y - l2;
	RmV3 v   = pa * l2 - vba * y;
	RmF  x2  = v.length_sq();
	RmF  y2  = y * y * l2;
	RmF  z2  = z * z * l2;
	RmF  k   = x2 * kk;

	RmM m1 = comp_gt(sgn(z) * z2 * a2, k);
	RmM m2 = comp_lt(sgn(y) * y2 * a2, k) & mask_not(m1);

	RmF res = (math::sqrt(x2 * a2 * il2) + y * rr) * il2 - r1;
	if (m1.or()) res = select(m1, math::sqrt(x2 + z2) * il2 - r2, res);
	if (m2.or()) res = select(m2, math::sqrt(x2 + y2) * il2 - r1, res);

	return to_sd(res);
}

fn RmSd sd_cone_infinite(RmV3 p, Vec2 c) @inline
{
	RmV2 q  = { p.xz().length(), -p.y };
	RmV2 vc = v2_from(c);

	RmF t   = math::max(q.dot(vc), (RmF)0.0f);
	RmF d   = (q - vc * t).length();
	RmM neg = comp_lt(q.cross(vc), 0.0f);

	return to_sd(select(neg, -d, d));
}

fn RmSd sd_capped_torus(RmV3 p, Vec2 sc, float ra, float rb) @inline
{
	RmV3 q  = { math::abs(p.x), p.y, p.z };
	RmV2 xy = q.xy();

	RmM m  = comp_gt(q.x * sc.y, q.y * sc.x);
	RmF k  = select(m, xy.dot(v2_from(sc)), xy.length());
	RmF d2 = q.length_sq() + ra * ra - k * (2.0f * ra);

	return to_sd(math::sqrt(max0(d2)) - rb);
}

fn RmSd sd_link(RmV3 p, float le, float r1, float r2) @inline
{
	RmV2 q = { p.x, max0(math::abs(p.y) - le) };
	RmV2 w = { q.length() - r1, p.z };
	return to_sd(w.length() - r2);
}

fn RmSd sd_cut_sphere(RmV3 p, float r, float h) @inline
{
	float w  = math::sqrt(r * r - h * h);
	float w2 = w * w;

	RmV2 q = { p.xz().length(), p.y };

	RmF s = math::max(q.x * q.x * (h - r) + ((RmF)(h + r) - q.y * 2.0f) * w2, q.x * h - q.y * w);

	RmM near = comp_lt(s, 0.0f);
	RmM side = comp_lt(q.x, w);
	RmM far  = mask_not(near) & mask_not(side);

	RmF l = select2(far, q - v2(w, h), q).length();

	return to_sd(select(near, l - r, select(side, (RmF)h - q.y, l)));
}

fn RmSd sd_cut_hollow_sphere(RmV3 p, float r, float h, float t) @inline
{
	float w = math::sqrt(r * r - h * h);

	RmV2 q = { p.xz().length(), p.y };

	RmM m = comp_lt(q.x * h, q.y * w);
	RmF l = select2(m, q - v2(w, h), q).length();

	return to_sd(select(m, l, math::abs(l - r)) - t);
}

fn RmSd sd_solid_angle(RmV3 p, Vec2 c, float ra) @inline
{
	RmV2 q  = { p.xz().length(), p.y };
	RmV2 vc = v2_from(c);

	RmF l = q.length() - ra;
	RmF t = math::clamp(q.dot(vc), (RmF)0.0f, (RmF)ra);
	RmF m = (q - vc * t).length();

	return to_sd(math::max(l, m * sgn(q.cross(vc))));
}

fn RmSd sd_vesica_segment(RmV3 p, Vec3 a, Vec3 b, float w) @inline
{
	Vec3  c  = (a + b) * 0.5f;
	Vec3  ba = b - a;
	float l  = ba.length();
	Vec3  v  = ba * (1.0f / l);
	float r  = l * 0.5f;
	float d  = 0.5f * (r * r - w * w) / w;

	RmV3 vv = v3_from(v);
	RmV3 pc = p.at(c);
	RmF  y  = pc.dot(vv);
	RmV3 e  = pc - vv * y;
	RmV2 q  = { e.length(), math::abs(y) };

	RmM  m  = comp_lt(q.x * r, (q.y - r) * d);
	RmV2 hv = { select(m, (RmF)0.0f, (RmF)(-d)), select(m, (RmF)r, (RmF)0.0f) };
	RmF  hz = select(m, (RmF)0.0f, (RmF)(d + w));

	return to_sd((q - hv).length() - hz);
}

fn RmSd sd_rhombus(RmV3 p, float la, float lb, float h, float ra) @inline
{
	float inv = 1.0f / (la * la + lb * lb);

	RmV3 a = p.abs();
	RmF  f = saturate((a.x * la - a.z * lb + lb * lb) * inv);

	RmV2 w = { a.x - f * la, a.z - ((RmF)1.0f - f) * lb };

	RmF qx = w.length() * sgn(w.x) - ra;
	RmF qy = a.y - h;

	return to_sd(rect2({ qx, qy }));
}

fn RmSd sd_death_star(RmV3 p, float ra, float rb, float d) @inline
{
	float a = (ra * ra - rb * rb + d * d) / (2.0f * d);
	float b = math::sqrt(math::max(ra * ra - a * a, 0.0f));

	RmV2 q = { p.x, p.yz().length() };

	RmM m = comp_gt(q.x * b - q.y * a, math::max((RmF)b - q.y, (RmF)0.0f) * d);

	RmF res = (RmF)0.0f;
	if (m.or()) res = (q - v2(a, b)).length();
	if (!m.and())
	{
		RmF l1 = q.length() - ra;
		RmF l2 = rb - (q - v2(d, 0.0f)).length();
		res = select(m, res, math::max(l1, l2));
	}
	return to_sd(res);
}

fn RmSd ud_triangle(RmV3 p, Vec3 a, Vec3 b, Vec3 c) @inline
{
	Vec3 ba  = b - a;
	Vec3 cb  = c - b;
	Vec3 ac  = a - c;
	Vec3 nor = ba.cross(ac);

	Vec3 u1 = ba.cross(nor);
	Vec3 u2 = cb.cross(nor);
	Vec3 u3 = ac.cross(nor);

	float inv_ba  = 1.0f / ba.dot(ba);
	float inv_cb  = 1.0f / cb.dot(cb);
	float inv_ac  = 1.0f / ac.dot(ac);
	float inv_nor = 1.0f / nor.dot(nor);

	RmV3 pa = p.at(a);
	RmV3 pb = p.at(b);
	RmV3 pc = p.at(c);

	RmF s = sgn(pa.dot(v3_from(u1))) + sgn(pb.dot(v3_from(u2))) + sgn(pc.dot(v3_from(u3)));
	RmM edge = comp_lt(s, 2.0f);

	RmF nd = pa.dot(v3_from(nor));
	RmF d2 = nd * nd * inv_nor;

	if (edge.or())
	{
		RmF e = math::min(math::min(edge_d2(pa, ba, inv_ba), edge_d2(pb, cb, inv_cb)), edge_d2(pc, ac, inv_ac));
		d2 = select(edge, e, d2);
	}
	return to_sd(math::sqrt(d2));
}

fn RmSd ud_quad(RmV3 p, Vec3 a, Vec3 b, Vec3 c, Vec3 d) @inline
{
	Vec3 ba  = b - a;
	Vec3 cb  = c - b;
	Vec3 dc  = d - c;
	Vec3 ad  = a - d;
	Vec3 nor = ba.cross(ad);

	Vec3 u1 = ba.cross(nor);
	Vec3 u2 = cb.cross(nor);
	Vec3 u3 = dc.cross(nor);
	Vec3 u4 = ad.cross(nor);

	float inv_ba  = 1.0f / ba.dot(ba);
	float inv_cb  = 1.0f / cb.dot(cb);
	float inv_dc  = 1.0f / dc.dot(dc);
	float inv_ad  = 1.0f / ad.dot(ad);
	float inv_nor = 1.0f / nor.dot(nor);

	RmV3 pa = p.at(a);
	RmV3 pb = p.at(b);
	RmV3 pc = p.at(c);
	RmV3 pd = p.at(d);

	RmF s = sgn(pa.dot(v3_from(u1))) + sgn(pb.dot(v3_from(u2))) + sgn(pc.dot(v3_from(u3))) + sgn(pd.dot(v3_from(u4)));
	RmM edge = comp_lt(s, 3.0f);

	RmF nd = pa.dot(v3_from(nor));
	RmF d2 = nd * nd * inv_nor;

	if (edge.or())
	{
		RmF e0 = math::min(edge_d2(pa, ba, inv_ba), edge_d2(pb, cb, inv_cb));
		RmF e1 = math::min(edge_d2(pc, dc, inv_dc), edge_d2(pd, ad, inv_ad));
		d2 = select(edge, math::min(e0, e1), d2);
	}
	return to_sd(math::sqrt(d2));
}

// ================================
// Public: SDF Operators
// ================================

fn RmSd RmSd.join(self, RmSd o) @operator(|) @inline
{
	RmM first = comp_le(self.d, o.d);
	return { math::min(self.d, o.d), select(first, self.mat, o.mat) };
}

fn RmSd RmSd.clip(self, RmSd o) @operator(&) @inline
{
	RmM first = comp_ge(self.d, o.d);
	return { select(first, self.d, o.d), select(first, self.mat, o.mat) };
}

// The carved surface takes the material of the cutter
fn RmSd RmSd.cut(self, RmSd o) @operator(-) @inline
{
	RmM cutter = comp_gt(-o.d, self.d);
	return { select(cutter, -o.d, self.d), select(cutter, o.mat, self.mat) };
}

fn RmSd RmSd.smooth_join(self, RmSd o, float k) @inline
{
	RmF kk = (RmF)(k * 4.0f);
	RmF h  = math::max(kk - math::abs(self.d - o.d), (RmF)0.0f);
	RmM first = comp_le(self.d, o.d);
	return { math::min(self.d, o.d) - h * h * 0.25f / kk, select(first, self.mat, o.mat) };
}

fn RmSd RmSd.smooth_clip(self, RmSd o, float k) @inline
{
	RmSd a = { -self.d, self.mat };
	RmSd b = { -o.d, o.mat };
	RmSd r = a.smooth_join(b, k);
	return { -r.d, r.mat };
}

fn RmSd RmSd.smooth_cut(self, RmSd o, float k) @inline
{
	RmSd a = { -self.d, self.mat };
	RmSd r = a.smooth_join(o, k);
	return { -r.d, r.mat };
}

// ================================
// Public: SDF Modifiers
// ================================

// Tag a shape with a material id: p.sphere(0.5f).material(RED)
fn RmSd RmSd.material(self, int id) @inline => { self.d, (RmI)id };

// Inflate (r > 0) or deflate (r < 0)
fn RmSd RmSd.grow(self, float r) @inline => { self.d - (RmF)r, self.mat };

// Hollow shell of the given thickness
fn RmSd RmSd.shell(self, float thickness) @inline => { math::abs(self.d) - (RmF)thickness, self.mat };

// ================================
// Local: Helper Functions
// ================================

fn RmSd to_sd(RmF d) @local @inline
{
	return { d, (RmI)0 };
}

fn RmF sgn(RmF x) @local @inline
{
	RmF one = (RmF)1.0f;
	return select(comp_gt(x, 0.0f), one, select(comp_lt(x, 0.0f), -one, (RmF)0.0f));
}

fn float sgnf(float x) @local @inline
{
	return x > 0.0f ? 1.0f : (x < 0.0f ? -1.0f : 0.0f);
}

fn RmF max0(RmF x) @local @inline
{
	return math::max(x, (RmF)0.0f);
}

// 2D box distance: d = (|p| - half extents)
fn RmF rect2(RmV2 d) @local @inline
{
	return math::min(d.max(), (RmF)0.0f) + d.max_each(0.0f).length();
}

// 3D box distance: q = (|p| - half extents)
fn RmF box3(RmV3 q) @local @inline
{
	return math::min(q.max(), (RmF)0.0f) + q.max_each(0.0f).length();
}

fn RmF edge_d2(RmV3 pe, Vec3 e, float inv_ee) @local @inline
{
	RmV3 ve = v3_from(e);
	RmF  t  = saturate(pe.dot(ve) * inv_ee);
	return (ve * t - pe).length_sq();
}
