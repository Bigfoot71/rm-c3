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
	return to_sd(p.length() - (RmF)r);
}

fn RmSd sd_box(RmV3 p, Vec3 half) @inline
{
	RmV3 q = p.abs() - v3_from(half);
	return to_sd(box3(q.x, q.y, q.z));
}

fn RmSd sd_round_box(RmV3 p, Vec3 half, float r) @inline
{
	RmV3 q = p.abs() - v3_from(half - (Vec3){ r, r, r });
	return to_sd(box3(q.x, q.y, q.z) - (RmF)r);
}

fn RmSd sd_plane(RmV3 p, Vec3 n, float h) @inline
{
	return to_sd(p.dot(v3_from(n)) + (RmF)h);
}

fn RmSd sd_plane_y(RmV3 p, float height) @inline
{
	return to_sd(p.y - (RmF)height);
}

fn RmSd sd_capsule(RmV3 p, Vec3 a, Vec3 b, float r) @inline
{
	Vec3  ba     = b - a;
	float inv_bb = 1.0f / ba.dot(ba);

	RmV3 pa = p - v3_from(a);
	RmF  h  = saturate(pa.dot(v3_from(ba)) * inv_bb);

	return to_sd((pa - v3_from(ba) * h).length() - (RmF)r);
}

fn RmSd sd_torus(RmV3 p, Vec2 t) @inline
{
	return to_sd(len2(len2(p.x, p.z) - (RmF)t.x, p.y) - (RmF)t.y);
}

fn RmSd sd_capped_cylinder(RmV3 p, float r, float h) @inline
{
	return to_sd(rect2(len2(p.x, p.z) - (RmF)r, math::abs(p.y) - (RmF)h));
}

fn RmSd sd_cylinder(RmV3 p, Vec3 c) @inline
{
	return to_sd(len2(p.x - (RmF)c.x, p.z - (RmF)c.y) - (RmF)c.z);
}

fn RmSd sd_ellipsoid(RmV3 p, Vec3 r) @inline
{
	Vec3 inv  = { 1.0f / r.x, 1.0f / r.y, 1.0f / r.z };
	Vec3 inv2 = { inv.x * inv.x, inv.y * inv.y, inv.z * inv.z };

	RmF k0 = p.mul(v3_from(inv)).length();
	RmF k1 = p.mul(v3_from(inv2)).length();
	return to_sd(k0 * (k0 - 1.0f) / k1);
}

fn RmSd sd_cone(RmV3 p, Vec2 c, float h) @inline
{
	float qx     = h * c.x / c.y;
	float qy     = -h;
	float inv_qq = 1.0f / (qx * qx + qy * qy);
	float inv_qx = 1.0f / qx;
	float k      = sgnf(qy);

	RmF wx = len2(p.x, p.z);
	RmF wy = p.y;

	RmF t  = saturate((wx * qx + wy * qy) * inv_qq);
	RmF ax = wx - t * qx;
	RmF ay = wy - t * qy;

	RmF u  = saturate(wx * inv_qx);
	RmF bx = wx - u * qx;
	RmF by = wy - qy;

	RmF d = math::min(ax * ax + ay * ay, bx * bx + by * by);
	RmF s = math::max((wx * qy - wy * qx) * k, (wy - qy) * k);

	return to_sd(math::sqrt(d) * sgn(s));
}

fn RmSd sd_capped_cone(RmV3 p, float h, float r1, float r2) @inline
{
	float k1x   = r2;
	float k1y   = h;
	float k2x   = r2 - r1;
	float k2y   = 2.0f * h;
	float inv_k = 1.0f / (k2x * k2x + k2y * k2y);

	RmF qx = len2(p.x, p.z);
	RmF qy = p.y;

	RmF cax = qx - math::min(qx, select(comp_lt(qy, 0.0f), (RmF)r1, (RmF)r2));
	RmF cay = math::abs(qy) - h;

	RmF t   = saturate(((k1x - qx) * k2x + (k1y - qy) * k2y) * inv_k);
	RmF cbx = qx - k1x + t * k2x;
	RmF cby = qy - k1y + t * k2y;

	RmM inside = comp_lt(cbx, 0.0f) & comp_lt(cay, 0.0f);
	RmF s      = select(inside, (RmF)(-1.0f), (RmF)1.0f);

	return to_sd(s * math::sqrt(math::min(cax * cax + cay * cay, cbx * cbx + cby * cby)));
}

fn RmSd sd_vertical_capsule(RmV3 p, float h, float r) @inline
{
	RmF  y = p.y - math::clamp(p.y, (RmF)0.0f, (RmF)h);
	RmV3 q = { p.x, y, p.z };
	return to_sd(q.length() - (RmF)r);
}

fn RmSd sd_rounded_cylinder(RmV3 p, float ra, float rb, float h) @inline
{
	float a = ra - rb;
	float b = h - rb;
	return to_sd(rect2(len2(p.x, p.z) - a, math::abs(p.y) - b) - (RmF)rb);
}

fn RmSd sd_octahedron(RmV3 p, float s) @inline
{
	RmV3 a = p.abs();
	RmF  m = a.x + a.y + a.z - s;

	RmM c1    = comp_lt(a.x * 3.0f, m);
	RmM c2    = comp_lt(a.y * 3.0f, m);
	RmM found = c1 | c2 | comp_lt(a.z * 3.0f, m);

	if (!found.or()) return to_sd(m * 0.57735027f);

	RmF qx = select(c1, a.x, select(c2, a.y, a.z));
	RmF qy = select(c1, a.y, select(c2, a.z, a.x));
	RmF qz = select(c1, a.z, select(c2, a.x, a.y));

	RmF k  = math::clamp((qz - qy + s) * 0.5f, (RmF)0.0f, (RmF)s);
	RmF dy = qy - s + k;
	RmF dz = qz - k;
	RmF r  = math::sqrt(qx * qx + dy * dy + dz * dz);

	return to_sd(select(found, r, m * 0.57735027f));
}

fn RmSd sd_octahedron_bound(RmV3 p, float s) @inline
{
	RmV3 a = p.abs();
	return to_sd((a.x + a.y + a.z - s) * 0.57735027f);
}

fn RmSd sd_pyramid(RmV3 p, float h) @inline
{
	float m2     = h * h + 0.25f;
	float inv_m2 = 1.0f / m2;
	float inv_t  = 1.0f / (m2 + 0.25f);

	RmF ax = math::abs(p.x);
	RmF az = math::abs(p.z);
	RmF px = math::max(ax, az) - 0.5f;
	RmF pz = math::min(ax, az) - 0.5f;

	RmF qx = pz;
	RmF qy = p.y * h - px * 0.5f;
	RmF qz = px * h + p.y * 0.5f;

	RmF s  = math::max(-qx, (RmF)0.0f);
	RmF t  = saturate((qy - pz * 0.5f) * inv_t);
	RmF a  = (qx + s) * (qx + s) * m2 + qy * qy;
	RmF bx = qx + t * 0.5f;
	RmF by = qy - t * m2;
	RmF b  = bx * bx * m2 + by * by;

	RmM inside = comp_gt(math::min(qy, -qx * m2 - qy * 0.5f), 0.0f);
	RmF d2     = select(inside, (RmF)0.0f, math::min(a, b));

	return to_sd(math::sqrt((d2 + qz * qz) * inv_m2) * sgn(math::max(qz, -p.y)));
}

// h.x = radius, h.y = half height
fn RmSd sd_hex_prism(RmV3 p, Vec2 h) @inline
{
	const float KX = -0.8660254f;
	const float KY = 0.5f;
	const float KZ = 0.57735f;

	RmF px = math::abs(p.x);
	RmF py = math::abs(p.y);
	RmF pz = math::abs(p.z);

	// p.xy -= 2 * min(dot(k.xy, p.xy), 0) * k.xy
	RmF d = math::min(px * KX + py * KY, (RmF)0.0f) * 2.0f;
	px = px - d * KX;
	py = py - d * KY;

	float lim = KZ * h.x;
	RmF   cx  = math::clamp(px, (RmF)(-lim), (RmF)lim);
	RmF   dx  = len2(px - cx, py - h.x) * sgn(py - h.x);
	RmF   dy  = pz - h.y;

	return to_sd(rect2(dx, dy));
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
	RmV3 q  = { math::abs(pp.x + e) - e, math::abs(pp.y + e) - e, math::abs(pp.z + e) - e };

	RmF d0 = box3(pp.x, q.y, q.z);
	RmF d1 = box3(q.x, pp.y, q.z);
	RmF d2 = box3(q.x, q.y, pp.z);

	return to_sd(math::min(math::min(d0, d1), d2));
}

fn RmSd sd_round_cone(RmV3 p, float r1, float r2, float h) @inline
{
	float b  = (r1 - r2) / h;
	float a  = math::sqrt(1.0f - b * b);
	float ah = a * h;

	RmF qx = len2(p.x, p.z);
	RmF qy = p.y;
	RmF k  = qy * a - qx * b;

	RmM lo = comp_lt(k, 0.0f);
	RmM hi = comp_gt(k, ah) & mask_not(lo);

	RmF l   = len2(qx, select(hi, qy - h, qy));
	RmF cap = l - select(hi, (RmF)r2, (RmF)r1);
	RmF mid = qx * a + qy * b - r1;

	return to_sd(select(lo | hi, cap, mid));
}

fn RmSd sd_capped_cylinder_ab(RmV3 p, Vec3 a, Vec3 b, float r) @inline
{
	Vec3  ba      = b - a;
	float baba    = ba.dot(ba);
	float inv     = 1.0f / baba;
	float half_bb = baba * 0.5f;

	RmV3 pa   = p - v3_from(a);
	RmF  paba = pa.dot(v3_from(ba));
	RmV3 v    = pa * (RmF)baba - v3_from(ba) * paba;

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

	RmV3 pa   = p - v3_from(a);
	RmF  papa = pa.dot(pa);
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

	RmV3 pa = p - v3_from(a);
	RmF  y  = pa.dot(v3_from(ba));
	RmF  z  = y - l2;
	RmV3 v  = pa * (RmF)l2 - v3_from(ba) * y;
	RmF  x2 = v.dot(v);
	RmF  y2 = y * y * l2;
	RmF  z2 = z * z * l2;
	RmF  k  = x2 * kk;

	RmM m1 = comp_gt(sgn(z) * z2 * a2, k);
	RmM m2 = comp_lt(sgn(y) * y2 * a2, k) & mask_not(m1);

	RmF res = (math::sqrt(x2 * a2 * il2) + y * rr) * il2 - r1;
	if (m1.or()) res = select(m1, math::sqrt(x2 + z2) * il2 - r2, res);
	if (m2.or()) res = select(m2, math::sqrt(x2 + y2) * il2 - r1, res);

	return to_sd(res);
}

fn RmSd sd_cone_infinite(RmV3 p, Vec2 c) @inline
{
	RmF qx = len2(p.x, p.z);
	RmF qy = -p.y;

	RmF t   = math::max(qx * c.x + qy * c.y, (RmF)0.0f);
	RmF d   = len2(qx - t * c.x, qy - t * c.y);
	RmM neg = comp_lt(qx * c.y - qy * c.x, 0.0f);

	return to_sd(select(neg, -d, d));
}

fn RmSd sd_capped_torus(RmV3 p, Vec2 sc, float ra, float rb) @inline
{
	RmF px = math::abs(p.x);
	RmM m  = comp_gt(px * sc.y, p.y * sc.x);
	RmF k  = select(m, px * sc.x + p.y * sc.y, len2(px, p.y));
	RmF d2 = px * px + p.y * p.y + p.z * p.z + ra * ra - k * (2.0f * ra);

	return to_sd(math::sqrt(max0(d2)) - (RmF)rb);
}

fn RmSd sd_link(RmV3 p, float le, float r1, float r2) @inline
{
	RmF qy = max0(math::abs(p.y) - le);
	RmF l  = len2(p.x, qy) - r1;
	return to_sd(len2(l, p.z) - (RmF)r2);
}

fn RmSd sd_cut_sphere(RmV3 p, float r, float h) @inline
{
	float w  = math::sqrt(r * r - h * h);
	float w2 = w * w;

	RmF qx = len2(p.x, p.z);
	RmF qy = p.y;

	RmF s = math::max(qx * qx * (h - r) + ((RmF)(h + r) - qy * 2.0f) * w2, qx * h - qy * w);

	RmM near = comp_lt(s, 0.0f);
	RmM side = comp_lt(qx, w);
	RmM far  = mask_not(near) & mask_not(side);

	RmF l = len2(select(far, qx - w, qx), select(far, qy - h, qy));

	return to_sd(select(near, l - r, select(side, (RmF)h - qy, l)));
}

fn RmSd sd_cut_hollow_sphere(RmV3 p, float r, float h, float t) @inline
{
	float w = math::sqrt(r * r - h * h);

	RmF qx = len2(p.x, p.z);
	RmF qy = p.y;

	RmM m = comp_lt(qx * h, qy * w);
	RmF l = len2(select(m, qx - w, qx), select(m, qy - h, qy));

	return to_sd(select(m, l, math::abs(l - r)) - (RmF)t);
}

fn RmSd sd_solid_angle(RmV3 p, Vec2 c, float ra) @inline
{
	RmF qx = len2(p.x, p.z);
	RmF qy = p.y;

	RmF l = len2(qx, qy) - ra;
	RmF t = math::clamp(qx * c.x + qy * c.y, (RmF)0.0f, (RmF)ra);
	RmF m = len2(qx - t * c.x, qy - t * c.y);

	return to_sd(math::max(l, m * sgn(qx * c.y - qy * c.x)));
}

fn RmSd sd_vesica_segment(RmV3 p, Vec3 a, Vec3 b, float w) @inline
{
	Vec3  c  = (a + b) * 0.5f;
	Vec3  ba = b - a;
	float l  = ba.length();
	Vec3  v  = ba * (1.0f / l);
	float r  = l * 0.5f;
	float d  = 0.5f * (r * r - w * w) / w;

	RmV3 pc = p - v3_from(c);
	RmF  y  = pc.dot(v3_from(v));
	RmV3 e  = pc - v3_from(v) * y;
	RmF  qx = e.length();
	RmF  qy = math::abs(y);

	RmM m  = comp_lt(qx * r, (qy - r) * d);
	RmF hx = select(m, (RmF)0.0f, (RmF)(-d));
	RmF hy = select(m, (RmF)r, (RmF)0.0f);
	RmF hz = select(m, (RmF)0.0f, (RmF)(d + w));

	return to_sd(len2(qx - hx, qy - hy) - hz);
}

fn RmSd sd_rhombus(RmV3 p, float la, float lb, float h, float ra) @inline
{
	float inv = 1.0f / (la * la + lb * lb);

	RmV3 a = p.abs();
	RmF  f = saturate((a.x * la - a.z * lb + lb * lb) * inv);

	RmF wx = a.x - f * la;
	RmF wz = a.z - ((RmF)1.0f - f) * lb;

	RmF qx = len2(wx, wz) * sgn(wx) - ra;
	RmF qy = a.y - h;

	return to_sd(rect2(qx, qy));
}

fn RmSd sd_death_star(RmV3 p, float ra, float rb, float d) @inline
{
	float a = (ra * ra - rb * rb + d * d) / (2.0f * d);
	float b = math::sqrt(math::max(ra * ra - a * a, 0.0f));

	RmF px = p.x;
	RmF py = len2(p.y, p.z);

	RmM m = comp_gt(px * b - py * a, math::max((RmF)b - py, (RmF)0.0f) * d);

	RmF res = (RmF)0.0f;
	if (m.or()) res = len2(px - a, py - b);
	if (!m.and())
	{
		RmF l1 = len2(px, py) - ra;
		RmF l2 = rb - len2(px - d, py);
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

	RmV3 pa = p - v3_from(a);
	RmV3 pb = p - v3_from(b);
	RmV3 pc = p - v3_from(c);

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

	RmV3 pa = p - v3_from(a);
	RmV3 pb = p - v3_from(b);
	RmV3 pc = p - v3_from(c);
	RmV3 pd = p - v3_from(d);

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

fn RmF len2(RmF x, RmF y) @local @inline
{
	return math::sqrt(x * x + y * y);
}

fn RmF rect2(RmF dx, RmF dy) @local @inline
{
	return math::min(math::max(dx, dy), (RmF)0.0f) + len2(max0(dx), max0(dy));
}

fn RmF box3(RmF a, RmF b, RmF c) @local @inline
{
	RmF ma = max0(a);
	RmF mb = max0(b);
	RmF mc = max0(c);
	return math::min(math::max(a, math::max(b, c)), (RmF)0.0f) + math::sqrt(ma * ma + mb * mb + mc * mc);
}

fn RmF edge_d2(RmV3 pe, Vec3 e, float inv_ee) @local @inline
{
	RmF  t = saturate(pe.dot(v3_from(e)) * inv_ee);
	RmV3 d = v3_from(e) * t - pe;
	return d.dot(d);
}
