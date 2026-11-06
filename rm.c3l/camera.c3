module rm;

import std::math;

// ================================
// Public: Camera Definitions
// ================================

const Vec3 WORLD_UP = { 0.0f, 1.0f, 0.0f };

struct RmCamera
{
	Vec3  position;
	Quat  rotation;
	float aspect;
	float fovy;
	float far;
}

// ================================
// Public: Camera Functions
// ================================

fn RmCamera camera_new(Vec3 position, Vec3 target, float aspect, float fovy = 60.0f, float far = 60.0f)
{
	RmCamera c = {
		.position = position,
		.rotation = { 0, 0, 0, 1 },
		.aspect   = aspect,
		.fovy     = fovy,
		.far      = far,
	};
	c.look_at(target);
	return c;
}

// ================================
// Public: Camera Axes
// ================================

fn Vec3 RmCamera.right(&self)   @inline => self.rotation.rotate_vec3({ 1.0f, 0.0f,  0.0f });
fn Vec3 RmCamera.up(&self)      @inline => self.rotation.rotate_vec3({ 0.0f, 1.0f,  0.0f });
fn Vec3 RmCamera.forward(&self) @inline => self.rotation.rotate_vec3({ 0.0f, 0.0f, -1.0f });
fn Vec3 RmCamera.back(&self)    @inline => self.rotation.rotate_vec3({ 0.0f, 0.0f,  1.0f });

fn float RmCamera.tan_half(&self) @inline => math::tan(self.fovy * (float)math::PI / 360.0f);

// ================================
// Public: Camera Movement
// ================================

fn void RmCamera.translate(&self, Vec3 delta) @inline
{
	self.position = self.position + delta;
}

fn void RmCamera.move(&self, Vec3 local)
{
	self.position = self.position
		+ self.right()   * local.x
		+ self.up()      * local.y
		+ self.forward() * local.z;
}

fn void RmCamera.walk(&self, float forward, float right, float up = 0.0f)
{
	Vec3 f = self.forward();
	f.y = 0.0f;
	if (f.length() < 1e-5f) f = self.up();
	f.y = 0.0f;
	f = f.normalize();

	Vec3 r = f.cross(WORLD_UP).normalize();

	self.position = self.position + f * forward + r * right + WORLD_UP * up;
}

fn void RmCamera.rotate_world(&self, Vec3 axis, float angle) @inline
{
	self.rotation = (quaternion::from_axis_angle{float}(axis, angle) * self.rotation).normalize();
}

fn void RmCamera.rotate_local(&self, Vec3 axis, float angle) @inline
{
	self.rotation = (self.rotation * quaternion::from_axis_angle{float}(axis, angle)).normalize();
}

fn void RmCamera.yaw(&self, float angle) @inline
{
	self.rotate_world(WORLD_UP, angle);
}

fn void RmCamera.pitch(&self, float angle) @inline
{
	self.rotate_local({ 1.0f, 0.0f, 0.0f }, angle);
}

fn void RmCamera.roll(&self, float angle) @inline
{
	self.rotate_local({ 0.0f, 0.0f, 1.0f }, angle);
}

fn void RmCamera.pitch_clamped(&self, float angle, float limit = 1.5533f) // ~89°
{
	float cur = math::asin(math::clamp(self.forward().y, -1.0f, 1.0f));
	float tgt = math::clamp(cur + angle, -limit, limit);
	self.pitch(tgt - cur);
}

fn void RmCamera.look_dir(&self, Vec3 dir, Vec3 up = WORLD_UP)
{
	Vec3 f = dir.normalize();
	Vec3 r = f.cross(up);

	if (r.length() < 1e-5f)
	{
		r = f.cross(math::abs(f.x) < 0.9f ? (Vec3) { 1.0f, 0.0f, 0.0f } : (Vec3) { 0.0f, 0.0f, 1.0f });
	}

	r = r.normalize();
	Vec3 u = r.cross(f);

	self.rotation = quat_from_basis(r, u, f * -1.0f);
}

fn void RmCamera.look_at(&self, Vec3 target, Vec3 up = WORLD_UP)
{
	self.look_dir(target - self.position, up);
}

fn void RmCamera.orbit(&self, Vec3 center, float yaw, float pitch)
{
	Vec3  off  = self.position - center;
	float dist = off.length();

	self.position = center;
	self.yaw(yaw);
	self.pitch_clamped(pitch);
	self.position = center - self.forward() * dist;
}

// ================================
// Local: Quaternion Helpers
// ================================

fn Quat quat_from_basis(Vec3 r, Vec3 u, Vec3 b) @local
{
	float m00 = r.x; float m10 = r.y; float m20 = r.z;
	float m01 = u.x; float m11 = u.y; float m21 = u.z;
	float m02 = b.x; float m12 = b.y; float m22 = b.z;

	float tr = m00 + m11 + m22;
	Quat  q;

	if (tr > 0.0f)
	{
		float s = math::sqrt(tr + 1.0f) * 2.0f;
		q = { (m21 - m12) / s, (m02 - m20) / s, (m10 - m01) / s, s * 0.25f };
	}
	else if (m00 > m11 && m00 > m22)
	{
		float s = math::sqrt(1.0f + m00 - m11 - m22) * 2.0f;
		q = { s * 0.25f, (m01 + m10) / s, (m02 + m20) / s, (m21 - m12) / s };
	}
	else if (m11 > m22)
	{
		float s = math::sqrt(1.0f + m11 - m00 - m22) * 2.0f;
		q = { (m01 + m10) / s, s * 0.25f, (m12 + m21) / s, (m02 - m20) / s };
	}
	else
	{
		float s = math::sqrt(1.0f + m22 - m00 - m11) * 2.0f;
		q = { (m02 + m20) / s, (m12 + m21) / s, s * 0.25f, (m10 - m01) / s };
	}

	return q.normalize();
}
