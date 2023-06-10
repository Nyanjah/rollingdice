#[derive(Copy, Clone, Debug)]
pub struct Quaternion {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub w: f64,
}

impl Default for Quaternion {
    fn default() -> Quaternion {
        Quaternion {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            w: 1.0,
        }
    }
}

impl Quaternion {

    pub fn from_axis_angle(mut axis: Vector3, angle: f64) -> Quaternion {
        let half_angle = angle / 2.0;
        
        axis = half_angle.sin() * axis.unit();
        
        Quaternion {
            x: axis.x,
            y: axis.y,
            z: axis.z,
            w: half_angle.cos(),
        }
    }

    pub fn to_axis_angle(mut self) -> (Vector3, f64) {
        self = self.unit();
        let angle = 2.0 * self.w.acos();
        let sin = (1.0 - self.w * self.w).sqrt();
        let axis = Vector3 {
            x: self.x / sin,
            y: self.y / sin,
            z: self.z / sin,
        };
        (axis.unit(), angle)
    }

    pub fn magnitude(self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z + self.w * self.w).sqrt()
    }

    pub fn unit(self) -> Quaternion {
        let magnitude = self.magnitude();
        Quaternion {
            x: self.x / magnitude,
            y: self.y / magnitude,
            z: self.z / magnitude,
            w: self.w / magnitude,
        }
    }

    pub fn inverse(self) -> Quaternion {
        Quaternion {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: self.w,
        }
    }

    pub fn basis_vectors(self) -> (Vector3, Vector3, Vector3) {
        let x = self.vector_to_world_space(Vector3::X_AXIS);
        let y = self.vector_to_world_space(Vector3::Y_AXIS);
        let z = self.vector_to_world_space(Vector3::Z_AXIS);
        (x, y, z)
    }

    pub fn vector_to_world_space(self, vector: Vector3) -> Vector3 {
        let q = self * Quaternion {
            x: vector.x,
            y: vector.y,
            z: vector.z,
            w: 0.0,
        } * self.inverse();
        Vector3 {
            x: q.x,
            y: q.y,
            z: q.z,
        }
    }

    pub fn vector_to_local_space(self, vector: Vector3) -> Vector3 {
        let q = self.inverse() * Quaternion {
            x: vector.x,
            y: vector.y,
            z: vector.z,
            w: 0.0,
        } * self;
        Vector3 {
            x: q.x,
            y: q.y,
            z: q.z,
        }
    }
}

impl std::ops::Mul for Quaternion {
    type Output = Quaternion;

    fn mul(self, rhs: Quaternion) -> Quaternion {
        Quaternion {
            x: self.w * rhs.x + self.x * rhs.w + self.y * rhs.z - self.z * rhs.y,
            y: self.w * rhs.y - self.x * rhs.z + self.y * rhs.w + self.z * rhs.x,
            z: self.w * rhs.z + self.x * rhs.y - self.y * rhs.x + self.z * rhs.w,
            w: self.w * rhs.w - self.x * rhs.x - self.y * rhs.y - self.z * rhs.z,
        }
    }
}

impl std::ops::MulAssign for Quaternion {
    fn mul_assign(&mut self, rhs: Quaternion) {
        *self = *self * rhs;
    }
}

#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Vector2 {
    pub x: f64,
    pub y: f64,
}

impl Vector2 {

    pub const ZERO: Vector2 = Vector2 { x: 0.0, y: 0.0 };

    pub const X_AXIS: Vector2 = Vector2 { x: 1.0, y: 0.0 };
    pub const Y_AXIS: Vector2 = Vector2 { x: 0.0, y: 1.0 };

    pub const fn new(x: f64, y: f64) -> Vector2 {
        Vector2 { x, y }
    }

    pub fn dot(self, rhs: Vector2) -> f64 {
        self.x * rhs.x + self.y * rhs.y
    }

    pub fn magnitude(self) -> f64 {
        self.dot(self).sqrt()
    }

    pub fn unit(self) -> Vector2 {
        let magnitude = self.magnitude();
        Vector2 {
            x: self.x / magnitude,
            y: self.y / magnitude,
        }
    }

    pub fn inverse(self) -> Vector2 {
        Vector2 {
            x: -self.x,
            y: -self.y
        }
    }

    pub fn lerp(self, rhs: Vector2, t: f64) -> Vector2 {
        Vector2 {
            x: self.x * (1.0 - t) + rhs.x * t,
            y: self.y * (1.0 - t) + rhs.y * t,
        }
    }
}

impl std::ops::Mul<f64> for Vector2 {
    type Output = Vector2;

    fn mul(self, rhs: f64) -> Vector2 {
        Vector2 {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

impl std::ops::Mul<Vector2> for f64 {
    type Output = Vector2;

    fn mul(self, rhs: Vector2) -> Vector2 {
        Vector2 {
            x: self * rhs.x,
            y: self * rhs.y
        }
    }
}

impl std::ops::Add for Vector2 {
    type Output = Vector2;

    fn add(self, rhs: Vector2) -> Vector2 {
        Vector2 {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl std::ops::AddAssign for Vector2 {
    fn add_assign(&mut self, rhs: Vector2) {
        *self = *self + rhs;
    }
}

impl std::ops::Sub for Vector2 {
    type Output = Vector2;

    fn sub(self, rhs: Vector2) -> Vector2 {
        Vector2 {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl std::ops::SubAssign for Vector2 {
    fn sub_assign(&mut self, rhs: Vector2) {
        *self = *self - rhs;
    }
}

impl std::ops::Mul for Vector2 {
    type Output = Vector2;

    fn mul(self, rhs: Vector2) -> Vector2 {
        Vector2 {
            x: self.x * rhs.x,
            y: self.y * rhs.y,
        }
    }
}

impl std::ops::Neg for Vector2 {
    type Output = Vector2;

    fn neg(self) -> Self::Output {
        Vector2 {
            x: -self.x,
            y: -self.y,
        }
    }
}

#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Vector3 {
    pub x: f64,
    pub y: f64,
    pub z: f64
}

impl Vector3 {

    pub const ZERO: Vector3 = Vector3 { x: 0.0, y: 0.0, z: 0.0 };

    pub const X_AXIS: Vector3 = Vector3 { x: 1.0, y: 0.0, z: 0.0 };
    pub const Y_AXIS: Vector3 = Vector3 { x: 0.0, y: 1.0, z: 0.0 };
    pub const Z_AXIS: Vector3 = Vector3 { x: 0.0, y: 0.0, z: 1.0 };

    pub const fn new(x: f64, y: f64, z: f64) -> Vector3 {
        Vector3 { x, y, z }
    }

    pub fn dot(self, rhs: Vector3) -> f64 {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }

    pub fn cross(self, rhs: Vector3) -> Vector3 {
        Vector3 {
            x: self.y * rhs.z - self.z * rhs.y,
            y: self.z * rhs.x - self.x * rhs.z,
            z: self.x * rhs.y - self.y * rhs.x,
        }
    }

    pub fn magnitude(self) -> f64 {
        self.dot(self).sqrt()
    }

    pub fn unit(self) -> Vector3 {
        let magnitude = self.magnitude();
        Vector3 {
            x: self.x / magnitude,
            y: self.y / magnitude,
            z: self.z / magnitude,
        }
    }

    pub fn inverse(self) -> Vector3 {
        Vector3 {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }

    pub fn lerp(self, rhs: Vector3, t: f64) -> Vector3 {
        Vector3 {
            x: self.x * (1.0 - t) + rhs.x * t,
            y: self.y * (1.0 - t) + rhs.y * t,
            z: self.z * (1.0 - t) + rhs.z * t,
        }
    }
}

impl std::ops::Mul<f64> for Vector3 {
    type Output = Vector3;

    fn mul(self, rhs: f64) -> Vector3 {
        Vector3 {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl std::ops::Mul<Vector3> for f64 {
    type Output = Vector3;

    fn mul(self, rhs: Vector3) -> Vector3 {
        Vector3 {
            x: self * rhs.x,
            y: self * rhs.y,
            z: self * rhs.z,
        }
    }
}

impl std::ops::Add for Vector3 {
    type Output = Vector3;

    fn add(self, rhs: Vector3) -> Vector3 {
        Vector3 {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl std::ops::AddAssign for Vector3 {
    fn add_assign(&mut self, rhs: Vector3) {
        *self = *self + rhs;
    }
}

impl std::ops::Sub for Vector3 {
    type Output = Vector3;

    fn sub(self, rhs: Vector3) -> Vector3 {
        Vector3 {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl std::ops::SubAssign for Vector3 {
    fn sub_assign(&mut self, rhs: Vector3) {
        *self = *self - rhs;
    }
}

impl std::ops::Mul for Vector3 {
    type Output = Vector3;

    fn mul(self, rhs: Vector3) -> Vector3 {
        Vector3 {
            x: self.x * rhs.x,
            y: self.y * rhs.y,
            z: self.z * rhs.z,
        }
    }
}

impl std::ops::Neg for Vector3 {
    type Output = Vector3;

    fn neg(self) -> Self::Output {
        Vector3 {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct Transform {
    pub position: Vector3,
    pub rotation: Quaternion,
}

impl Transform {
    pub fn new(position: Vector3, rotation: Quaternion) -> Transform {
        Transform { position, rotation: rotation.unit() }
    }

    pub fn inverse(&self) -> Transform {
        Transform {
            position: self.rotation.vector_to_local_space(self.position.inverse()),
            rotation: self.rotation.inverse().unit(),
        }
    }

    pub fn point_to_world_space(&self, point: Vector3) -> Vector3 {
        self.rotation.vector_to_world_space(point) + self.position
    }

    pub fn point_to_local_space(&self, point: Vector3) -> Vector3 {
        self.rotation.vector_to_local_space(point - self.position)
    }
}

impl std::ops::Mul for Transform {
    type Output = Transform;

    fn mul(self, rhs: Transform) -> Transform {
        Transform {
            position: self.rotation.vector_to_world_space(rhs.position.inverse()),
            rotation: (self.rotation * rhs.rotation).unit(),
        }
    }
}

pub trait Transformable {
    fn transform(&self) -> &Transform;
    fn transform_mut(&mut self) -> &mut Transform;
}