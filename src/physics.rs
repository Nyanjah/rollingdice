#![allow(dead_code)]

use std::iter::zip;
#[derive(Copy, Clone)]
pub struct Quaternion {
    pub angle: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}
pub struct Cube {
    side_length: f64,
    transform: Transform,
}

pub struct Transform {
    pub position: [f64; 3],     // Position vector
    pub quaternion: Quaternion, // Rotation quaternion to track orientation in world
}

pub trait Transformable {
    fn transform(&self) -> &Transform;
    fn transform_mut(&mut self) -> &mut Transform;

    fn translate(&mut self, vector: &[f64; 3]) {
        for (mut pos, translation) in zip(self.transform_mut().position, vector) {
            pos += translation;
        }
    }

    fn rotate(&mut self, quaternion: &Quaternion) {
        self.transform_mut().quaternion = self.transform_mut().quaternion * *quaternion;
    }
}

impl std::ops::Mul for Quaternion {
    type Output = Quaternion;
    fn mul(self, rhs: Self) -> Self::Output {
        let product = Quaternion {
            angle: self.angle * rhs.angle - (self.x * rhs.x + self.y * rhs.y + self.z * rhs.z),
            x: self.angle * rhs.x + self.x * rhs.angle + self.y * rhs.z - self.z * rhs.y,
            y: self.angle * rhs.y - self.x * rhs.z + self.y * rhs.angle + self.z * rhs.x,
            z: self.angle * rhs.z + self.x * rhs.y - self.y * rhs.x + self.z * rhs.angle,
        };
        return product;
    }
}

impl Quaternion {
    pub fn new(angle: f64, axis: &[f64; 3]) -> Self {
        let magnitude = (axis[0].powf(2.0) + axis[1].powf(2.0) + axis[2].powf(2.0)).sqrt();
        Quaternion {
            angle: (angle / 2.0).cos(),
            x: axis[0] * (angle / 2.0).sin() / magnitude,
            y: axis[1] * (angle / 2.0).sin() / magnitude,
            z: axis[2] * (angle / 2.0).sin() / magnitude,
        }
    }

    pub fn get_inverse(&self) -> Quaternion {
        Quaternion {
            angle: self.angle,
            x: -1.0 * self.x,
            y: -1.0 * self.y,
            z: -1.0 * self.z,
        }
    }
}

impl From<&[f64; 3]> for Quaternion {
    fn from(vector: &[f64; 3]) -> Self {
        return Quaternion {
            angle: 0.0,
            x: vector[0],
            y: vector[1],
            z: vector[2],
        };
    }
}

impl Transformable for Cube {
    fn transform(&self) -> &Transform {
        &self.transform
    }
    fn transform_mut(&mut self) -> &mut Transform {
        &mut self.transform
    }
}

impl Cube {
    pub fn new(side_length: f64, position: &[f64; 3], quaternion: Quaternion) -> Self {
        if side_length <= 0.0 {
            panic!("Attempted to instantiate cube of side length <= 0.");
        }
        Cube {
            side_length: side_length,
            transform: Transform {
                position: *position,
                quaternion: quaternion,
            },
        }
    }

    pub fn get_vertices(&self) -> Vec<[f64; 3]> {
        let mut vertices: Vec<[f64; 3]> = Vec::new();
        let length = self.side_length / 2.0;
        let (x_0, y_0, z_0) = (
            self.transform.position[0],
            self.transform.position[1],
            self.transform.position[2],
        );
        // Initial state -> Translation -> Rotation -> Resulting Vertices
        for x in [1.0, -1.0] {
            for y in [-1.0, 1.0] {
                for z in [-1.0, 1.0] {
                    // Initial vertices after the translation step
                    vertices.push([length * x + x_0, length * y + y_0, length * z + z_0]);
                }
            }
        }
        // Applying the cube's stored rotation quaternion
        let vertices = vertices
            .iter()
            .map(|vertex| {
                //p -> q * p * q^-1
                let rotation_results = (self.transform.quaternion * Quaternion::from(vertex))
                    * self.transform.quaternion.get_inverse();
                [rotation_results.x, rotation_results.y, rotation_results.z]
            })
            .collect();

        return vertices;
    }
}
