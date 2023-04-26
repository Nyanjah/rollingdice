#![allow(dead_code)]
use tobj::*;

#[derive(Copy, Clone)]
#[derive(Debug)]
pub struct Quaternion {
    pub angle: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}
pub struct Object {
    scale: f64,
    transform: Transform,
    surfaces: Vec<[[f64;3];3]>
}

pub struct Transform {
    pub position: [f64; 3],     // Position 
    pub quaternion: Quaternion, // Rotation
}

pub trait Transformable {
    fn transform(&self) -> &Transform;
    fn transform_mut(&mut self) -> &mut Transform;

    fn translate(&mut self, vector: &[f64; 3]) {
        for i in 0..3{
            // Increment the position by the vector
            self.transform_mut().position[i] += vector[i];
        }
    }

    fn rotate_globally(&mut self, angle:f64, axis:[f64;3]) {
        let mut quaternion = Quaternion::new(angle,&axis);
        quaternion.normalize();
        self.transform_mut().quaternion =  quaternion * self.transform_mut().quaternion.normalize() ;
    }

    fn rotate(&mut self, angle:f64, axis:[f64;3]) {
        let mut quaternion = Quaternion::new(angle,&axis);
        quaternion.normalize();
        self.transform_mut().quaternion = self.transform_mut().quaternion.normalize() * quaternion;
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
    pub fn normalize(&mut self) -> Quaternion{
        let magnitude = (self.angle.powf(2.0) + self.x.powf(2.0) + self.y.powf(2.0) + self.z.powf(2.0)).sqrt();
            self.angle =  self.angle / magnitude;
            self.x = self.x / magnitude;
            self.y = self.y / magnitude;
            self.z = self.z / magnitude;
        return *self
    }
    pub fn normalize_as_vector(&mut self) {
        let magnitude = (self.x.powf(2.0) + self.y.powf(2.0) + self.z.powf(2.0)).sqrt();
            self.x = self.x / magnitude;
            self.y = self.y / magnitude;
            self.z = self.z / magnitude;
    }

    pub fn get_inverse(&self) -> Quaternion {
        Quaternion {
            angle: self.angle,
            x: -1.0 * self.x,
            y: -1.0 * self.y,
            z: -1.0 * self.z,
        }
    }
    pub fn to_vector(&self) -> [f64;3]{
        let mut vector = [0.0;3];
        vector[0] = self.x;
        vector[1] = self.y;
        vector[2] = self.z;
        return vector
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

impl Transformable for Object {
    fn transform(&self) -> &Transform {
        &self.transform
    }
    fn transform_mut(&mut self) -> &mut Transform {
        &mut self.transform
    }
}

impl Object {

    // pub fn plane() -> Self {

    //     let mut surfaces = Vec::new();
    //     // Creating the grid of triangles
    //     let SIDE_LENGTH = 200.0;
    //     for mut x in -250..250{
    //         for mut z in -250..250{
    //             let x = x as f64  * SIDE_LENGTH ;
    //             let z = z as f64  * SIDE_LENGTH ;
    //             surfaces.push([[x,0.0,z],[x + SIDE_LENGTH ,0.0,z] ,[x,0.0,z+SIDE_LENGTH]]);
    //             surfaces.push([[x+SIDE_LENGTH,0.0,z+SIDE_LENGTH],[x + SIDE_LENGTH ,0.0,z] ,[x,0.0,z+SIDE_LENGTH]]);

    //         }
    //     }
    //     return Object {
    //         scale: 1.0,
    //         transform: Transform {
    //             position: [0.0,0.0,0.0],
    //             quaternion: Quaternion::new(0.0,&[1.0,1.0,1.0]),
    //         },
    //         surfaces: surfaces,
    //     }
    // }



    pub fn new(scale: f64, position: &[f64; 3], quaternion: Quaternion, path:String) -> Self {
        if scale <= 0.0 {
            panic!("Attempted to instantiate with size <= 0.");
        }
        // Loading the surfaces from the object data   
        let (models,_materials) = load_obj(path, &LoadOptions {triangulate: true, ..Default::default()}).unwrap();
        let mut initial_surfaces: Vec<[[f64;3];3]> = Vec::new();
        for k in 0..models.len() {
            let mesh = &models[k].mesh;
            let mut index = 0;
            while index < mesh.indices.len(){

                let vertex1 = [
                mesh.positions[3 * mesh.indices[index] as usize ] as f64,
                mesh.positions[3 * mesh.indices[index] as usize+ 1] as f64,
                mesh.positions[3 * mesh.indices[index] as usize+ 2] as f64
                ];
                let vertex2 = [
                mesh.positions[3 * mesh.indices[index + 1] as usize ] as f64,
                mesh.positions[3 * mesh.indices[index + 1] as usize  + 1] as f64,
                mesh.positions[3 * mesh.indices[index + 1] as usize  + 2] as f64
                ];
                let vertex3 = [
                mesh.positions[3 * mesh.indices[index + 2] as usize ] as f64,
                mesh.positions[3 * mesh.indices[index + 2] as usize  + 1] as f64,
                mesh.positions[3 * mesh.indices[index + 2] as usize  + 2] as f64
                ];
                initial_surfaces.push([vertex1,vertex2,vertex3]);

                index = index + 3;
            }
    }
        return Object {
            scale: scale,
            transform: Transform {
                position: *position,
                quaternion: quaternion,
            },
            surfaces: initial_surfaces,
        }
    }

    pub fn get_surfaces(&self) -> Vec<[[f64;3];3]> {
        
        let (x, y, z) = (
            self.transform.position[0],
            self.transform.position[1],
            self.transform.position[2],
        );

        let scale = self.scale;
        let quat = self.transform.quaternion;
        let surfaces = self
            .surfaces
            .iter()
            .map(|surface| {
                let mut rotated_surface = [
                    // v1 -> v1'
                    ( (quat * Quaternion::from(&surface[0])) *quat.get_inverse() ).to_vector(),
                    // v2 -> v2'
                    ( (quat * Quaternion::from(&surface[1])) *quat.get_inverse() ).to_vector(),
                    // v3 -> v3'
                    ( (quat * Quaternion::from(&surface[2])) *quat.get_inverse() ).to_vector()
                ];
                // Translating
                for vertex in rotated_surface.iter_mut() {
                    vertex[0] += x;
                    vertex[1] += y;
                    vertex[2] += z;
                }
                // Scaling
                for vertex in rotated_surface.iter_mut() {
                    vertex[0] *= scale;
                    vertex[1] *= scale;
                    vertex[2] *= scale;
                }
                return rotated_surface;
            })
            .collect();
        return surfaces;
    }
    
}
