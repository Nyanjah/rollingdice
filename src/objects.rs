#![allow(dead_code)]
use tobj::*;

use crate::vertex::Vertex;

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
    triangles: Vec<[usize;3]>,
    vertices: Vec<Vertex>
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

    pub fn rotate_vector(&self, vector:&[f64;3]) -> [f64;3]{
       return ( (*self * Quaternion::from(vector)) *self.get_inverse() ).to_vector()
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

    pub fn new(scale: f64, position: &[f64; 3], quaternion: Quaternion, path:String) -> Self {
        if scale <= 0.0 {
            panic!("Attempted to instantiate with size <= 0.");
        }
        let (triangles,vertices) = load_obj_file(&path).unwrap();
        
        return Object {
            scale: scale,
            transform: Transform {
                position: *position,
                quaternion: quaternion,
            },
            // Triangles are a triplet of indices which are used to access vertices 
            triangles: triangles,
            // Vertices have positions, normals, and texture coordinates
            vertices: vertices
        }
    }

    pub fn get_triangles(&self) -> &Vec<[usize;3]>{
        return &self.triangles
    }
    pub fn get_transformed_vertices(&self) -> Vec<Vertex> {
        // x,y,z translation to be applied to each vertex
        let (x, y, z) = (
            self.transform.position[0],
            self.transform.position[1],
            self.transform.position[2],
        );
        // scaling factor to be applied to each vertex
        let scale = self.scale;
        // rotation to be applied to each vertex represented as a quaternion
        let quat = self.transform.quaternion;
    
        // Create a vector to store the transformed vertices
        let mut transformed_vertices = Vec::with_capacity(self.vertices.len());
    
        // Loop through each vertex, apply the transformation and store it in the vector
        // The new vertex is built piece-by-piece using the transformations stored in the object
        for vertex in &self.vertices {
            // Translate the vertex
            let translated_vertex = Vertex {
                pos: [
                    vertex.pos[0] + x,
                    vertex.pos[1] + y,
                    vertex.pos[2] + z,
                ],
                ..*vertex
            };
    
            // Scale the vertex
            let scaled_vertex = Vertex {
                pos: [
                    translated_vertex.pos[0] * scale,
                    translated_vertex.pos[1] * scale,
                    translated_vertex.pos[2] * scale,
                ],
                ..translated_vertex
            };
    
            // Rotate the vertex
            let mut rotated_vertex = Vertex {
                pos: quat.rotate_vector(&scaled_vertex.pos),
                normal: quat.rotate_vector(&vertex.normal),
                ..scaled_vertex
            };

            rotated_vertex.original_pos = rotated_vertex.pos;
    
            // Add the transformed vertex to the vector
            transformed_vertices.push(rotated_vertex);
        }
        // Return the new transformed vertices
        return transformed_vertices
    }
    

}

fn load_obj_file(path: &str) -> Result<(Vec<[usize; 3]>, Vec<Vertex>), String> {
    let mut triangles = Vec::new();
    let mut vertices = Vec::new();

    let options = LoadOptions::default();
    let (models, _) = load_obj(path, &options).unwrap();

    for model in models.iter() {
        let mesh = &model.mesh;

        // Iterate over the face indices and add triangles to the vector
        for i in (0..mesh.indices.len()).step_by(3) {
            let i1 = mesh.indices[i] as usize;
            let i2 = mesh.indices[i + 1] as usize;
            let i3 = mesh.indices[i + 2] as usize;
            triangles.push([i1, i2, i3]);
        }

        // Iterate over the vertex positions and create a Vertex struct for each vertex
        for i in 0..mesh.positions.len() / 3 {
            // Vertex positions
            let x = mesh.positions[i * 3]     as f64;
            let y = mesh.positions[i * 3 + 1] as f64;
            let z = mesh.positions[i * 3 + 2] as f64;

            // Vertex normals
            //let nx = mesh.normals[i * 3]      as f64;
            //let ny = mesh.normals[i * 3 + 1]  as f64;
            //let nz = mesh.normals[i * 3 + 2]  as f64;

            // Vertex U,V texture coordinates
            //let u = mesh.texcoords[i * 2]     as f64;
            //let v = mesh.texcoords[i * 2 + 1] as f64;

            let vertex = Vertex {
                original_pos: [x, y, z],
                pos: [x, y, z],
                // normal: [nx, ny, nz],
                // tex_coords: [u, v],
                normal: [0.0, 0.0, 0.0],
                tex_coords: [0.0, 0.0],
            };

            vertices.push(vertex);
        }
    }

    Ok((triangles, vertices))
}
