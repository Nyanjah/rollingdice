use std::iter::zip;
#[derive(Copy,Clone)]
struct Quaternion{
    angle:f32,
    x:f32,
    y:f32,
    z:f32
}

impl std::ops::Mul for Quaternion{
    type Output = Quaternion;
    fn mul(self, rhs: Self) -> Self::Output {
        let product = Quaternion {
            angle: self.angle * rhs.angle - (self.x*rhs.x + self.y*rhs.y + self.z*rhs.z),
            x: self.angle * rhs.x + self.x * rhs.angle + self.y * rhs.z - self.z * rhs.y,
            y: self.angle * rhs.y - self.x * rhs.z + self.y * rhs.angle + self.z*rhs.x,
            z: self.angle * rhs.z + self.x * rhs.y - self.y * rhs.x + self.z * rhs.angle
        };
        return product;
    }
}

// impl std::ops::Add for Quaternion{
//     type Output = Quaternion;
//     fn add(self, rhs: Self) -> Self::Output {
        
//     }
// }

struct Cube {
    side_length: f32,   
    position: [f32;3],     // Position vector tracks center of cube
    quaternion: Quaternion // Rotation quaternion to track orientation in world
}

impl Cube {
    pub fn new(side_length: f32, position: [f32;3], quaternion:Quaternion) -> Self {
        if side_length <= 0.0 {
            panic!("Attempted to instantiate cube of side length <= 0.");
        }
        Cube {
            side_length:side_length,
            position:(position),
            quaternion: quaternion,
        }
    }

    pub fn translate(&mut self, vector:&[f32;3]){
        for (mut pos,translation) in zip(self.position,vector){
            pos += translation;
        }
    }

    pub fn rotate(&mut self,quaterion: Quaternion){          
        self.quaternion = self.quaternion * quaterion;
    }
}