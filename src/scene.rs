
use minifb::*;
use crate::objects::*;
use crate::camera::*;
use std::f64::consts::PI;


struct LightSource {
    transform: Transform, 
    intensity: f64,     
    color: [u8;3]
}

struct Scene {
    pub objects: Vec<Object>,
    pub light_sources: Vec<LightSource>,
    pub cameras: Vec<Camera>,
}

impl Scene {

    pub fn new() -> Scene{
        return Scene {
            objects: Vec::new(),
            light_sources : Vec::new(),
            cameras: Vec::new()
        }
    }

}