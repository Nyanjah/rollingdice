

#[derive(Copy,Clone)]
#[derive(Debug)]
pub struct Vertex {
    pub original_pos: [f64;3],
    pub pos: [f64;3],
    pub normal :[f64;3],
    pub tex_coords:[f64;2]
    
}

impl Vertex {

    pub fn new(pos:[f64;3] ,normal: [f64;3], tex_coords:[f64;2]) -> Vertex{
        return Vertex {
            original_pos: pos,
            pos: pos,
            normal: normal,
            tex_coords : tex_coords
        }
    }

    pub fn interpolate_attributes(self, vertex: &Vertex, val: f64) -> Vertex {
        return Vertex {
            original_pos: [lerp(self.original_pos[0],vertex.original_pos[0],val), lerp(self.original_pos[1],vertex.original_pos[1],val),lerp(self.original_pos[2],vertex.original_pos[2],val)],
            pos: [lerp(self.pos[0],vertex.pos[0],val), lerp(self.pos[1],vertex.pos[1],val),lerp(self.pos[2],vertex.pos[2],val)],
            normal: [lerp(self.normal[0],vertex.normal[0],val), lerp(self.normal[1],vertex.normal[1],val),lerp(self.normal[2],vertex.normal[2],val)],
            tex_coords: [lerp(self.tex_coords[0],vertex.tex_coords[0],val) , lerp(self.tex_coords[1],vertex.tex_coords[1],val)]
        }
    }  
    
    // Takes in the vertex and outputs its RGBA color encoded as a u32.
    pub fn shader(&self) -> u32 {

        let light: [f64;3] = [0.0,-100.0,-50.0];
        // difference  = light_source - pos
        let mut difference =  [light[0] - self.original_pos[0] ,light[1] - self.original_pos[1], light[2] - self.original_pos[2]];
        // normalizing the difference vector
        let magnitude = (difference[0].powi(2) + difference[1].powi(2) + difference[2].powi(2)).sqrt();
        difference = [difference[0] / magnitude, difference[1]/magnitude,difference[2]/magnitude];
        
        let dot = difference[0] * self.normal[0] + difference[1] * self.normal[1] + difference[2] * self.normal[2];
        if dot > 0.0 {
            let value =  (((1000000.0 / (magnitude).powi(2))) % 255.0 * dot) as usize;
            return (value << 16 | value << 8 | value  ) as u32;
        }
        else {
            let value = 0 as usize;
            return (value << 16 | value  ) as u32;
        }
        
        


    }

}

fn lerp(a: f64, b: f64, t: f64) -> f64 {
    return a * t + b * (1.0 - t);
}