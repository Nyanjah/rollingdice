

#[derive(Copy,Clone)]
#[derive(Debug)]
pub struct Vertex {
    pub pos: [f64;3],
    pub normal :[f64;3],
    pub tex_coords:[f64;2]
    
}

impl Vertex {

    pub fn new(pos:[f64;3] ,normal: [f64;3], tex_coords:[f64;2]) -> Vertex{
        return Vertex {
            pos: pos,
            normal: normal,
            tex_coords : tex_coords
        }
    }


    pub fn interpolate_attributes(self, vertex: &Vertex, val: f64) -> Vertex {
        return Vertex {
            pos: [lerp(self.pos[0],vertex.pos[0],val), lerp(self.pos[1],vertex.pos[1],val),lerp(self.pos[2],vertex.pos[2],val)],
            normal: [lerp(self.normal[0],vertex.normal[0],val), lerp(self.normal[1],vertex.normal[1],val),lerp(self.normal[2],vertex.normal[2],val)],
            tex_coords: [lerp(self.tex_coords[0],vertex.tex_coords[0],val) , lerp(self.tex_coords[1],vertex.tex_coords[1],val)]
        }
    }  
    
    // Takes in the vertex and outputs its RGBA color encoded as a u32.
    pub fn shader(&self) -> u32 {

        // Return (255,255,255)
        return 255 << 16 | 255 << 8 | 255 ;

    }

}

fn lerp(a: f64, b: f64, t: f64) -> f64 {
    return a * t + b * (1.0 - t);
}