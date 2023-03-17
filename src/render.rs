use crate::objects::*;
pub struct Camera {
    transform: Transform,
    buffer: Vec<Vec<u32>>,
    z_buffer: Vec<Vec<u8>>,
    width: usize,
    height: usize,
}

// The properties of the fragments are determined by their corresponding geometric primitives
// In this case, they are determined by their corresponding triangular surface loaded from the .obj file.
pub struct Fragment{
    x: usize,       // x-coord in the raster
    y: usize,       // y-cord in the raster
    distance: f64,  // Distance to be used in the z-buffer to solve the visibility problem
    color: [f64;4]  // RGBA color
}
pub struct Surface{
    vertices: [Vertex;3], // Contains three vertices
    // Todo: Add more properties for improved graphics (diffuse color, specular color, other material properties, etc...)
}

pub struct Vertex{
    x:f64, // x-coord
    y:f64, // y-coord
    z:f64,  // z-coord
    normal: [f64;3],   // vertex-normal ( calculated by summing the surface normals, then normalizing to a unit vector)
    color: [f64;4]     // vertex-color ( from loaded .obj file )
    // Todo: This would be a good place to store texture coordinates
}

impl Transformable for Camera {
    fn transform(&self) -> &Transform {
        &self.transform
    }
    fn transform_mut(&mut self) -> &mut Transform {
        &mut self.transform
    }
}

pub fn lerp(a: f64, b: f64, t: f64) -> f64 {
    return a * t + b * (1.0 - t);
}

impl Camera {
    pub fn new(position: &[f64; 3], width: &usize, height: &usize) -> Self {
        return Camera {
            transform: Transform {
                position: *position,
                quaternion: Quaternion::new(0.0, &[1.0, 1.0, 1.0]),
            },
            buffer: vec![vec![0; *height]; *width],
            z_buffer: vec![vec![0; *height]; *width],
            width: *width,
            height: *height,
        };
    }

    pub fn in_view(&self, point: &[f64; 3]) -> bool {
        let (x,y,z) = (point[0],point[1],point[2]);
        // x and y values should be withing [-1,1]
        if 
        x >= -1.0 && x <= 1.0 && y >= -1.0 && y <= 1.0 &&
        // z coord should be negative, otherwise it is behind the camera
        //&& 
        z < 0.0 {true}
        else {
        false} // Otherwise they are outside the field of view
    }

    pub fn update_buffer_with_surfaces(&mut self, world: &Vec<Object>) {
        // Make sure the camera's buffer is empty
        self.clear_buffer();
        // Getting the (x,y) position of the top-left most point of the camera
        for object in &*world {

            // Getting the surfaces of the cube
            let surfaces = object.get_surfaces();

            // Getting the position of the viewing frustum
            let (x_0, y_0, z_0) = (self.transform.position[0], self.transform.position[1], self.transform.position[2]);

            // Getting the basis vectors of the camera's coordinate system
            let mut x_basis = (self.transform.quaternion * Quaternion::from(&[1.0,0.0,0.0])) * self.transform.quaternion.get_inverse();
            x_basis.normalize_as_vector();
            let mut y_basis = (self.transform.quaternion * Quaternion::from(&[0.0,1.0,0.0])) * self.transform.quaternion.get_inverse();
            y_basis.normalize_as_vector();
            let mut z_basis = (self.transform.quaternion * Quaternion::from(&[0.0,0.0,-1.0])) * self.transform.quaternion.get_inverse();
            z_basis.normalize_as_vector();
            // Note: The raster's normal vector is the negative of the z_basis because the camera's line of sight is along it's -z axis.
            // Getting the top-left most vertex of the raster after applying it's stored translation and rotation

            // Plane's equation: A(x-x0)+B(y-y0)+C(z-z0) = 0 where unit_norm = <A,B,C>  
            let surfaces: Vec<[[usize;3];3]> = surfaces
                .iter()
                // Get point coords w.r.t to the camera's position 
                // ( It has an inverted z-axis because its line of sight is positioned along the -z axis...)
                // So Vertex' = ( Vertex - V_0 ) with an inverted z-component
                .map(|surface| -> [[f64;3];3] {
                    let mut plane_surface:[[f64;3];3] = [[0.0;3];3];
                    let mut i = 0;
                    for point in surface{
                        let temp_point = [point[0]-x_0,point[1]-y_0,point[2]-z_0];
                        // Getting the points relative to the camera's position
                        let mut new_point = [0.0;3];
                        // Swapping from world coords to camera coords

                        // x = V' * x_basis
                        new_point[0] = x_basis.x * temp_point[0] + x_basis.y * temp_point[1] + x_basis.z * temp_point[2];
                        // y = V' * y_basis
                        new_point[1] = y_basis.x * temp_point[0] + y_basis.y * temp_point[1] + y_basis.z * temp_point[2];
                        // z = V' * z_basis
                        new_point[2] = z_basis.x * temp_point[0] + z_basis.y * temp_point[1] + z_basis.z * temp_point[2];

                        // PERSPECTIVE DIVIDE STEP
                        new_point[0] = new_point[0] / ( -1.0 * new_point[2]);
                        new_point[1] = new_point[1] / ( -1.0 * new_point[2]);
                        
                        // Creating the new surface using the new points
                        plane_surface[i] = [new_point[0], new_point[1], new_point[2]];
                        // Note: All points in the field of view should lie withing a [-1,1] square.
                        i = i + 1;
                    }
                    return plane_surface
                })

                // filtering out surfaces which are outside the field of view
                .filter(|surface| {
                    // Checking if all points in the surface are within the grid
                           self.in_view(&surface[0])
                        && self.in_view(&surface[1])
                        && self.in_view(&surface[2])
                })
                // Converting from NDC (Normalized Device Coordinates) to screen coordinates
                // and setting a temporary alpha value to the pixels in the surface directly
                .map(|surface| -> [[usize;3];3]{
                    let mut grid_surface:[[usize;3];3] = [[0;3];3];
                    //println!("{:?}",surface); 
                    let mut i = 0;
                    for point in surface{
                        // X-values
                        grid_surface[i][0] = ((point[0] * self.width as f64/2.0) + self.width as f64/2.0).round() as usize;
                        // Y-values
                        grid_surface[i][1] = ((point[1] * self.height as f64/2.0) + self.height as f64/2.0).round() as usize;
                        // orthogonal distance from raster mapped to a brightness 0-255
                        // grid_surface[i][2] = 255 - (point[2].abs()).clamp(0.0,255.0) as usize; 
                        grid_surface[i][2] = 255;  
                        i = i + 1;
                    }
                    return grid_surface
                })                
                .collect();
            // Drawing the triangle:
            for surface in surfaces {   
                self.draw_triangle(&surface);
            }
            // Clear out the z-buffer for the frame
            self.z_buffer =  vec![vec![0; self.height]; self.width];
        }
    }

    pub fn buffer(&self) -> &Vec<Vec<u32>> {
        return &self.buffer;
    }
    pub fn clear_buffer(&mut self) {
        self.buffer = vec![vec![0; self.height]; self.width];
    }
    fn draw_triangle(&mut self, surface: &[[usize; 3]; 3]) {
        let mut point1 = surface[0];
        let mut point2 = surface[1];
        let mut point3 = surface[2];

        // Getting bounds for the subsection of the screen the triangle is drawn to
        let x_min = ((point1[0]).min(point2[0])).min(point3[0]);
        let x_max = ((point1[0]).max(point2[0])).max(point3[0]);
        let y_min = ((point1[1]).min(point2[1])).min(point3[1]);
        let y_max = ((point1[1]).max(point2[1])).max(point3[1]);
        // Clamping the x and y ranges to limit the buffer size to the size of the screen
        let y_range = (y_max - y_min).clamp(0,self.height);
        let x_range = (x_max - x_min).clamp(0,self.width);

        for point in [&mut point1,&mut point2,&mut point3]{
            point[0] -= x_min;
            point[1] -= y_min;
        }
        // Making a temporary per-triangle buffer for intermediate processing
        // Stores the u8 (0-255) alpha value and a flag for if there is a pixel present at the location
        let mut triangle_buffer: Vec<Vec<(u8,bool)>> = vec![vec![(0,false); (y_range) + 1]; (x_range) + 1];

        // Plotting the three lines that make up the surfaces boundary
        
        
        // Inserting the triangle into the camera's buffer
    






    }

    // Takes in two fragments and returns a vector of fragments which make up the line
    fn get_line_fragments(
        &mut self,
        fragment_0: &Fragment,
        fragment_1: &Fragment,
    ) -> Vec<Fragment>
    {
        let (x_0,y_0) = (fragment_0.x,fragment_0.y);
        let (x_1,y_1) = (fragment_1.x,fragment_1.y);
        let mut fragments = Vec::new();
        if (y_1 as i32 - y_0 as i32).abs() < (x_1 as i32 - x_0 as i32).abs() {
            if x_0 > x_1 {
                fragments = plot_line_low(fragment_0,fragment_1);
            } else {
                fragments = plot_line_low(fragment_0,fragment_1);
            }
        } else {
            if y_0 > y_1 {
                fragments = plot_line_high(fragment_0,fragment_1);
            } else {
                fragments = plot_line_high(fragment_0,fragment_1);
            }
        };
        return fragments
        }
    }

// Takes in two points and returns a vector of points to color in the raster
fn plot_line_low(
fragment_0:&Fragment,
fragment_1:&Fragment
) -> Vec<Fragment> {

    let (x_0,y_0) = (fragment_0.x,fragment_0.y);
    let (x_1,y_1) = (fragment_1.x,fragment_1.y);

    let mut output = Vec::new();
    let dx: i32 = x_1 as i32 - x_0 as i32;
    let mut dy: i32 = y_1 as i32 - y_0 as i32;
    let mut y_i: i32 = 1;
    if dy < 0 {
        y_i = -1;
        dy = -1 * dy as i32;
    }
    let mut d = 2 * dy - dx;
    let mut y = y_0 as i32;
    for x in x_0..=x_1 {
        let difference = ((x - x_0) as f64).abs() / ((x_1 - x_0) as f64).abs();
        output.push(
            Fragment { x: x as usize, y: y as usize,
                distance: lerp(fragment_0.distance,fragment_1.distance,difference),
                color: [
                    lerp(fragment_0.color[0],fragment_1.color[0],difference),
                    lerp(fragment_0.color[1],fragment_1.color[1],difference),
                    lerp(fragment_0.color[2],fragment_1.color[2],difference),
                    0.0
                ] 
            }
        );
        if d > 0 {
            y = y + y_i;
            d = d + (2 * (dy - dx));
        } else {
            d = d + 2 * dy;
        }
    }
    return output;
}

// Takes in two points and returns a vector of points to color in the raster
fn plot_line_high(
fragment_0: &Fragment,
fragment_1: &Fragment
) -> Vec<Fragment> {

    let (x_0,y_0) = (fragment_0.x,fragment_0.y);
    let (x_1,y_1) = (fragment_1.x,fragment_1.y);

    let mut output = Vec::new();
    let mut dx = x_1 as i32 - x_0 as i32;
    let dy = y_1 as i32 - y_0 as i32;
    let mut x_i = 1 as i32;
    if dx < 0 {
        x_i = -1;
        dx = -dx;
    }
    let mut d = (2 * dx) - dy;
    let mut x: i32 = x_0 as i32;
    for y in y_0..=y_1 {
        let difference = ((x - x_0 as i32) as f64).abs() / ((x_1 - x_0) as f64).abs();
        output.push(
            Fragment { x: x as usize, y: y as usize,
                distance: lerp(fragment_0.distance,fragment_1.distance,difference),
                color: [
                    lerp(fragment_0.color[0],fragment_1.color[0],difference),
                    lerp(fragment_0.color[1],fragment_1.color[1],difference),
                    lerp(fragment_0.color[2],fragment_1.color[2],difference),
                    0.0
                ] 
            }
        );
        if d > 0 {
            x = x + x_i;
            d = d + (2 * (dx - dy));
        } else {
            d = d + 2 * dx;
        }
    }
    return output;
}
