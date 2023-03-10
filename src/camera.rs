use crate::objects::*;
pub struct Camera {
    transform: Transform,
    buffer: Vec<Vec<u8>>,
    z_buffer: Vec<Vec<u8>>,
    width: usize,
    height: usize,
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
        //x >= -1.0 && x <= 1.0 && y >= -1.0 && y <= 1.0 
        // z coord should be negative, otherwise it is behind the camera
        //&& 
        z < -0.0 {true}
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

    pub fn buffer(&self) -> &Vec<Vec<u8>> {
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
        // Line P1->P2
        self.plot_line(
            &point1[0],
            &point1[1],
            &point2[0],
            &point2[1],
            &point1[2],
            &point2[2],
            &mut triangle_buffer,
        );
        // Line P2->P3
        self.plot_line(
            &point2[0],
            &point2[1],
            &point3[0],
            &point3[1],
            &point2[2],
            &point3[2],
            &mut triangle_buffer,
        );
        // Line P3->P1
        self.plot_line(
            &point3[0],
            &point3[1],
            &point1[0],
            &point1[1],
            &point3[2],
            &point1[2],
            &mut triangle_buffer,
        );
        // Iterating over the sub-section the screen is drawn to and applying the scanline algorithm
        // for y in 0..=(y_range) {
        // 'draw_loop:for x in 0..(x_range) {
        //         // If there is a pixel drawn at [x][y] & not at [x+1][y]
        //         if triangle_buffer[x][y].1 == true && triangle_buffer[x+1][y].1 == false {
        //             let first_pixel = (x, y);
        //             for x_i in ((first_pixel.0 + 1)..(x_range)).rev() {
        //                 // If there is a not a pixel drawn at [x_i][y] and one drawn at [x_i+1][y]
        //                 if triangle_buffer[x_i][y].1 == false && triangle_buffer[x_i+1][y].1 == true {
        //                     let second_pixel = (x_i+1, y);
        //                     // Plot the line between the two pixels
        //                     self.plot_line(
        //                         &first_pixel.0,
        //                         &first_pixel.1,
        //                         &second_pixel.0,
        //                         &second_pixel.1,
        //                         &(triangle_buffer[first_pixel.0][first_pixel.1].0 as usize),
        //                         &(triangle_buffer[second_pixel.0][second_pixel.1].0 as usize),
        //                         &mut triangle_buffer,
        //                     );
        //                     // Exit the loop for the current y-value
        //                     break 'draw_loop;
        //                 }
        //             }
        //         }
        //     }
        // }
        // Inserting the triangle into the camera's buffer
        for x in 1..(x_range) {
            if x + x_min >= self.width {break}
            for y in 1..(y_range){
            if y + y_min >= self.height {break}
            // If the pixel is closer (brighter) than the one in the z_buffer, add it to the z_buffer and draw it
                if triangle_buffer[x][y].0 > self.z_buffer[x+x_min][y+y_min]{
                    self.z_buffer[x+x_min][y+y_min] = triangle_buffer[x][y].0;
                    self.buffer[x+x_min][y+y_min] = triangle_buffer[x][y].0;
                }
           }
        }
    }

    fn plot_line(
        &mut self,
        x_0: &usize,
        y_0: &usize,
        x_1: &usize,
        y_1: &usize,
        alpha_1: &usize,
        apha_0: &usize,
        buffer: &mut Vec<Vec<(u8,bool)>>,
    ) {

        let pixels_to_plot:Vec<(usize, usize, usize)> = {
            
        if (*y_1 as i32 - *y_0 as i32).abs() < (*x_1 as i32 - *x_0 as i32).abs() {
            if x_0 > x_1 {
                plot_line_low(x_1, y_1, x_0, y_0, alpha_1, apha_0)
            } else {
                plot_line_low(x_0, y_0, x_1, y_1, apha_0, alpha_1)
            }
        } else {
            if y_0 > y_1 {
                plot_line_high(x_1, y_1, x_0, y_0, alpha_1, apha_0)
            } else {
                plot_line_high(x_0, y_0, x_1, y_1, apha_0, alpha_1)
            }
        }};
        let x_bound = buffer.len();
        let y_bound = buffer[0].len();
        for pixel in pixels_to_plot.iter().filter(|pixel| {
            if pixel.0 > x_bound - 1 || pixel.1 > y_bound -1 {
                false
            }
            else {true}
        })  {
            
            // Drawing the pixel to the buffer
            buffer[pixel.0 as usize][pixel.1 as usize].0 = pixel.2 as u8;
            // Setting the pixel_present flag to True
            buffer[pixel.0 as usize][pixel.1 as usize].1 = true;

        }
    }
}

// Takes in two points and returns a vector of points to color in the raster
fn plot_line_low(
    x_0: &usize,
    y_0: &usize,
    x_1: &usize,
    y_1: &usize,
    dist_0: &usize,
    dist_1: &usize,
) -> Vec<(usize, usize, usize)> {
    let mut output = Vec::new();
    let dx: i32 = *x_1 as i32 - *x_0 as i32;
    let mut dy: i32 = *y_1 as i32 - *y_0 as i32;
    let mut y_i: i32 = 1;
    if dy < 0 {
        y_i = -1;
        dy = -1 * dy as i32;
    }
    let mut d = 2 * dy - dx;
    let mut y = *y_0 as i32;
    for x in *x_0..=*x_1 {
        output.push((
            x as usize,
            y as usize,
            lerp(
                *dist_0 as f64,
                *dist_1 as f64,
                ((x - x_0) as f64).abs() / ((x_1 - x_0) as f64).abs(),
            ) as usize,
        ));
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
    x_0: &usize,
    y_0: &usize,
    x_1: &usize,
    y_1: &usize,
    dist_0: &usize,
    dist_1: &usize,
) -> Vec<(usize, usize, usize)> {
    let mut output = Vec::new();
    let mut dx = *x_1 as i32 - *x_0 as i32;
    let dy = *y_1 as i32 - *y_0 as i32;
    let mut x_i = 1 as i32;
    if dx < 0 {
        x_i = -1;
        dx = -dx;
    }
    let mut d = (2 * dx) - dy;
    let mut x: i32 = *x_0 as i32;
    for y in *y_0..=*y_1 {
        output.push((
            x as usize,
            y as usize,
            lerp(
                *dist_0 as f64,
                *dist_1 as f64,
                ((y - *y_0) as f64).abs() / ((y_1 - y_0) as f64).abs(),
            ) as usize,
        ));
        if d > 0 {
            x = x + x_i;
            d = d + (2 * (dx - dy));
        } else {
            d = d + 2 * dx;
        }
    }
    return output;
}
