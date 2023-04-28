use crate::objects::*;
use crate::vertex::*;
pub struct Camera {
    transform: Transform,
    frame_buffer: Vec<Vec<u32>>,
    triangle_buffer: Vec<Vec<(Vertex, bool)>>,
    z_buffer: Vec<Vec<f64>>,
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

impl Camera {
    pub fn new(position: &[f64; 3], width: &usize, height: &usize) -> Self {
        return Camera {
            transform: Transform {
                position: *position,
                quaternion: Quaternion::new(0.0, &[1.0, 1.0, 1.0]),
            },
            frame_buffer: vec![vec![0; *height]; *width],
            z_buffer: vec![vec![f64::MAX; *height]; *width],
            triangle_buffer: vec![
                vec![
                    (
                        Vertex::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0]),
                        false
                    );
                    *height+1
                ];
                *width+1
            ],
            width: *width,
            height: *height,
        };
    }

    pub fn in_view(&self, vertex: &Vertex) -> bool {
        let (x, y, z) = (vertex.pos[0], vertex.pos[1], vertex.pos[2]);
        // x and y values should be withing [-1,1]
        if x >= -1.0 && x <= 1.0 && y >= -1.0 && y <= 1.0 &&
        // z coord should be negative, otherwise it is behind the camera
        // && 
        z < 0.0
        {
            true
        } else {
            false
        } // Otherwise they are outside the field of view
    }

    pub fn update_buffer_with_surfaces(&mut self, world: &Vec<Object>) {
        // Make sure the camera's buffer is empty
        self.clear_buffer();
        // Getting the (x,y) position of the top-left most point of the camera
        for object in world.iter() {
            // Getting the surfaces of the cube
            let triangles = object.get_triangles();
            let mut vertices = object.get_transformed_vertices();

            // Getting the position of the viewing frustum
            let (x_0, y_0, z_0) = (
                self.transform.position[0],
                self.transform.position[1],
                self.transform.position[2],
            );
            let quat = self.transform.quaternion;
            // Getting the basis vectors of the camera's coordinate system
            let mut x_basis = (quat * Quaternion::from(&[1.0, 0.0, 0.0])) * quat.get_inverse();
            x_basis.normalize_as_vector();
            let mut y_basis = (quat * Quaternion::from(&[0.0, 1.0, 0.0])) * quat.get_inverse();
            y_basis.normalize_as_vector();
            let mut z_basis = (quat * Quaternion::from(&[0.0, 0.0, -1.0])) * quat.get_inverse();
            z_basis.normalize_as_vector();

            // Note: The raster's normal vector is the negative of the z_basis because the camera's line of sight is along it's -z axis.
            // Getting the top-left most vertex of the raster after applying it's stored translation and rotation

            // Plane's equation: A(x-x0)+B(y-y0)+C(z-z0) = 0 where unit_norm = < A, B, C >

            'vertex_projection_loop: for mut vertex in vertices.iter_mut() {
                let temp_point = [
                    vertex.pos[0] - x_0,
                    vertex.pos[1] - y_0,
                    vertex.pos[2] - z_0,
                ];
                // Getting the points relative to the camera's position
                // Swapping from world coords to camera coords

                // x = V' * x_basis
                vertex.pos[0] = x_basis.x * temp_point[0]
                    + x_basis.y * temp_point[1]
                    + x_basis.z * temp_point[2];
                // y = V' * y_basis
                vertex.pos[1] = y_basis.x * temp_point[0]
                    + y_basis.y * temp_point[1]
                    + y_basis.z * temp_point[2];
                // z = V' * z_basis
                vertex.pos[2] = z_basis.x * temp_point[0]
                    + z_basis.y * temp_point[1]
                    + z_basis.z * temp_point[2];

                // PERSPECTIVE DIVIDE STEP
                vertex.pos[0] = vertex.pos[0] / (-1.0 * vertex.pos[2] as f64);
                vertex.pos[1] = vertex.pos[1] / (-1.0 * vertex.pos[2] as f64);

                // If the vertex would not be in view, discard that vertex and move on
                if !self.in_view(&vertex) {
                    // We need to set the vertex to something special to signify that is has been discarded.
                    vertex.pos[0] = 1234567890.123456789012345678901234567890;
                    // Skip to the next vertex
                    continue 'vertex_projection_loop;
                }
                // X-values
                vertex.pos[0] =
                    ((vertex.pos[0] * self.width as f64 / 2.0) + (self.width as f64 / 2.0)).round();
                // Y-values
                vertex.pos[1] =
                    ((vertex.pos[1] * self.height as f64 / 2.0) + (self.height as f64 / 2.0)).round();
            }

            // This should ONLY BE CALLED if ALL of the vertices the current triangle's indices correspond to are still valid
            for triangle in triangles {
                // We only draw the triangle if it doesn't contain any vertices which are marked as discarded with the special value
                if !triangle
                    .iter()
                    .any(|i| vertices[*i].pos[0] == 1234567890.123456789012345678901234567890)
                {
                    self.draw_triangle(
                        vertices[triangle[0]],
                        vertices[triangle[1]],
                        vertices[triangle[2]],
                    );
                }
            }
        }
        // Clear out the z-buffer for the frame
        for vector in &mut self.z_buffer {
            vector.fill(f64::MAX);
        }
        
        // Clear out the triangle buffer
        // self.triangle_buffer = vec![
        //     vec![
        //         (
        //             Vertex::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0]),
        //             false
        //         );
        //         self.height + 1
        //     ];
        //     self.width + 1
        // ]
        
    }

    pub fn export_frame(&self) -> &Vec<Vec<u32>> {
        return &self.frame_buffer;
    }

    pub fn clear_buffer(&mut self) {
        self.frame_buffer = vec![vec![0; self.height]; self.width];
    }

    fn draw_triangle(&mut self, vertex_1: Vertex, vertex_2: Vertex, vertex_3: Vertex) {
        // Extracting the three vertices which make up the triangle

        // Getting bounds for the subsection of the screen the triangle is drawn to
        let x_min = ((vertex_1.pos[0]).min(vertex_2.pos[0])).min(vertex_3.pos[0]) as usize;
        let x_max = ((vertex_1.pos[0]).max(vertex_2.pos[0])).max(vertex_3.pos[0]).min((self.width - 1) as f64) as usize;
        let y_min = ((vertex_1.pos[1]).min(vertex_2.pos[1])).min(vertex_3.pos[1]) as usize;
        let y_max = ((vertex_1.pos[1]).max(vertex_2.pos[1])).max(vertex_3.pos[1]).min((self.height - 1) as f64) as usize;

        // Plotting the three lines into the triangle buffer that make up the triangle's perimeter

        self.plot_line(&vertex_1, &vertex_2); // Line P1->P2
        self.plot_line(&vertex_2, &vertex_3); // Line P2->P3
        self.plot_line(&vertex_3, &vertex_1); // Line P3->P1

        // // Iterating over the sub-section the screen is drawn to and applying the scanline algorithm
        for y in y_min..=y_max {
            'draw_loop: for x in x_min..=x_max {
                // If there is a pixel drawn at [x][y] & not at [x+1][y]
                if self.triangle_buffer[x][y].1 && !self.triangle_buffer[x + 1][y].1 {
                    let first_pixel = (x, y);
                    for x_i in ((first_pixel.0 + 1)..(x_max)).rev() {
                        // If there is a not a pixel drawn at [x_i][y] and one drawn at [x_i+1][y]
                        if !self.triangle_buffer[x_i][y].1 && self.triangle_buffer[x_i + 1][y].1 {
                            let second_pixel = (x_i + 1, y);
                            // Plot the line between the two pixels
                            let first_point =  self.triangle_buffer[first_pixel.0][first_pixel.1].0;
                            let second_point = self.triangle_buffer[second_pixel.0][second_pixel.1].0;

                            self.plot_line(&first_point, &second_point);
                            // Exit the loop for the current y-value
                            break 'draw_loop;
                        }
                    }
                }
            }
        }

        // Inserting the finished triangle into the frame buffer using the triangle in the triangle buffer
        for x in x_min..=x_max {
            for y in y_min..=y_max {
                // If there is a pixel in the triangle buffer at location (x,y)
                if self.triangle_buffer[x][y].1 == true {
                    let vertex = self.triangle_buffer[x][y].0;
                    // Setting the pixel_present flag in the corresponding entry of triangle buffer to false
                    self.triangle_buffer[x][y].1 = false;
                    // If the pixel is closer than the one in the z_buffer
                    if vertex.pos[2].abs() < self.z_buffer[x][y] {
                        // Overwrite the value stored in the z-buffer with the pixel's distance
                        self.z_buffer[x][y] = vertex.pos[2].abs();
                        // Call the vertex shader to get the u32 RGB value to be placed in the final frame
                        self.frame_buffer[x][y] = vertex.shader();
                    }
                }
            }
        }
    }


// TODO: Switch to barycentric coodrinates for smoother 3-points interpolation.
    fn plot_line(&mut self, vertex_0: &Vertex, vertex_1: &Vertex) {
        // Extracting the vertex coordinates
        let (x_0, y_0) = (vertex_0.pos[0] as i32, vertex_0.pos[1] as i32);
        let (x_1, y_1) = (vertex_1.pos[0] as i32, vertex_1.pos[1] as i32);
        // plot_line_low and plot_line_high should place the pixels.

        if (y_1 as i32 - y_0 as i32).abs() < (x_1 as i32 - x_0 as i32).abs() {
            if x_0 > x_1 {
                self.plot_line_low(vertex_1, vertex_0);
            } else {
                self.plot_line_low(vertex_0, vertex_1);
            }
        } else {
            if y_0 > y_1 {
                self.plot_line_high(vertex_1, vertex_0);
            } else {
                self.plot_line_high(vertex_0, vertex_1);
            }
        };
    }

    fn plot_line_low(&mut self, vertex_0: &Vertex, vertex_1: &Vertex) {
        let (x_0, y_0) = (vertex_0.pos[0] as i32, vertex_0.pos[1] as i32);
        let (x_1, y_1) = (vertex_1.pos[0] as i32, vertex_1.pos[1] as i32);

        let dx: i32 = x_1 - x_0;
        let mut dy: i32 = y_1 - y_0;
        let mut y_i: i32 = 1;
        if dy < 0 {
            y_i = -1;
            dy = -1 * dy as i32;
        }
        let mut d = 2 * dy - dx;
        let mut y = y_0 as i32;
        for x in x_0..=x_1 {
           
            // Place the interpolated vertex directly into the triangle buffer
            let interpolation_value = 
            (((x - x_0) as f64).powi(2) + ((x - x_0) as f64).powi(2)) / (((x_1 - x_0) as f64).powi(2) + ((y_1 - y_0) as f64).powi(2));
            
            self.triangle_buffer[x as usize][y as usize].0 =
                vertex_0.interpolate_attributes(vertex_1, interpolation_value);

            self.triangle_buffer[x as usize][y as usize].0.pos[0] = x as f64;
            self.triangle_buffer[x as usize][y as usize].0.pos[1] = y as f64;

            // Setting the pixel-present flag to true
            self.triangle_buffer[x as usize][y as usize].1 = true;

            if d > 0 {
                y = y + y_i;
                d = d + (2 * (dy - dx));
            } else {
                d = d + 2 * dy;
            }
        }
    }

  
    fn plot_line_high(&mut self, vertex_0: &Vertex, vertex_1: &Vertex) {
        let (x_0, y_0) = (vertex_0.pos[0] as i32, vertex_0.pos[1] as i32);
        let (x_1, y_1) = (vertex_1.pos[0] as i32, vertex_1.pos[1] as i32);

        let mut dx = x_1 - x_0;
        let dy = y_1 - y_0;
        let mut x_i = 1 as i32;
        if dx < 0 {
            x_i = -1;
            dx = -dx;
        }
        let mut d = (2 * dx) - dy;
        let mut x: i32 = x_0;
        for y in y_0..=y_1 {
          
            // Place the interpolated vertex directly into the triangle buffer
            let interpolation_value = 
            (((x - x_0) as f64).powi(2) + ((x - x_0) as f64).powi(2)) / (((x_1 - x_0) as f64).powi(2) + ((y_1 - y_0) as f64).powi(2));


            
            self.triangle_buffer[x as usize][y as usize].0 =
                vertex_0.interpolate_attributes(vertex_1, interpolation_value);
            self.triangle_buffer[x as usize][y as usize].0.pos[0] = x as f64;
            self.triangle_buffer[x as usize][y as usize].0.pos[1] = y as f64;


            // Setting the pixel-present flag to true
            self.triangle_buffer[x as usize][y as usize].1 = true;

            if d > 0 {
                x = x + x_i;
                d = d + (2 * (dx - dy));
            } else {
                d = d + 2 * dx;
            }
        }
    }
}
