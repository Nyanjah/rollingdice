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

//NOTE: Assuming camera is NOT rotated for current implementation...
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

    pub fn update_buffer_with_vertices(&mut self, world: &Vec<Cube>) {
        // Make sure the camera's buffer is empty
        self.clear_buffer();
        // Getting the (x,y) position of the top-left most point of the camera
        for cube in &*world {
            // Getting the vertices of the cube
            let vertices = cube.get_vertices();
            // Getting the positions of the center, top left, and bottom right of the camera grid
            let (x_0, y_0) = (self.transform.position[0], self.transform.position[1]);
            let top_left = (
                x_0 - 0.5 * self.width as f64,
                y_0 + 0.5 * self.height as f64,
            );
            let bottom_right = (
                x_0 + 0.5 * self.width as f64,
                y_0 - 0.5 * self.height as f64,
            );
            // Bounds checking vertices to make sure they land in the grid
            for vertex in vertices
                .iter()
                .filter(|vertex| {
                    vertex[0] < bottom_right.0 - 2.0
                        && vertex[0] > top_left.0 + 2.0
                        && vertex[1] < top_left.1 - 2.0
                        && vertex[1] > bottom_right.1 + 2.0
                })
                .map(|vertex| -> (usize, usize) {
                    // Mapping points from world to local screen layout (0,0) in top-left
                    return (
                        (vertex[0] - (x_0 - 0.5 * self.width as f64)).abs() as usize,
                        (vertex[1] - (y_0 + 0.5 * self.height as f64)).abs() as usize,
                    );
                })
            {
                self.buffer[vertex.0 as usize][vertex.1 as usize] = 255 as u8;
            }
        }
    }

    pub fn in_grid(&self, point: &[f64; 3]) -> bool {
        let (x_0, y_0) = (self.transform.position[0], self.transform.position[1]);
        let top_left = (
            x_0 - 0.5 * self.width as f64,
            y_0 + 0.5 * self.height as f64,
        );
        let bottom_right = (
            x_0 + 0.5 * self.width as f64,
            y_0 - 0.5 * self.height as f64,
        );
        if point[0] < bottom_right.0
            && point[0] > top_left.0
            && point[1] < top_left.1
            && point[1] > bottom_right.1
        {
            true
        } else {
            false
        }
    }

    pub fn update_buffer_with_surfaces(&mut self, world: &Vec<Cube>) {
        // Make sure the camera's buffer is empty
        self.clear_buffer();
        // Getting the (x,y) position of the top-left most point of the camera
        for cube in &*world {
            // Getting the vertices of the cube
            let mut surfaces = cube.get_tesselation();
            // Getting the position of the center of the camera grid
            let (x_0, y_0) = (self.transform.position[0], self.transform.position[1]);
            let surfaces: Vec<[[usize; 3]; 3]> = surfaces
                .iter_mut()
                // filtering out surfaces which shouldnt be rendered
                .filter(|surface| {
                    // Checking if all points in the surface are within the grid
                    self.in_grid(&surface[0])
                        && self.in_grid(&surface[1])
                        && self.in_grid(&surface[2])
                })
                // Mapping the coordinates of the points which make up the surfaces to pixel-grid coords
                .map(|surface| -> [[usize; 3]; 3] {
                    let mut grid_surface: [[usize; 3]; 3] = [[0; 3]; 3];
                    for i in 0..grid_surface.len() {
                        grid_surface[i] = [
                            (surface[i][0] as f64 - (x_0 - 0.5 * self.width as f64)).abs() as usize,
                            (surface[i][1] as f64 - (y_0 + 0.5 * self.height as f64)).abs()
                                as usize,
                            255 - (self.transform.position[2] - surface[i][2]).abs()
                                as usize,
                        ]
                    }
                    return grid_surface;
                })
                .collect();
            // Applying the Bresenham algorithm to plot lines at the pixel granularity
            // Drawing the triangle:
            for surface in surfaces {
                // Line P0->P1
                self.plot_line(
                    &surface[0][0],
                    &surface[0][1],
                    &surface[1][0],
                    &surface[1][1],
                    &surface[0][2],
                    &surface[1][2],
                );
                // Line P1->P2
                self.plot_line(
                    &surface[1][0],
                    &surface[1][1],
                    &surface[2][0],
                    &surface[2][1],
                    &surface[1][2],
                    &surface[2][2],
                );
                // Line P2->P0
                self.plot_line(
                    &surface[2][0],
                    &surface[2][1],
                    &surface[0][0],
                    &surface[0][1],
                    &surface[2][2],
                    &surface[0][2],
                );
            }
        }
    }

    pub fn buffer(&self) -> &Vec<Vec<u8>> {
        return &self.buffer;
    }
    pub fn clear_buffer(&mut self) {
        self.buffer = vec![vec![0; self.height]; self.width];
    }

    fn plot_line(
        &mut self,
        x_0: &usize,
        y_0: &usize,
        x_1: &usize,
        y_1: &usize,
        dist_0: &usize,
        dist_1: &usize,
    ) {
        let mut pixels_to_plot: Vec<(usize, usize, usize)> = Vec::new();

        if (*y_1 as i32 - *y_0 as i32).abs() < (*x_1 as i32 - *x_0 as i32).abs() {
            if x_0 > x_1 {
                pixels_to_plot = plot_line_low(x_1, y_1, x_0, y_0, dist_0, dist_1);
            } else {
                pixels_to_plot = plot_line_low(x_0, y_0, x_1, y_1, dist_0, dist_1);
            }
        } else {
            if y_0 > y_1 {
                pixels_to_plot = plot_line_high(x_1, y_1, x_0, y_0, dist_0, dist_1);
            } else {
                pixels_to_plot = plot_line_high(x_0, y_0, x_1, y_1, dist_0, dist_1);
            }
        }
        for pixel in pixels_to_plot {
            // TODO: Implement z-buffer for depth and lerp between colors
            self.buffer[pixel.0 as usize][pixel.1 as usize] = pixel.2 as u8;
            
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
    let mut D = 2 * dy - dx;
    let mut y = *y_0 as i32;
    for x in *x_0..*x_1 {
        output.push((
            x as usize,
            y as usize,
            lerp(
                *dist_0 as f64,
                *dist_1 as f64,
                (x - x_0) as f64 / (x_1 - x_0) as f64,
            ) as usize,
        ));
        if D > 0 {
            y = y + y_i;
            D = D + (2 * (dy - dx));
        } else {
            D = D + 2 * dy;
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
    for y in *y_0..*y_1 {
        output.push((
            x as usize,
            y as usize,
            lerp(
                *dist_0 as f64,
                *dist_1 as f64,
                (y - *y_0) as f64 / (y_1 - y_0) as f64,
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
