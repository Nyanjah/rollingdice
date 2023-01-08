mod physics;
mod rasterize;
use minifb::*;
use physics::*;
use std::f64::consts::PI;

const WIDTH: usize = 100;
const HEIGHT: usize = 100;
const SECONDS_PER_FRAME: f32 = 0.2;

struct Camera {
    transform: Transform,
    buffer: Vec<Vec<u8>>,
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
//NOTE: Assuming camera is NOT rotated for current implementation...
impl Camera {
    fn new(position: &[f64; 3], width: &usize, height: &usize) -> Self {
        return Camera {
            transform: Transform {
                position: *position,
                quaternion: Quaternion::new(0.0, &[1.0, 1.0, 1.0]),
            },
            buffer: vec![vec![0;*height];*width],
            width: *width,
            height: *height,
        };
    }

    fn update_buffer(&mut self, world: &Vec<Cube>) {
        // Getting the (x,y) position of the top-left most point of the camera
        let (x_0, y_0) = (
            self.transform.position[0] - 0.5 * self.width as f64,
            self.transform.position[1] + 0.5 * self.height as f64,
        );
        for cube in &*world{
        let vertices = cube.get_vertices();

        // For each (x,y) pixel in the camera's grid:
        for y in 0..self.height {
            for x in 0..self.width {
                // Calculate where the current sample point lies in relation to the world
                let (x_i, y_i) = (x_0 + x as f64 + 0.5, y_0 - y as f64 - 0.5);
                for &vertex in &vertices {
                    // If the vertex is within 1-pixel radius of the sample grid point:
                    if ((vertex[0] - x_i).powi(2) + (vertex[1] - y_i).powi(2)).sqrt() <= 1.0 {
                        // Add the grayscale pixel value to the buffer
                        self.buffer[x][y] = 255 as u8;
                        break
                    }
                    else {
                        //self.buffer[x][y] = 0 as u8;
                    }
                    
                }
            }
        }
    }
    }
    
    pub fn buffer(&self) -> &Vec<Vec<u8>>{
        return &self.buffer;
    }
    pub fn clear_buffer(&mut self){
        self.buffer = vec![vec![0;self.height];self.width];
        
    }
}


fn main() {
    // Setting up the window
    let mut window = Window::new(
        "Rasterization",
        WIDTH,
        HEIGHT,
        WindowOptions {
            borderless: false,
            title: true,
            resize: true,
            scale: Scale::X4,
            scale_mode: ScaleMode::Stretch,
            transparency: false,
            none: false,
            topmost: false,
        },
    )



    .unwrap_or_else(|e| {
        panic!("{}", e);
    });
    let initial_rotation = Quaternion::new(0.0, &[1.0, 1.0, 1.0]);
    let cube1 = Cube::new(35.0, &[0.0, 0.0, 0.0], initial_rotation);
    let cube2 = Cube::new(35.0/2.0, &[40.0, -50.0, 0.0], initial_rotation);
    let cube3 = Cube::new(35.0/3.0, &[-40.0, 0.0, 0.0], initial_rotation);
    let mut camera = Camera::new(&[0.0, 0.0, 10.0], &WIDTH, &HEIGHT);
    let rotation1 = Quaternion::new(PI / 10.0, &[1.0, 1.0, 1.0]);
    let rotation2 = Quaternion::new(PI / 20.0, &[-1.0, -1.0, 1.0]);
    let mut world = Vec::new();
    world.push(cube1);
    world.push(cube2);
    world.push(cube3);
   
    let mut window_buffer: Vec<u32> = Vec::new();
    window.limit_update_rate(Some(std::time::Duration::from_secs_f32(SECONDS_PER_FRAME)));
    while window.is_open() && !window.is_key_down(Key::Escape) {
        for cube in &mut world {
            cube.rotate(&rotation1);
            cube.rotate(&rotation2);
        }
        camera.update_buffer(&world);
        
        for x in 0..WIDTH{
            for y in 0..HEIGHT{
                let pixel = camera.buffer()[x][y] as u32;
                window_buffer.push(
                    (pixel | pixel << 8 | pixel << 16) as u32
                )
            }
        }

        window.update_with_buffer(&window_buffer, WIDTH, HEIGHT).unwrap();
        camera.clear_buffer();
        window_buffer.clear();
    }
}
