mod physics;
mod rasterize;
use minifb::*;
use physics::*;
use std::f64::consts::PI;

const WIDTH: usize = 800;
const HEIGHT: usize = 800;
const SECONDS_PER_FRAME: f32 = 0.02;

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
            buffer: vec![vec![0; *height]; *width],
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
        for cube in &*world {
            let vertices = cube.get_vertices();


            let (x_0, y_0) = (self.transform.position[0], self.transform.position[1]);
            let top_left = (x_0 - 0.5*self.width as f64,y_0 + 0.5*self.height as f64);
            let bottom_right = (x_0 + 0.5*self.width as f64, y_0 - 0.5 * self.height as f64);
            // Bounds checking 
            for vertex in vertices.iter().filter(|vertex| {
                vertex[0] < bottom_right.0 - 1.0 && vertex[0] > top_left.0 + 1.0
                && vertex[1] < top_left.1 - 1.0 && vertex[1] > bottom_right.1 + 1.0
            }) {
                // Rounding down to get [x][y] pos on screen
                let (x, y) = (
                    (vertex[0] - (x_0 - 0.5 * self.width as f64)).abs(),
                    (vertex[1] - (y_0 + 0.5 * self.height as f64)).abs(),
                );
                self.buffer[x as usize][y as usize] = 255 as u8;
            }
        }
    }

    pub fn buffer(&self) -> &Vec<Vec<u8>> {
        return &self.buffer;
    }
    pub fn clear_buffer(&mut self) {
        self.buffer = vec![vec![0; self.height]; self.width];
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
            scale: Scale::X1,
            scale_mode: ScaleMode::Stretch,
            transparency: false,
            none: false,
            topmost: false,
        },
    )
    .unwrap_or_else(|e| {
        panic!("{}", e);
    });
  

    let mut camera1 = Camera::new(&[0.0, 0.0, 10.0], &WIDTH, &HEIGHT);
   

    let mut world = Vec::new();

    for i in 1..=10000{
        world.push(Cube::new(0.1*i as f64, &[0.0, 1.0, 0.0]
            ,Quaternion::new(i as f64, &[0.0, 1.0, 0.0])))
    }

    let mut window_buffer: Vec<u32> = Vec::new();
    window.limit_update_rate(Some(std::time::Duration::from_secs_f32(SECONDS_PER_FRAME)));
    while window.is_open() && !window.is_key_down(Key::Escape) {
        for cube in &mut world {
            let (x, y, z) = (
                cube.transform().position[0],
                cube.transform().position[1],
                cube.transform().position[2],
            );
            cube.rotate(&Quaternion::new(PI / 60.0, &[0.0+x, 1.0+y, 0.0+z]));
            //cube.rotate(&Quaternion::new(PI / 20.0, &[-1.0+x, -1.0+y, 1.0+z]));
        }
        camera1.update_buffer(&world);

        for x in 0..WIDTH {
            for y in 0..HEIGHT {
                let pixel = camera1.buffer()[x][y] as u32;
                window_buffer.push((pixel | pixel << 8 | pixel << 16) as u32)
            }
        }

        window
            .update_with_buffer(&window_buffer, WIDTH, HEIGHT)
            .unwrap();
        camera1.clear_buffer();
        window_buffer.clear();
    }
}
