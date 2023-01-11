mod objects;
mod camera;
use minifb::*;
use objects::*;
use camera::*;
use std::f64::consts::PI;

const WIDTH: usize = 600;
const HEIGHT: usize = 600;
const SECONDS_PER_FRAME: f32 = 0.02;


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

    let mut camera1 = Camera::new(&[0.0, 0.0, 0.0], &WIDTH, &HEIGHT);
    let mut world = Vec::new(); 

    // world.push(Cube::new(
    //     10.0,
    //     &[-50.0, 50.0, 0.0],
    //     Quaternion::new(0.0, &[1.0, 1.0, 1.0]),
    // ));  
        world.push(Cube::new(
            300.0,
            &[0.0,0.0, 0.0],
            Quaternion::new(5.0, &[1.0, 1.0, 0.0])));
    


    // Creating an empty window buffer for minifb to update the window with
    let mut window_buffer: Vec<u32> = Vec::new();
    // (Optional) Limit the window update rate to control CPU usage
    window.limit_update_rate(Some(std::time::Duration::from_secs_f32(SECONDS_PER_FRAME)));

    // Main window loop
    while window.is_open() && !window.is_key_down(Key::Escape) {
        // Rotate every cube +PI/200 radians about the vector <0,1,1>
        for cube in &mut world {
            cube.rotate(&mut Quaternion::new(PI / 200.0, &[0.5, 0.5, 0.5]));
        }
        // Take a snapshot with the camera
        camera1.update_buffer_with_surfaces(&world);
        //camera1.update_buffer_with_vertices(&world);
        // Convert the snapshot to a window_buffer to be used with minifb crate
        window_buffer.clear(); 
        for y in 0..HEIGHT {
            for x in 0..WIDTH {
                let pixel = camera1.buffer()[x][y] as u32;
                window_buffer.push((pixel | pixel << 8 | pixel << 16) as u32)
            }
        }
        // Update the window with the prepared snapshot 
        window.update_with_buffer(&window_buffer, WIDTH, HEIGHT).unwrap();
    }
}
