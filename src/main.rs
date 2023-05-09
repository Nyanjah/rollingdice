mod objects;
mod camera;
mod vertex;
mod scene;
use minifb::*;
use objects::*;
use camera::*;
use std::f64::consts::PI;

const WIDTH: usize = 1080;
const HEIGHT: usize = 720;
const SECONDS_PER_FRAME: f32 = 1.0/60.0; // MAX 60 FPS

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

    let mut camera1 = Camera::new(&[0.0, -10.0, -30.0], &WIDTH, &HEIGHT);
   
    let mut world = Vec::new(); 


    world.push(Object::new(
        1.0,
        &[0.0 as f64 ,0.0 , 0.0 as f64],
        Quaternion::new(PI, &[1.0, 0.0, 0.0]),
        "./models/skull.obj".to_string()
    ));


    for object in &mut world{
        object.load_texture("textures./skull_texture.jpg");
        object.rotate(-1.0 * PI/2.0 ,[1.0,0.0,0.0]);
    }

    // Creating an empty window buffer for minifb to update the window with
    let mut window_buffer: Vec<u32> = Vec::new();
    
    // (Optional) Limit the window update rate to control CPU usage
    window.limit_update_rate(Some(std::time::Duration::from_secs_f32(SECONDS_PER_FRAME)));

    // Main window loop
    let mut t : f64 = 0.0;
    while window.is_open() && !window.is_key_down(Key::Escape) {
        // Rotate every cube +PI/200 radians about the vector <0,1,1>
        for object in &mut world {
            object.rotate(PI / 200.0, [0.0, 0.0, 1.0]);
            object.translate(&[0.0,0.0,t.sin()/5.0]);
        }
        t = t + 0.3;
        // Take a snapshot with the camera
        camera1.update_buffer_with_surfaces(&world);
     
        // Convert the snapshot to a window_buffer to be used with minifb crate
        window_buffer.clear(); 
        for y in 0..HEIGHT {
            for x in 0..WIDTH {
                let pixel = camera1.export_frame()[x][y];
                window_buffer.push(pixel)
            }
        }

        // Update the window with the prepared frame
        window.update_with_buffer(&window_buffer, WIDTH, HEIGHT).unwrap();
  
        if window.is_key_down(Key::D){
            // Go left
            camera1.translate(&[-1.0,0.0,0.0]);
        }
        if window.is_key_down(Key::A){
            // Go right
            camera1.translate(&[1.0,0.0,0.0]);
        }
        if window.is_key_down(Key::W){
            // Go forwards
            camera1.translate(&[0.0,0.0,-1.0]);
        }
        if window.is_key_down(Key::S){
            // Go backwards
            camera1.translate(&[0.0,0.0,1.0]);
        }
        if window.is_key_down(Key::Up){
            // Look up
            camera1.rotate(PI/200.0,[1.0, 0.0, 0.0]);
        }
        if window.is_key_down(Key::Down){
            // Look down
            camera1.rotate(-PI/200.0,[1.0, 0.0, 0.0]);
        }
        if window.is_key_down(Key::Left){
            // Look left
            camera1.rotate_globally(-PI/200.0,[0.0, 1.0, 0.0]);
        }
        if window.is_key_down(Key::Right){
            // Look right
            camera1.rotate_globally(PI/200.0,[0.0, 1.0, 0.0]);
        }
        if window.is_key_down(Key::LeftShift){
            // Go down
            camera1.translate(&[0.0,1.0,0.0]);
        }
        if window.is_key_down(Key::Space){
            // Go up
            camera1.translate(&[0.0,-1.0,0.0]);
        }
    }
}
