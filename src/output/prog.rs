use std::time::Instant;

use minifb::*;

use crate::world::*;
use crate::raster::*;
use crate::transform::*;
use crate::camera::*;

const WIDTH: u16 = 1080;
const HEIGHT: u16 = 720;
const SECONDS_PER_FRAME: f32 = 1.0/60.0; // MAX 60 FPS
const CAMERA_LINEAR_MIN_SPEED: f64 = 50.0;
const CAMERA_LINEAR_MAX_SPEED: f64 = 200.0;
const CAMERA_LINEAR_SPEED_TRANSITION_TIME: f64 = 5.0;
const CAMERA_INPUT_DIRECTION_MAP: [(Key, Vector3); 6] = [
    (Key::W, Vector3::new( 0.0,  0.0,  1.0)),
    (Key::A, Vector3::new(-1.0,  0.0,  0.0)),
    (Key::S, Vector3::new( 0.0,  0.0, -1.0)),
    (Key::D, Vector3::new( 1.0,  0.0,  0.0)),
    (Key::Q, Vector3::new( 0.0, -1.0,  0.0)),
    (Key::E, Vector3::new( 0.0,  1.0,  0.0))
];
//const CMAERA_BASE_ANGULAR_SPEED: f64 = 

//const CAMERA_LINEAR_SPEED: f64 = 10.0;

fn ease_in_cubic(alpha: f64) -> f64 {
    alpha.powi(3)
}

pub fn run() {
    
    // Setting up the window
    let mut window = Window::new(
        "Rasterization",
        WIDTH as usize,
        HEIGHT as usize,
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
    ).unwrap();

    let mut test_world = test_world::TestWorld::new();
    let mut buffer = Buffer2D::<Color>::default();
    let mut raster = Raster::default();
    let mut camera = PerspectiveCamera {
        transform: Transform::new(
            Vector3 {
                x: 0.0,
                y: 30.0,
                z: 30.0,
            },
            Quaternion::from_axis_angle(
                Vector3 {
                    x: 1.0,
                    y: 0.0,
                    z: 0.0,
                },
                20.0_f64.to_radians(),
            ),
        ), // ..Default::default()
    };

    // Creating an empty window buffer for minifb to update the window with
    let mut window_buffer: Vec<u32> = Vec::with_capacity(WIDTH as usize * HEIGHT as usize);
    
    // (Optional) Limit the window update rate to control CPU usage
    window.limit_update_rate(Some(std::time::Duration::from_secs_f32(SECONDS_PER_FRAME)));

    let mut then = Instant::now();

    let mut translating_camera_start = None;

    while window.is_open() && !window.is_key_down(Key::Escape) {

        let now = Instant::now();
        let dt = now.duration_since(then).as_secs_f64();
        test_world.update(dt);
        then = now;
        
        let mut camera_translation_dir = Vector3::ZERO;

        let mut translating_camera = false;

        for (key, dir) in CAMERA_INPUT_DIRECTION_MAP {
            if window.is_key_down(key) {
                if translating_camera_start.is_none() {
                    translating_camera_start = Some(now);
                }
                translating_camera = true;
                camera_translation_dir += dir;
            }
        }

        if !translating_camera {
            translating_camera_start = None;
        }

        if let Some(translating_camera_start) = translating_camera_start {
            if camera_translation_dir != Vector3::ZERO {
                let translation_duration = now.saturating_duration_since(translating_camera_start).as_secs_f64();
                let speed_alpha = ease_in_cubic((translation_duration/CAMERA_LINEAR_SPEED_TRANSITION_TIME).min(1.0));
                let speed = CAMERA_LINEAR_MAX_SPEED*speed_alpha + CAMERA_LINEAR_MIN_SPEED*(1.0 - speed_alpha);
                let translation = camera.transform().rotation.vector_to_world_space(camera_translation_dir.unit()*speed*dt);
                camera.transform_mut().position += translation;
            }
        }

        // Rasterize the world
        raster.rasterize(&test_world, &camera, &mut buffer, WIDTH, HEIGHT);

        // Update the window with the prepared frame
        window_buffer.clear();
        window_buffer.extend(buffer.data.iter().map(|pixel| pixel.u32()));
        window.update_with_buffer(&window_buffer, WIDTH as usize, HEIGHT as usize).unwrap();
        // if window.is_key_down(Key::D){
        //     // Go left
        //     camera1.translate(&[-1.0,0.0,0.0]);
        // }
        // if window.is_key_down(Key::A){
        //     // Go right
        //     camera1.translate(&[1.0,0.0,0.0]);
        // }
        // if window.is_key_down(Key::W){
        //     // Go forwards
        //     camera1.translate(&[0.0,0.0,-1.0]);
        // }
        // if window.is_key_down(Key::S){
        //     // Go backwards
        //     camera1.translate(&[0.0,0.0,1.0]);
        // }
        // if window.is_key_down(Key::Up){
        //     // Look up
        //     camera1.rotate(PI/200.0,[1.0, 0.0, 0.0]);
        // }
        // if window.is_key_down(Key::Down){
        //     // Look down
        //     camera1.rotate(-PI/200.0,[1.0, 0.0, 0.0]);
        // }
        // if window.is_key_down(Key::Left){
        //     // Look left
        //     camera1.rotate_globally(-PI/200.0,[0.0, 1.0, 0.0]);
        // }
        // if window.is_key_down(Key::Right){
        //     // Look right
        //     camera1.rotate_globally(PI/200.0,[0.0, 1.0, 0.0]);
        // }
        // if window.is_key_down(Key::LeftShift){
        //     // Go down
        //     camera1.translate(&[0.0,1.0,0.0]);
        // }
        // if window.is_key_down(Key::Space){
        //     // Go up
        //     camera1.translate(&[0.0,-1.0,0.0]);
        // }
    }
}
