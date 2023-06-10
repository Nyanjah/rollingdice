use std::time::Instant;

use minifb::*;

use crate::camera::*;
use crate::raster::*;
use crate::transform::*;
use crate::world::*;

const WIDTH: u16 = 1080;
const HEIGHT: u16 = 720;
const SECONDS_PER_FRAME: f32 = 1.0 / 60.0; // MAX 60 FPS

const CAMERA_ROTATION_RADIANS_PER_PIXEL: f64 = 0.01;

const CAMERA_LINEAR_MIN_SPEED: f64 = 50.0;
const CAMERA_LINEAR_MAX_SPEED: f64 = 500.0;
const CAMERA_LINEAR_SPEED_TRANSITION_TIME: f64 = 10.0;
const CAMERA_LINEAR_SPEED_INPUT_DIRECTION_MAP: [(Key, Vector3); 6] = [
    (Key::W, Vector3::new(0.0, 0.0, 1.0)),
    (Key::A, Vector3::new(-1.0, 0.0, 0.0)),
    (Key::S, Vector3::new(0.0, 0.0, -1.0)),
    (Key::D, Vector3::new(1.0, 0.0, 0.0)),
    (Key::Q, Vector3::new(0.0, -1.0, 0.0)),
    (Key::E, Vector3::new(0.0, 1.0, 0.0)),
];

const CAMERA_LINEAR_SPEED_SLOWDOWN_INPUT: Key = Key::LeftShift;
const CAMERA_LINEAR_SPEED_SLOWDOWN_MULTIPLIER: f64 = 0.1;

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
    )
    .unwrap();

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
        ),
    };

    // Creating an empty window buffer for minifb to update the window with
    let mut window_buffer: Vec<u32> = Vec::with_capacity(WIDTH as usize * HEIGHT as usize);

    // (Optional) Limit the window update rate to control CPU usage
    window.limit_update_rate(Some(std::time::Duration::from_secs_f32(SECONDS_PER_FRAME)));

    let mut then = Instant::now();

    let mut translating_camera_start = None;
    let mut mouse_position = Some(Vector2::ZERO);

    while window.is_open() && !window.is_key_down(Key::Escape) {
        let now = Instant::now();
        let dt = now.duration_since(then).as_secs_f64();
        test_world.update(dt);
        then = now;

        // Camera rotation
        
        if window.get_mouse_down(MouseButton::Right) {
            let new_mouse_position = window
                .get_mouse_pos(MouseMode::Pass)
                .map(|(mx, my)| Vector2::new(mx as f64, my as f64));

            if let Some(mouse_position) = mouse_position {
                if let Some(new_mouse_position) = new_mouse_position {
                    let mouse_delta = new_mouse_position - mouse_position;
                    let mouse_delta_y_rotation = Quaternion::from_axis_angle(
                        Vector3::X_AXIS,
                        CAMERA_ROTATION_RADIANS_PER_PIXEL * mouse_delta.y,
                    );
                    let mouse_delta_x_rotation = Quaternion::from_axis_angle(
                        camera
                            .transform()
                            .rotation
                            .vector_to_local_space(Vector3::Y_AXIS),
                        CAMERA_ROTATION_RADIANS_PER_PIXEL * mouse_delta.x,
                    );

                    camera.transform_mut().rotation *=
                        mouse_delta_x_rotation * mouse_delta_y_rotation;
                }
            }

            mouse_position = new_mouse_position;
        } else {
            mouse_position = None;
        }

        // Camera translation

        let mut camera_translation_dir = Vector3::ZERO;

        let mut translating_camera = false;

        for (key, dir) in CAMERA_LINEAR_SPEED_INPUT_DIRECTION_MAP {
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
                let translation_duration = now
                    .saturating_duration_since(translating_camera_start)
                    .as_secs_f64();
                let speed_alpha = ease_in_cubic(
                    (translation_duration / CAMERA_LINEAR_SPEED_TRANSITION_TIME).min(1.0),
                );
                let speed = (CAMERA_LINEAR_MAX_SPEED * speed_alpha
                    + CAMERA_LINEAR_MIN_SPEED * (1.0 - speed_alpha))
                    * if window.is_key_down(CAMERA_LINEAR_SPEED_SLOWDOWN_INPUT) {
                        CAMERA_LINEAR_SPEED_SLOWDOWN_MULTIPLIER
                    } else {
                        1.0
                    };
                let translation = camera
                    .transform()
                    .rotation
                    .vector_to_world_space(camera_translation_dir.unit() * speed * dt);
                camera.transform_mut().position += translation;
            }
        }

        // Rasterize the world
        raster.rasterize(&test_world, &camera, &mut buffer, WIDTH, HEIGHT);

        // Update the window with the prepared frame
        window_buffer.clear();
        window_buffer.extend(buffer.data.iter().map(|pixel| pixel.u32()));
        window
            .update_with_buffer(&window_buffer, WIDTH as usize, HEIGHT as usize)
            .unwrap();
    }
}
