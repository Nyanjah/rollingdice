use crate::raster::Scene;

pub trait World: Scene {
    fn update(&mut self, dt: f32);
}

pub mod test_world {

    use super::World;

    use crate::body::*;
    use crate::mesh::Mesh;
    use crate::raster::*;
    use crate::transform::*;

    const FLOOR_BOUND: f32 = 1_000_000.0;

    const ANGULAR_VELOCITY: f32 = -2.5;
    const ANGULAR_ROTATION_AXIS: Vector3 = Vector3::new(1.0, 1.0, 1.0);
    const INITIAL_ANGLE_RADIANS: f32 = 2.0;

    pub struct TestWorld {
        //pub light: Vector3,
        // pub camera_transform: Transform,
        pub bodies: Vec<Body>,
        pub is_colliding: bool,
        pub meshes: Vec<Mesh>,
    }

    impl TestWorld {
        pub fn new() -> Self {
            TestWorld {
                is_colliding: false,
                //light: Vector3::new(0.0, -1.0, 0.0).unit(),
                // camera_transform: Transform::new(
                //     Vector3 {
                //         x: 0.0,
                //         y: 10.0 + 256.0,
                //         z: 0.0 - 256.0/1.0_f32.to_radians().tan(),
                //     },
                //     Quaternion::from_axis_angle(
                //         Vector3 {
                //             x: 1.0,
                //             y: 0.0,
                //             z: 0.0,
                //         },
                //         1.0_f32.to_radians(),
                //     ),
                // ),
                bodies: vec![
                    Body {
                        transform: Transform::new(
                            Vector3::new(0.0, 0.0, 75.0),
                            Quaternion::from_axis_angle(
                                ANGULAR_ROTATION_AXIS,
                                INITIAL_ANGLE_RADIANS,
                            ),
                        ),
                        half_size: Vector3::new(10.0, 10.0, 10.0),
                        ..Default::default()
                    },
                    Body {
                        transform: Transform::new(
                            Vector3 {
                                x: 25.0,
                                y: 0.0,
                                z: 75.0,
                            },
                            Quaternion::from_axis_angle(
                                ANGULAR_ROTATION_AXIS,
                                INITIAL_ANGLE_RADIANS,
                            ),
                        ),
                        half_size: Vector3::new(10.0, 10.0, 10.0),
                        ..Default::default()
                    },
                ],

                meshes: vec![
                    Mesh::load("models/skull.obj").unwrap()
                ]
            }
        }
    }

    impl World for TestWorld {
        fn update(&mut self, dt: f32) {
            for body in self.bodies.iter_mut() {
                body.transform.rotation *=
                    Quaternion::from_axis_angle(ANGULAR_ROTATION_AXIS, ANGULAR_VELOCITY * dt);
            }

            self.is_colliding = is_colliding(&self.bodies[0], &self.bodies[1]);
        }
    }

    impl Scene for TestWorld {
        fn light(&self) -> Light {
            Light {
                position: Vector3::new(0.0, 100.0, 0.0),
                offset: 0.0,
                intensity: 2000.0,
            }
        }

        fn update_geometry<T: ViewportProjector>(&self, projector: &mut RasterProjector<'_, T>) {
            let floor_color = if self.is_colliding {
                Color::RED
            } else {
                Color::GREEN
            };

            // green "floor"
            projector.project(Triangle {
                normal: Vector3::Y_AXIS,
                points: [
                    Vector3::new(-FLOOR_BOUND, 0.0, FLOOR_BOUND - FLOOR_BOUND + 1000.0),
                    Vector3::new(FLOOR_BOUND, 0.0, FLOOR_BOUND - FLOOR_BOUND + 1000.0),
                    Vector3::new(-FLOOR_BOUND, 0.0, -FLOOR_BOUND - FLOOR_BOUND + 1000.0),
                ],
                color: floor_color,
            });
            projector.project(Triangle {
                normal: Vector3::Y_AXIS,
                points: [
                    Vector3::new(-FLOOR_BOUND, 0.0, -FLOOR_BOUND - FLOOR_BOUND + 1000.0),
                    Vector3::new(FLOOR_BOUND, 0.0, FLOOR_BOUND - FLOOR_BOUND + 1000.0),
                    Vector3::new(FLOOR_BOUND, 0.0, -FLOOR_BOUND - FLOOR_BOUND + 1000.0),
                ],
                color: floor_color,
            });

            for body in self.bodies.iter() {
                body.geometry().for_each(|triangle| {
                    projector.project(triangle);
                });
            }

            for mesh in self.meshes.iter() {
                for tri in mesh.geometry.iter().cloned() {
                    projector.project(tri);
                }
            }

            // bodies
            // buf.extend(
            //     self.bodies.iter().flat_map(|body| {
            //         body.geometry().map(|mut geometry| {
            //             // let brightness = (0.5*(1.0 - geometry.normal.dot(self.light))).clamp(0.0, 1.0);
            //             // geometry.color =
            //             //  Color::from_rgb(
            //             //     (255.0*brightness).round() as u8,
            //             //     (255.0*brightness).round() as u8,
            //             //     (255.0*brightness).round() as u8
            //             // );
            //             geometry
            //         })
            //     })
            // );
        }
        //fn update_geometry(&self, consume: fn(Triangle)) {
        //buf.clear();

        //}
    }
}
