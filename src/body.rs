use crate::{
    raster::{Color, Triangle},
    transform::*,
};

#[derive(Default)]
pub struct Body {
    pub transform: Transform,
    pub half_size: Vector3,
    pub linear_velocity: Vector3,
    pub angular_velocity: Vector3,
}

impl Body {
    pub fn geometry(&self) -> impl Iterator<Item = Triangle> + '_ {
        BODY_UNIT_GEOMETRY.iter().map(|geometry| {
            let normal = self
                .transform
                .rotation
                .vector_to_world_space(geometry.normal);
            Triangle {
                normal,
                points: geometry
                    .points
                    .map(|point| self.transform.point_to_world_space(point * self.half_size)),
                color: geometry.color,
            }
        })
    }
}

// SAT collision test by checking for separating planes in all 15 axes
pub fn is_colliding(box1: &Body, box2: &Body) -> bool {
    let r_pos = box2.transform.position - box1.transform.position;

    let (box1_x_axis, box1_y_axis, box1_z_axis) = box1.transform.rotation.basis_vectors();
    let (box2_x_axis, box2_y_axis, box2_z_axis) = box2.transform.rotation.basis_vectors();

    !(get_separating_plane(r_pos, box1_x_axis, box1, box2)
        || get_separating_plane(r_pos, box1_y_axis, box1, box2)
        || get_separating_plane(r_pos, box1_z_axis, box1, box2)
        || get_separating_plane(r_pos, box2_x_axis, box1, box2)
        || get_separating_plane(r_pos, box2_y_axis, box1, box2)
        || get_separating_plane(r_pos, box2_z_axis, box1, box2)
        || get_separating_plane(r_pos, box1_x_axis.cross(box2_x_axis), box1, box2)
        || get_separating_plane(r_pos, box1_x_axis.cross(box2_y_axis), box1, box2)
        || get_separating_plane(r_pos, box1_x_axis.cross(box2_z_axis), box1, box2)
        || get_separating_plane(r_pos, box1_y_axis.cross(box2_x_axis), box1, box2)
        || get_separating_plane(r_pos, box1_y_axis.cross(box2_y_axis), box1, box2)
        || get_separating_plane(r_pos, box1_y_axis.cross(box2_z_axis), box1, box2)
        || get_separating_plane(r_pos, box1_z_axis.cross(box2_x_axis), box1, box2)
        || get_separating_plane(r_pos, box1_z_axis.cross(box2_y_axis), box1, box2)
        || get_separating_plane(r_pos, box1_z_axis.cross(box2_z_axis), box1, box2))
}

fn get_separating_plane(r_pos: Vector3, plane: Vector3, box1: &Body, box2: &Body) -> bool {
    let (box1_x_axis, box1_y_axis, box1_z_axis) = box1.transform.rotation.basis_vectors();
    let (box2_x_axis, box2_y_axis, box2_z_axis) = box2.transform.rotation.basis_vectors();

    return r_pos.dot(plane).abs()
        > (box1_x_axis * box1.half_size.x).dot(plane).abs()
            + (box1_y_axis * box1.half_size.y).dot(plane).abs()
            + (box1_z_axis * box1.half_size.z).dot(plane).abs()
            + (box2_x_axis * box2.half_size.x).dot(plane).abs()
            + (box2_y_axis * box2.half_size.y).dot(plane).abs()
            + (box2_z_axis * box2.half_size.z).dot(plane).abs();
}

pub const BODY_UNIT_GEOMETRY: [Triangle; 12] = [
    // Front 0
    Triangle {
        normal: Vector3 {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        },
        points: [
            Vector3 {
                x: 1.000,
                y: 1.000,
                z: 1.0,
            },
            Vector3 {
                x: -1.000,
                y: 1.000,
                z: 1.0,
            },
            Vector3 {
                x: -1.000,
                y: -1.000,
                z: 1.0,
            },
        ],
        color: Color::BLUE,
    },
    Triangle {
        // Front 1
        normal: Vector3 {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        },
        points: [
            Vector3 {
                x: 1.000,
                y: 1.000,
                z: 1.0,
            },
            Vector3 {
                x: -1.000,
                y: -1.000,
                z: 1.0,
            },
            Vector3 {
                x: 1.000,
                y: -1.000,
                z: 1.0,
            },
        ],
        color: Color::BLUE,
    },
    Triangle {
        // Back 0
        normal: Vector3 {
            x: 0.0,
            y: 0.0,
            z: -1.0,
        },
        points: [
            Vector3 {
                x: -1.000,
                y: 1.000,
                z: -1.0,
            },
            Vector3 {
                x: 1.000,
                y: 1.000,
                z: -1.0,
            },
            Vector3 {
                x: -1.000,
                y: -1.000,
                z: -1.0,
            },
        ],
        color: Color::RED,
    },
    Triangle {
        // Back 1
        normal: Vector3 {
            x: 0.0,
            y: 0.0,
            z: -1.0,
        },
        points: [
            Vector3 {
                x: -1.000,
                y: -1.000,
                z: -1.0,
            },
            Vector3 {
                x: 1.000,
                y: 1.000,
                z: -1.0,
            },
            Vector3 {
                x: 1.000,
                y: -1.000,
                z: -1.0,
            },
        ],
        color: Color::RED,
    },
    Triangle {
        // Right 0
        normal: Vector3 {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        },
        points: [
            Vector3 {
                x: 1.0,
                y: 1.000,
                z: 1.000,
            },
            Vector3 {
                x: 1.0,
                y: -1.000,
                z: 1.000,
            },
            Vector3 {
                x: 1.0,
                y: -1.000,
                z: -1.000,
            },
        ],
        color: Color::GREEN,
    },
    Triangle {
        // Right 1
        normal: Vector3 {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        },
        points: [
            Vector3 {
                x: 1.0,
                y: 1.000,
                z: 1.000,
            },
            Vector3 {
                x: 1.0,
                y: -1.000,
                z: -1.000,
            },
            Vector3 {
                x: 1.0,
                y: 1.000,
                z: -1.000,
            },
        ],
        color: Color::GREEN,
    },
    Triangle {
        // Left 0
        normal: Vector3 {
            x: -1.0,
            y: 0.0,
            z: 0.0,
        },
        points: [
            Vector3 {
                x: -1.0,
                y: -1.000,
                z: 1.000,
            },
            Vector3 {
                x: -1.0,
                y: 1.000,
                z: 1.000,
            },
            Vector3 {
                x: -1.0,
                y: -1.000,
                z: -1.000,
            },
        ],
        color: Color::from_rgb(255, 000, 255),
    },
    Triangle {
        // Left 1
        normal: Vector3 {
            x: -1.0,
            y: 0.0,
            z: 0.0,
        },
        points: [
            Vector3 {
                x: -1.0,
                y: -1.000,
                z: -1.000,
            },
            Vector3 {
                x: -1.0,
                y: 1.000,
                z: 1.000,
            },
            Vector3 {
                x: -1.0,
                y: 1.000,
                z: -1.000,
            },
        ],
        color: Color::from_rgb(255, 000, 255),
    },
    Triangle {
        // Up 0
        normal: Vector3 {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        },
        points: [
            Vector3 {
                x: -1.000,
                y: 1.0,
                z: 1.000,
            },
            Vector3 {
                x: 1.000,
                y: 1.0,
                z: 1.000,
            },
            Vector3 {
                x: -1.000,
                y: 1.0,
                z: -1.000,
            },
        ],
        color: Color::WHITE,
    },
    Triangle {
        // Up 1
        normal: Vector3 {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        },
        points: [
            Vector3 {
                x: -1.000,
                y: 1.0,
                z: -1.000,
            },
            Vector3 {
                x: 1.000,
                y: 1.0,
                z: 1.000,
            },
            Vector3 {
                x: 1.000,
                y: 1.0,
                z: -1.000,
            },
        ],
        color: Color::WHITE,
    },
    Triangle {
        // Down 0
        normal: Vector3 {
            x: 0.0,
            y: -1.0,
            z: 0.0,
        },
        points: [
            Vector3 {
                x: 1.000,
                y: -1.0,
                z: 1.000,
            },
            Vector3 {
                x: -1.000,
                y: -1.0,
                z: 1.000,
            },
            Vector3 {
                x: -1.000,
                y: -1.0,
                z: -1.000,
            },
        ],
        color: Color::from_rgb(000, 255, 255),
    },
    Triangle {
        // Down 1
        normal: Vector3 {
            x: 0.0,
            y: -1.0,
            z: 0.0,
        },
        points: [
            Vector3 {
                x: 1.000,
                y: -1.0,
                z: 1.000,
            },
            Vector3 {
                x: -1.000,
                y: -1.0,
                z: -1.000,
            },
            Vector3 {
                x: 1.000,
                y: -1.0,
                z: -1.000,
            },
        ],
        color: Color::from_rgb(000, 255, 255),
    },
];
