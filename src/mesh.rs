use std::path::Path;

use tobj::{load_obj, LoadOptions, LoadError};

use crate::{raster::{Triangle, Color}, transform::Vector3};

pub struct Mesh {
    pub geometry: Vec<Triangle>
}

impl Mesh {
    pub fn load<P: AsRef<Path> + std::fmt::Debug>(path: P) -> Result<Self, LoadError> {
        let (models, _) = load_obj(path, &LoadOptions {
            triangulate: true,
            single_index: true,
            ..Default::default()
        })?;

        let geometry = models.iter().map(|model| {
            let indices = &model.mesh.indices;
            let positions = &model.mesh.positions;

            indices.chunks_exact(3).map(|tri_indicies| {
                let &[i0, i1, i2] = tri_indicies else {
                    unreachable!("Indices are not an exact chunk of 3")
                };

                let points = [
                    Vector3::new(positions[i0 as usize * 3], positions[i0 as usize * 3 + 1], positions[i0 as usize * 3 + 2]),
                    Vector3::new(positions[i1 as usize * 3], positions[i1 as usize * 3 + 1], positions[i1 as usize * 3 + 2]),
                    Vector3::new(positions[i2 as usize * 3], positions[i2 as usize * 3 + 1], positions[i2 as usize * 3 + 2]),
                ];

                let normal = (points[1] - points[0]).cross(points[2] - points[0]).unit();

                Triangle {
                    points,
                    normal,
                    color: Color::random(),
                }
            })
        }).flatten().collect();

        Ok(Mesh { geometry })
    }
}