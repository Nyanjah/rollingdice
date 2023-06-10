use crate::{
    transform::*, 
    raster::{Triangle, TriangleProjection, Viewport, ViewportProjector}
};

#[derive(Clone, Default)]
pub struct OrthographicCamera {
    pub transform: Transform
}

impl Viewport for OrthographicCamera {

    type Projector = OrthographicProjector;

    fn projector(&self, screen_width: u16, screen_height: u16) -> Self::Projector {
        OrthographicProjector {
            size: Vector2::new(screen_width as f64, screen_height as f64),
            ..Default::default()
        }
    }
}

impl Transformable for OrthographicCamera {
    fn transform(&self) -> &Transform {
        &self.transform
    }

    fn transform_mut(&mut self) -> &mut Transform {
        &mut self.transform
    }
}

#[derive(Default)]
pub struct  OrthographicProjector {
    size: Vector2,
    z_origin: Vector3,
    z_change: Vector3,
    z_y: f64,
}

impl ViewportProjector for OrthographicProjector {
    fn project(&self, tri: Triangle, geometry: &mut Vec<(Triangle, TriangleProjection)>) {
        let tri_proj = TriangleProjection { 
            points: tri.points.map(|p| {
                Vector3::new(0.5 + p.x / self.size.x, 0.5 - p.y / self.size.y, p.z)
            }) 
        };

        geometry.push((tri, tri_proj));
    }

    fn prepare_z_compute(&mut self, scene_tri: &Triangle, tri_proj: &TriangleProjection) {
        // Calculate the distance the plane containing the triangle recedes from the camera plane
        // when traversing the camera's width and height
        self.z_origin = tri_proj.points[0];
        self.z_change = Vector3 {
            x: -scene_tri.normal.x * self.size.x / scene_tri.normal.z,
            y: scene_tri.normal.y * self.size.y / scene_tri.normal.z,
            z: 0.0,
        };
    }

    fn set_y(&mut self, y: f64) {
        self.z_y = self.z_origin.z + (y - self.z_origin.y) * self.z_change.y
    }

    fn compute_z(&self, x: f64) -> f64 {
        self.z_y + (x - self.z_origin.x) * self.z_change.x
    }
}

pub struct PerspectiveCamera {
    pub transform: Transform
}

impl Transformable for PerspectiveCamera {
    fn transform(&self) -> &Transform {
        &self.transform
    }

    fn transform_mut(&mut self) -> &mut Transform {
        &mut self.transform
    }
}

impl Viewport for PerspectiveCamera {

    type Projector = PerspectiveProjector;

    fn projector(&self, screen_width: u16, screen_height: u16) -> Self::Projector {

        let vertical_tan_half_fov = 35.0f64.to_radians().tan();

        PerspectiveProjector {
            tan_half_fov: Vector2::new((screen_width as f64 / screen_height as f64) * vertical_tan_half_fov, vertical_tan_half_fov),
            z_cutoff: 0.1,
            ..Default::default()
        }
    }
}

#[derive(Default)]
pub struct PerspectiveProjector {
    tan_half_fov: Vector2,
    z_cutoff: f64,
    z_plane_origin: Vector3,
    z_plane_normal: Vector3,
    y: f64,
}

impl PerspectiveProjector {
    fn project(&self, tri: Triangle, geometry: &mut Vec<(Triangle, TriangleProjection)>) {
        let tri_proj = TriangleProjection { 
            points: tri.points.map(|point| {
                Vector3::new(
                    0.5 * (1.0 + point.x / (self.tan_half_fov.x * point.z.abs())), 
                    0.5 * (1.0 - point.y / (self.tan_half_fov.y * point.z.abs())), 
                    point.z - self.z_cutoff
                )
            })
        };

        geometry.push((tri, tri_proj));
    }
}

impl ViewportProjector for PerspectiveProjector {
    fn project(&self, tri: Triangle, geometry: &mut Vec<(Triangle, TriangleProjection)>) {

        let clip_plane_z = self.z_cutoff;

        for i in 0..3 {

            let p0 = tri.points[i];
            let p1 = tri.points[(i + 1) % 3];

            if p0.z > clip_plane_z && p1.z < clip_plane_z {

                let p2 = tri.points[(i + 2) % 3];
                
                if p2.z < clip_plane_z {

                    let Triangle { normal, color, .. } = tri;

                    let v01 = p1 - p0;
                    let v02 = p2 - p0;
                    
                    self.project(
                        Triangle { 
                            // TODO: Enforce winding order
                            points: [
                                p0, 
                                p0.lerp(p1, (clip_plane_z - p0.z) / (v01.unit().z * v01.magnitude())),
                                p0.lerp(p2, (clip_plane_z - p0.z) / (v02.unit().z * v02.magnitude())),
                            ], 
                            normal, color 
                        }, 
                        geometry
                    );
                } else {
                    let v01 = p1 - p0;
                    let v21 = p1 - p2;

                    let p0z = p0.lerp(p1, (clip_plane_z - p0.z) / (v01.unit().z * v01.magnitude()));
                    let p2z = p2.lerp(p1, (clip_plane_z - p2.z) / (v21.unit().z * v21.magnitude()));

                    let Triangle { normal, color, .. } = tri;

                    // TODO: Enforce winding order
                    self.project(Triangle { points: [p0z, p0, p2], normal, color }, geometry);
                    self.project(Triangle { points: [p0z, p2z, p2], normal, color }, geometry);
                }

                return;
            }
        }

        self.project(tri, geometry);
    }

    fn prepare_z_compute(&mut self, scene_tri: &Triangle, _: &TriangleProjection) {
        self.z_plane_origin = scene_tri.points[0];
        self.z_plane_normal = scene_tri.normal;
    }

    fn set_y(&mut self, y: f64) {
        self.y = y;
    }

    fn compute_z(&self, x: f64) -> f64 {
        let y = self.y;
        let point = Vector3::new(
            self.z_cutoff * self.tan_half_fov.x * (2.0 * x - 1.0), 
            self.z_cutoff * self.tan_half_fov.y * (1.0 - 2.0 * y), 
            self.z_cutoff
        );

        self.z_plane_origin.dot(self.z_plane_normal) / point.dot(self.z_plane_normal)
    }
}

