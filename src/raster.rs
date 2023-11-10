use crate::transform::*;

pub trait Scene {
    fn light(&self) -> Light;
    fn update_geometry<T: ViewportProjector>(&self, projector: &mut RasterProjector<'_, T>);
}

pub struct Light {
    pub position: Vector3,
    pub offset: f32,
    pub intensity: f32,
}

// TODO: Decouple color from triangle in scene (not even sure what the best way to architect computing color)
#[derive(Debug, Clone, Copy)]
pub struct Triangle {
    pub points: [Vector3; 3],
    pub normal: Vector3,
    pub color: Color,
}

#[derive(Debug)]
pub struct TriangleProjection {
    pub points: [Vector3; 3],
}

pub trait Viewport: Transformable  {

    type Projector: ViewportProjector;

    fn projector(&self, screen_width: u16, screen_height: u16) -> Self::Projector;
}

pub trait ViewportProjector {
    fn prepare_z_compute(&mut self, tri: &Triangle, tri_proj: &TriangleProjection);
    fn compute_z(&self, x: f32) -> f32;
    fn set_y(&mut self, y: f32);
    fn project(&self, tri: Triangle, geometry: &mut Vec<(Triangle, TriangleProjection)>);
    fn screen_space_to_world_space(&self, point: Vector3) -> Vector3;
}

pub struct RasterProjector<'a, T: ViewportProjector> {
    transform: &'a Transform,
    geometry: &'a mut Vec<(Triangle, TriangleProjection)>,
    projector: &'a T,
}

impl<'a, T: ViewportProjector> RasterProjector<'a, T> {
    // TODO: ensure winding order is correct
    pub fn project(&mut self, mut tri: Triangle) {
        tri.points = tri.points.map(|point| self.transform.point_to_local_space(point));
        tri.normal = self.transform.rotation.vector_to_local_space(tri.normal);
        self.projector.project(tri, self.geometry);
    }
}

pub trait Rasterable {
    fn rasterize(
        &mut self,
        scene: &impl Scene,
        viewport: &impl Viewport,
        screen_buffer: &mut Buffer2D<Color>,
        screen_width: u16,
        screen_height: u16,
    );

    fn antialias(self, aliasing: u8) -> Antialias<Self> where Self: Sized {
        Antialias {
            raster: self,
            aliasing,
        }
    }
}

fn edge(e0: Vector2, e1: Vector2, p: Vector2) -> f32 {
    -((p.x - e0.x) * (e1.y - e0.y) - (p.y - e0.y) * (e1.x - e0.x))
}

#[derive(Debug, Default)]
pub struct Raster {
    z_buffer: Buffer2D<f32>,
    geometry_buffer: Vec<(Triangle, TriangleProjection)>,
}

impl Rasterable for Raster {
    fn rasterize(
        &mut self,
        scene: &impl Scene,
        viewport: &impl Viewport,
        screen_buffer: &mut Buffer2D<Color>,
        screen_width: u16,
        screen_height: u16,
    ) {
        screen_buffer.clear_and_resize(screen_width as usize, screen_height as usize, Color::default());

        if screen_width == 0 || screen_height == 0 {
            return;
        }

        self.geometry_buffer.clear();

        self.z_buffer
            .clear_and_resize(screen_width as usize, screen_height as usize, f32::INFINITY);

        let mut projector = viewport.projector(screen_width, screen_height);

        scene.update_geometry(&mut RasterProjector {
            transform: viewport.transform(),
            geometry: &mut self.geometry_buffer,
            projector: &projector,
        });

        let light = scene.light();
        let light_position = viewport.transform().point_to_local_space(light.position);
        let light_intensity_recip = light.intensity.recip();

        for (tri, tri_proj) in self.geometry_buffer.iter() {
            // Backface culling
            if 
                (tri_proj.points[0].x - tri_proj.points[1].x) * (tri_proj.points[0].y - tri_proj.points[2].y) - 
                (tri_proj.points[0].y - tri_proj.points[1].y) * (tri_proj.points[0].x - tri_proj.points[2].x) < 0.0 
            {
                continue;
            }

            // Cull triangles behind the camera
            if tri_proj.points[0].z.max(tri_proj.points[1].z).max(tri_proj.points[2].z) <= 0.0 {
                continue;
            }

            // TODO: support UV texture
            let tri_color_r = tri.color.r as f32;
            let tri_color_g = tri.color.g as f32;
            let tri_color_b = tri.color.b as f32;

            let tri_screen_float_points = tri_proj
                .points
                .map(|point| Vector2 {
                    x: point.x * screen_width as f32,
                    y: point.y * screen_height as f32,
                });

            let tri_bounding_x0 = tri_screen_float_points[0].x.min(tri_screen_float_points[1].x).min(tri_screen_float_points[2].x).floor().max(0.0) as usize;
            let tri_bounding_y0 = tri_screen_float_points[0].y.min(tri_screen_float_points[1].y).min(tri_screen_float_points[2].y).floor().max(0.0) as usize;
            let tri_bounding_x1 = tri_screen_float_points[0].x.max(tri_screen_float_points[1].x).max(tri_screen_float_points[2].x).ceil().min(screen_width as f32) as usize;
            let tri_bounding_y1 = tri_screen_float_points[0].y.max(tri_screen_float_points[1].y).max(tri_screen_float_points[2].y).ceil().min(screen_height as f32) as usize;

            projector.prepare_z_compute(tri, tri_proj);

            for screen_y in tri_bounding_y0..tri_bounding_y1 {
                let screen_float_y = screen_y as f32 + 0.5;
                let viewport_y = screen_float_y / (screen_height as f32);
                projector.set_y(viewport_y);

                let z_buffer_slice = self.z_buffer.get_row_mut(screen_y);
                let pixel_buffer_slice = screen_buffer.get_row_mut(screen_y);

                for screen_x in tri_bounding_x0..tri_bounding_x1 {
                    let screen_float_x = screen_x as f32 + 0.5;
                    let pixel_screen_float_point = Vector2 {
                        x: screen_float_x,
                        y: screen_float_y,
                    };

                    let w0 = edge(tri_screen_float_points[0], tri_screen_float_points[1], pixel_screen_float_point);
                    let w1 = edge(tri_screen_float_points[1], tri_screen_float_points[2], pixel_screen_float_point);
                    let w2 = edge(tri_screen_float_points[2], tri_screen_float_points[0], pixel_screen_float_point);

                    if !(w0 >= 0.0 && w1 >= 0.0 && w2 >= 0.0) {
                        continue;
                    }

                    let viewport_x = screen_float_x / (screen_width as f32);
                    let viewport_z = projector.compute_z(viewport_x);
                    if viewport_z < 0.0 || viewport_z > z_buffer_slice[screen_x] {
                        continue;
                    }

                    z_buffer_slice[screen_x] = viewport_z;

                    let pixel_point = light_position - projector.screen_space_to_world_space(Vector3::new(
                        viewport_x,
                        viewport_y,
                        viewport_z,
                    ));

                    // shading
                    
                    // TODO: optimize this

                    let pixel_point_distance = pixel_point.magnitude();

                    let light_direction_factor = tri.normal.dot(pixel_point) / pixel_point_distance;
                    let light_intensity_factor = ((pixel_point_distance - light.offset).max(0.0) * light_intensity_recip + 1.0).powi(2).recip();
                    let pixel_brightness = light_direction_factor * light_intensity_factor;

                    pixel_buffer_slice[screen_x] = Color::new(
                        (tri_color_r * pixel_brightness) as u8,
                        (tri_color_g * pixel_brightness) as u8,
                        (tri_color_b * pixel_brightness) as u8,
                    );
                }
            }
        }
    }
}

pub struct Antialias<T: Rasterable> {
    aliasing: u8,
    raster: T,
}

impl<T: Rasterable> Rasterable for Antialias<T> {
    fn rasterize(
        &mut self,
        scene: &impl Scene,
        viewport: &impl Viewport,
        screen_buffer: &mut Buffer2D<Color>,
        screen_width: u16,
        screen_height: u16,
    ) {
        self.raster.rasterize(
            scene,
            viewport,
            screen_buffer,
            screen_width * self.aliasing as u16,
            screen_height * self.aliasing as u16,
        );
        let sample_size = self.aliasing as usize;

        screen_buffer.condense(screen_width as usize, screen_height as usize, |buf, x, y| {
            let (sr, sg, sb) = buf
                .area(sample_size * x, sample_size * y, sample_size, sample_size)
                .flatten()
                .fold((0_u16, 0_u16, 0_u16), |(sr, sg, sb), color| {
                    
                    (sr + color.r as u16, sg + color.g as u16, sb + color.b as u16)
                });
            let sample_area = (sample_size * sample_size) as u16;
            Color::new(
                (sr / sample_area) as u8,
                (sg / sample_area) as u8,
                (sb / sample_area) as u8,
            )
        })
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub struct Color {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}

impl Color {

    pub const BLACK: Color = Color::from_u32(0x000000);
    pub const WHITE: Color = Color::from_u32(0xFFFFFF);
    pub const RED: Color = Color::from_u32(0xFF0000);
    pub const GREEN: Color = Color::from_u32(0x00FF00);
    pub const BLUE: Color = Color::from_u32(0x0000FF);

    pub const fn new(r: u8, g: u8, b: u8) -> Color {
        Color { r, g, b }
    }

    pub const fn from_u32(color: u32) -> Color {
        Color {
            r: ((color >> 16) & 0xFF) as u8,
            g: ((color >> 8) & 0xFF) as u8,
            b: (color & 0xFF) as u8,
        }
    }

    pub fn random() -> Color {
        Color {
            r: rand::random(),
            g: rand::random(),
            b: rand::random(),
        }
    }

    pub const fn u32(self) -> u32 {
        (self.r as u32) << 16 | (self.g as u32) << 8 | self.b as u32
    }
}

#[derive(Debug, Default)]
pub struct Buffer2D<T> {
    pub width: usize,
    pub height: usize,
    pub data: Vec<T>,
}

impl<T> Buffer2D<T> {
    pub fn clear_and_resize(&mut self, width: usize, height: usize, default: T)
    where
        T: Copy,
    {
        let current_len = self.data.len();
        let desired_len = width * height;
        self.data.resize(desired_len, default);
        self.data[..current_len.min(desired_len)].fill(default);

        self.width = width;
        self.height = height;
    }

    pub fn get(&self, x: usize, y: usize) -> &T {
        &self.data[y * self.width + x]
    }

    pub fn get_mut(&mut self, x: usize, y: usize) -> &mut T {
        &mut self.data[y * self.width + x]
    }

    pub fn get_row(&self, y: usize) -> &[T] {
        &self.data[y * self.width..(y + 1) * self.width]
    }

    pub fn get_row_mut(&mut self, y: usize) -> &mut [T] {
        &mut self.data[y * self.width..(y + 1) * self.width]
    }

    pub fn area(
        &self,
        x: usize,
        y: usize,
        width: usize,
        height: usize,
    ) -> impl Iterator<Item = impl Iterator<Item = &T>> {
        (y..y + height).map(move |y| self.get_row(y)[x..=x + width - 1].iter())
    }

    pub fn condense<F>(&mut self, width: usize, height: usize, f: F)
    where
        F: Fn(&Buffer2D<T>, usize, usize) -> T,
    {
        for y in 0..height {
            for x in 0..width {
                self.data[y * width + x] = f(&self, x, y);
            }
        }

        self.width = width;
        self.height = height;
        self.data.truncate(width * height);
    }
}
