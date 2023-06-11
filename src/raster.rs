use crate::transform::*;

use std::ops::RangeInclusive;

pub trait Scene {
    fn update_geometry<T: ViewportProjector>(&self, projector: &mut RasterProjector<'_, T>);
}

// TODO: Decouple color from triangle in scene (not even sure what the best way to architect computing color)
#[derive(Debug)]
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
    fn compute_z(&self, x: f64) -> f64;
    fn set_y(&mut self, y: f64);
    fn project(&self, tri: Triangle, geometry: &mut Vec<(Triangle, TriangleProjection)>);
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

#[derive(Debug, Default)]
pub struct Raster {
    z_buffer: Buffer2D<f64>,
    geometry_buffer: Vec<(Triangle, TriangleProjection)>,
    horizonal_line_buffer: Vec<(i64, i64)>,
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
            .clear_and_resize(screen_width as usize, screen_height as usize, f64::INFINITY);

        self.horizonal_line_buffer
            .resize(screen_height as usize, (i64::MAX, i64::MAX));

        let mut projector = viewport.projector(screen_width, screen_height);

        scene.update_geometry(&mut RasterProjector {
            transform: viewport.transform(),
            geometry: &mut self.geometry_buffer,
            projector: &projector,
        });

        for (tri, tri_proj) in self.geometry_buffer.iter() {
            // Backface culling
            if 
                (tri_proj.points[0].x - tri_proj.points[1].x) * (tri_proj.points[0].y - tri_proj.points[2].y) - 
                (tri_proj.points[0].y - tri_proj.points[1].y) * (tri_proj.points[0].x - tri_proj.points[2].x) < 0.0 
            {
                continue;
            }

            // Cull triangles that are completely off screen
            {
                let x0 = tri_proj
                    .points
                    .iter()
                    .fold(f64::MAX, |x0, point| x0.min(point.x));
                let y0 = tri_proj
                    .points
                    .iter()
                    .fold(f64::MAX, |y0, point| y0.min(point.y));
                let x1 = tri_proj
                    .points
                    .iter()
                    .fold(f64::MIN, |x1, point| x1.max(point.x));
                let y1 = tri_proj
                    .points
                    .iter()
                    .fold(f64::MIN, |y1, point| y1.max(point.y));
                let z1 = tri_proj
                    .points
                    .iter()
                    .fold(f64::MIN, |z1, point| z1.max(point.z));

                if !(z1 > 0.0 && x0 < 1.0 && x1 > 0.0 && y0 < 1.0 && y1 > 0.0) {
                    continue;
                }
            }

            let screen_points = tri_proj
                .points
                .map(|point| {
                    (
                        (point.x * screen_width as f64).round() as i64,
                        (point.y * screen_height as f64).round() as i64,
                    )
                });

            for line in [
                [&screen_points[0], &screen_points[1]],
                [&screen_points[1], &screen_points[2]],
                [&screen_points[2], &screen_points[0]],
            ] {
                let &(x1, y1) = line[0];
                let &(x2, y2) = line[1];

                let m = (y2 - y1) as f64 / (x2 - x1) as f64;

                for y in y1.min(y2).max(0)..=y1.max(y2).min(screen_height as i64 - 1) {
                    let x = ((y - y1) as f64 / m).round() as i64 + x1;
                    let row = &mut self.horizonal_line_buffer[y as usize];
                    if row.0 == i64::MAX {
                        *row = (x, x)
                    } else if x < row.0 {
                        row.0 = x
                    } else if x > row.1 {
                        row.1 = x
                    }
                }

                // // if dx > dy {
                // //     self.bressenham_line(x1, y1, x2, y2, dx, dy, 0);
                // // } else {
                // //     self.bressenham_line(y1, x1, y2, x2, dy, dx, 1);
                // // }
            }

            let y_min = screen_points
                .iter()
                .fold(i64::MAX, |min, &(_, y)| min.min(y))
                .clamp(0, screen_height.saturating_sub(1) as i64) as usize;
            let y_max = screen_points
                .iter()
                .fold(i64::MIN, |max, &(_, y)| max.max(y))
                .clamp(0, screen_height.saturating_sub(1) as i64) as usize;

            projector.prepare_z_compute(tri, tri_proj);

            for y in y_min..=y_max {
                let (x_min, x_max) = {
                    let (min, max) = self.horizonal_line_buffer[y];
                    if min < 0 && max < 0 || min >= screen_width as i64 && max >= screen_width as i64
                    {
                        continue;
                    }
                    (
                        min.clamp(0, screen_width as i64 - 1) as usize,
                        max.clamp(0, screen_width as i64 - 1) as usize,
                    )
                };

                projector.set_y((y as f64 + 0.5) / (screen_height as f64));

                let z_slice = self.z_buffer.get_range_mut(y, x_min..=x_max);
                let pixel_slice = screen_buffer.get_range_mut(y, x_min..=x_max);

                for (i, (last_z, pixel)) in
                    z_slice.iter_mut().zip(pixel_slice.iter_mut()).enumerate()
                {
                    let z = projector.compute_z(((x_min + i) as f64 + 0.5) / (screen_width as f64));
                    if z > *last_z {
                        continue;
                    }

                    *last_z = z;
                    *pixel = tri.color;
                }
            }

            self.horizonal_line_buffer[y_min..=y_max].fill((i64::MAX, i64::MAX));
        }
    }

    // fn bressenham_line(
    //     &mut self,
    //     mut x1: i64,
    //     mut y1: i64,
    //     x2: i64,
    //     y2: i64,
    //     dx: i64,
    //     dy: i64,
    //     decide: i64,
    // ) {
    //     let mut pk = 2 * dy - dx;

    //     for _ in 0..dx {
    //         if x1 < x2 {
    //             x1 += 1;
    //         } else {
    //             x1 -= 1;
    //         }

    //         if pk < 0 {
    //             let (x1, y1) = if decide == 0 { (x1, y1) } else { (y1, x1) };
    //             if y1 >= 0 && y1 < self.screen_buffer.height as i64 {
    //                 self.insert_bressenham_point(x1, y1 as usize);
    //             }
    //             pk += 2 * dy;
    //         } else {
    //             if y1 < y2 {
    //                 y1 += 1;
    //             } else {
    //                 y1 -= 1;
    //             }
    //             let (x1, y1) = if decide == 0 { (x1, y1) } else { (y1, x1) };
    //             if y1 >= 0 && y1 < self.screen_buffer.height as i64 {
    //                 self.insert_bressenham_point(x1, y1 as usize);
    //             }
    //             pk += 2 * dy - 2 * dx;
    //         }
    //     }
    // }

    // fn insert_bressenham_point(&mut self, x1: i64, y1: usize) {
    //     let row = &mut self.horizonal_line_buffer[y1];
    //     if row.0 == i64::MAX {
    //         *row = (x1, x1)
    //     } else if x1 < row.0 {
    //         row.0 = x1
    //     } else if x1 > row.1 {
    //         row.1 = x1
    //     }
    // }
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
                    let (r, g, b) = (*color).rgb();
                    (sr + r as u16, sg + g as u16, sb + b as u16)
                });
            let sample_area = (sample_size * sample_size) as u16;
            Color::from_rgb(
                (sr / sample_area) as u8,
                (sg / sample_area) as u8,
                (sb / sample_area) as u8,
            )
        })
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub struct Color(u32);

impl Color {

    pub const BLACK: Color = Color(0x000000);
    pub const WHITE: Color = Color(0xFFFFFF);
    pub const RED: Color = Color(0xFF0000);
    pub const GREEN: Color = Color(0x00FF00);
    pub const BLUE: Color = Color(0x0000FF);

    pub const fn new(color: u32) -> Color {
        Color(color)
    }

    pub const fn from_rgb(r: u8, g: u8, b: u8) -> Color {
        Color((r as u32) << 16 | (g as u32) << 8 | b as u32)
    }

    pub const fn r(self) -> u8 {
        (self.0 >> 16) as u8
    }

    pub const fn g(self) -> u8 {
        (self.0 >> 8) as u8
    }

    pub const fn b(self) -> u8 {
        self.0 as u8
    }

    pub const fn u32(self) -> u32 {
        self.0
    }

    pub const fn rgb(self) -> (u8, u8, u8) {
        (self.r(), self.g(), self.b())
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

    pub fn get_range(&self, y: usize, r: RangeInclusive<usize>) -> &[T] {
        &self.data[y * self.width + r.start()..=y * self.width + r.end()]
    }

    pub fn get_range_mut(&mut self, y: usize, r: RangeInclusive<usize>) -> &mut [T] {
        &mut self.data[y * self.width + r.start()..=y * self.width + r.end()]
    }

    pub fn area(
        &self,
        x: usize,
        y: usize,
        width: usize,
        height: usize,
    ) -> impl Iterator<Item = impl Iterator<Item = &T>> {
        (y..y + height).map(move |y| self.get_range(y, x..=x + width - 1).iter())
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
