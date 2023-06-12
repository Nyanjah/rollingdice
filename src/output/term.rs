extern crate log;

use tui_logger::{TuiLoggerWidget, TuiLoggerLevelOutput};
use tui::{
    backend::CrosstermBackend, 
    layout::{Layout, Direction, Constraint}, 
    widgets::{Borders, Block}, 
    style::{Style, Color as TerminalColor}
};
use crossterm::{
    event::{poll, read, Event, KeyCode, KeyModifiers},
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};

use std::{io::stdout, time::{Duration, Instant}};

use crate::world::*;
use crate::raster::*;
use crate::transform::*;
use crate::camera::*;

pub fn run() {
    tui_logger::init_logger(log::LevelFilter::Trace).unwrap();
    tui_logger::set_default_level(log::LevelFilter::Trace);

    enable_raw_mode().expect("Failed to enable terminal raw mode");

    let mut stdout = stdout();
    execute!(stdout, EnterAlternateScreen).expect("Failed to enter alternate terminal screen");

    let mut terminal = tui::Terminal::new(CrosstermBackend::new(stdout))
        .expect("Failed to create interface to terminal backend");

    let mut test_world = test_world::TestWorld::new();
    let mut buffer = Buffer2D::<Color>::default();
    let mut raster = Raster::default().antialias(4);
    let camera = PerspectiveCamera {
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
                20.0_f32.to_radians(),
            ),
        ), // ..Default::default()
    };
    
    let mut then = Instant::now();

    loop {
        if poll(Duration::from_micros(16666)).unwrap() {
            match read().unwrap() {
                Event::Key(key_event) => {
                    // Esc or Crtl+C interrupt handler
                    if key_event.code == KeyCode::Esc // Esc is an exit if debugger isnt sinking keys
                        || key_event.modifiers.contains(KeyModifiers::CONTROL) // Ctrl+C is a hard exit
                            && (key_event.code == KeyCode::Char('c')
                                || key_event.code == KeyCode::Char('C'))
                    {
                        disable_raw_mode().expect("Failed to disable terminal raw mode");
                        execute!(terminal.backend_mut(), LeaveAlternateScreen)
                            .expect("Failed to leave alternate terminal screen");
                        terminal
                            .show_cursor()
                            .expect("Failed to show terminal cursor");
                        return;
                    }
                }
                _ => (),
            }
        }

        let now = Instant::now();
        test_world.update(now.duration_since(then).as_secs_f32());
        then = now;

        terminal
            .draw(|frame| {
                let [raster_area, logger_area, ..] = Layout::default()
                    .constraints([
                        Constraint::Percentage(50), 
                        Constraint::Percentage(50)
                    ])
                    .direction(Direction::Horizontal)
                    .split(frame.size())[..] else { unreachable!() }; 

                frame.render_widget(
                    RasterWidget {
                        scene: &test_world,
                        camera: &camera,
                        raster: &mut raster,
                        buffer: &mut buffer,
                    },
                    raster_area,
                );

                frame.render_widget(
                    TuiLoggerWidget::default()
                        .block(
                            Block::default()
                                .title(" Log ")
                                .border_style(Style::default().fg(TerminalColor::White))
                                .borders(Borders::ALL),
                        )
                        .output_separator('|')
                        .output_timestamp(Some("%H:%M:%S%.3f".to_string()))
                        .output_level(Some(TuiLoggerLevelOutput::Abbreviated))
                        .output_target(false)
                        .output_file(false)
                        .output_line(false)
                        .style_error(Style::default().fg(TerminalColor::Red))
                        .style_debug(Style::default().fg(TerminalColor::Cyan))
                        .style_warn(Style::default().fg(TerminalColor::Yellow))
                        .style_trace(Style::default().fg(TerminalColor::White))
                        .style_info(Style::default().fg(TerminalColor::Green)), 
                    logger_area
                );
            })
            .expect("Failed to draw to terminal");
    }
}

pub struct RasterWidget<'a, R: Rasterable, S: Scene, V: Viewport> {
    pub scene: &'a S,
    pub camera: &'a V,
    pub raster: &'a mut R,
    pub buffer: &'a mut Buffer2D<Color>,
}

impl<'a, R: Rasterable, S: Scene, V: Viewport> tui::widgets::Widget for RasterWidget<'a, R, S, V> {
    fn render(self, area: tui::layout::Rect, buf: &mut tui::buffer::Buffer) {
        if area.area() == 0 {
            return;
        }

        self.raster.rasterize(
            self.scene,
            self.camera,
            self.buffer,
            area.width,
            area.height * 2,
        );

        for y in 0..area.height {
            for x in 0..area.width {
                let top_color = tui_color(*self.buffer.get(x as usize, y as usize * 2));
                let bot_color = tui_color(*self.buffer.get(x as usize, y as usize * 2 + 1));

                let cell = buf.get_mut(x, y);
                cell.set_symbol("▄");
                cell.set_bg(top_color);
                cell.set_fg(bot_color);
            }
        }
    }
}

fn tui_color(color: Color) -> TerminalColor {
    TerminalColor::Rgb(color.r, color.g, color.b)
}