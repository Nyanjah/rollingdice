mod body;
mod camera;
mod raster;
mod transform;
mod world;
mod output;

fn main() {
   if let Some("term") = std::env::args().nth(1).as_ref().map(|x| x.as_str()) {
      output::term::run();
   } else {
      output::prog::run();
   }
}