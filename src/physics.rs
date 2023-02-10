mod objects;
use objects::*;

struct RigidBody {
    cube:Cube,
    mass:f64,
    moment_of_intertia:[f64;3],
    velocity:[f64;3],
    angular_velocity:[f64;3],
}

// Detect the collisions between cubes and return the pairs of cubes colliding
fn detect_collisions(world:&Vec<Cube>) -> Vec<(Cube,Cube)>{

}