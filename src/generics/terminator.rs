pub trait Terminator {
    fn check(&self) -> bool;
    fn update(&mut self, fitnesses: Vec<f32>, diversity: f32, age: f32) -> ();
}