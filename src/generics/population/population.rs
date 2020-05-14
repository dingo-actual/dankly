use solution::{Solution};

pub trait Population<G: Clone + Copy, P: Clone + Copy> {
    fn new(n: &usize) -> Self;
    fn get(&self, n: &usize) -> &Solution<G,P>;
    fn set(&mut self, n: &usize, soln: Solution<G,P>) -> ();
    fn len(&self) -> usize;
    fn get_fitnesses(&self) -> Vec<f32>;
    fn get_diversity(&self) -> f32;
    fn update_fitness(&mut self, n: &usize) -> ();
    fn update_diversity(&mut self);
    fn add(&mut self, soln: Solution<G,P>) -> ();
    fn rem(&mut self, n: &usize) -> ();
    fn inc_age(&mut self) -> ();
}