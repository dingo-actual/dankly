pub trait Solution<G: Clone + Copy, P: Clone + Copy> {
    fn new() -> Self;
    fn set_gtype(&mut self, gtype_new: &G) -> ();
    fn set_ptype(&mut self, ptype_new: &P) -> ();
    fn get_gtype(&self) -> &G;
    fn get_ptype(&self) -> &P;
    fn induce_gtype(&mut self) -> ();
    fn induce_ptype(&mut self) -> ();
    fn update_fitness(&mut self) -> ();
    fn get_fitness(&self) -> f32;
}

pub trait Population<G: Clone + Copy, P: Clone + Copy> {
    fn new(n: &usize) -> Self;
    fn get(&self, n: &usize) -> &Solution<G,P>;
    fn set(&mut self, n: &usize, soln: Solution<G,P>) -> ();
    fn len(&self) -> usize;
    fn get_fitnesses(&self) -> Vec<f32>;
    fn get_fitness(&self) -> f32;
    fn get_diversity(&self) -> f32;
    fn get_age(&self) -> f32;
    fn update_fitness(&mut self, n: &usize) -> ();
    fn update_diversity(&mut self);
    fn add(&mut self, soln: Solution<G,P>) -> ();
    fn rem(&mut self, n: &usize) -> ();
    fn inc_age(&mut self) -> ();
}