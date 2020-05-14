use crate::population::population::Population;

pub trait Diversity<G: Clone + Copy, P: Clone + Copy> {
    fn eval(&self, pop: &Population<G,P>) -> f32;
}