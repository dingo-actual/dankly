use crate::population::population::Population;

pub trait Fitness<G: Clone + Copy, P: Clone + Copy> {
    fn eval(&self, pop: Population<G,P>) -> f32;
}