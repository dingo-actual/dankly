use crate::generics::population::*;

pub trait Diversity<G: Clone + Copy, P: Clone + Copy> {
    fn eval(&self, pop: &Population<G,P>) -> f32;
}

pub trait Fitness<G: Clone + Copy, P: Clone + Copy> {
    fn eval(&self, soln: &Solution<G,P>) -> f32;
}