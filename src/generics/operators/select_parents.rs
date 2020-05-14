use crate::population::population::Population;
use crate::metrics::fitness::Fitness;

pub trait SelectParents<G: Clone + Copy, P: Clone + Copy> {
    fn choose(&self, pop: &Population<G,P>, fit: &Fitness<G,P>) -> Vec<usize>;
}