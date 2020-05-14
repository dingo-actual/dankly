use crate::population::population::Population;
use crate::metrics::fitness::Fitness;

pub trait SelectSurvivors<G: Clone + Copy, P: Clone + Copy> {
    fn choose(&self, pop: &mut Population<G,P>, fit: &Fitness<G,P>) -> Population<G,P>;
}