use crate::generics::population::solution::Solution;

pub trait Crossover<G: Clone + Copy, P: Clone + Copy> {
    fn crossover(&self, parents: Vec<&Solution<G,P>>) -> Vec<Solution<G,P>>;
}