use crate::generics::population::solution::Solution;

pub trait LocalSearch<G: Clone + Copy, P: Clone + Copy> {
    fn search(&self, soln: &mut Solution<G,P>) -> ();
}