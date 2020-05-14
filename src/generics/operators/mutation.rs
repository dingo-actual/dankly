use crate::generics::population::solution::Solution;

pub trait Mutation<G: Clone + Copy, P: Clone + Copy> {
    fn mutate(&self, soln: &mut Solution<G,P>) -> ();
}