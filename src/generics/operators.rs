use crate::generics::population::{Population, Solution};

pub trait Crossover<G: Clone + Copy, P: Clone + Copy> {
    fn crossover(&self, parents: Vec<&Solution<G,P>>) -> Vec<Solution<G,P>>;
}

pub trait LocalSearch<G: Clone + Copy, P: Clone + Copy> {
    fn search(&self, soln: &mut Solution<G,P>) -> ();
}

pub trait Mutation<G: Clone + Copy, P: Clone + Copy> {
    fn mutate(&self, soln: &mut Solution<G,P>) -> ();
}

pub trait SelectParents<G: Clone + Copy, P: Clone + Copy> {
    fn choose(&self, pop: &Population<G,P>) -> Vec<Vec<usize>>;
}

pub trait SelectSurvivors<G: Clone + Copy, P: Clone + Copy> {
    fn choose(&self, pop: &mut Population<G,P>) -> ();
}