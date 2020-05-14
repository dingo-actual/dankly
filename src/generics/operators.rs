use crate::generics::population::Solution;

pub trait Crossover<G: Clone + Copy, P: Clone + Copy> {
    fn crossover(&self, parents: Vec<&Solution<G,P>>) -> Vec<Solution<G,P>>;
}

pub trait LocalSearch<G: Clone + Copy, P: Clone + Copy> {
    fn search(&self, soln: &mut Solution<G,P>) -> ();
}

pub trait Mutation<G: Clone + Copy, P: Clone + Copy> {
    fn mutate(&self, soln: &mut Solution<G,P>) -> ();
}

pub trait SelectParents {
    fn choose(&self, fitnesses: Vec<f32>) -> Vec<Vec<usize>>;
}

pub trait SelectSurvivors {
    fn choose(&self, fitnesses: Vec<f32>) -> Vec<usize>;
}