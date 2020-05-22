use crate::generics::population::{Solution, Population};

pub trait Crossover<G: Clone + Copy, P: Clone + Copy> {
    fn choose(&self, fitnesses: &Vec<f32>) -> ();
    fn crossover(&self, parents: Vec<&Solution<G,P>>) -> Vec<Solution<G,P>>;
    fn crossover_all(&self, pop: &Population<G,P>, parent_ixs: &Vec<Vec<usize>>) -> Vec<Vec<Solution<G,P>>>;
}

pub trait LocalSearch<G: Clone + Copy, P: Clone + Copy> {
    fn choose(&self, fitnesses: &Vec<f32>) -> Vec<Vec<usize>>;
    fn search(&self, soln: &mut Solution<G,P>) -> ();
    fn search_all(&self, pop: &mut Population<G,P>, search_ixs: &Vec<Vec<usize>>) -> ();
}

pub trait Mutation<G: Clone + Copy, P: Clone + Copy> {
    fn mutate_solution(&self, soln: &mut Solution<G,P>) -> ();
}

pub trait Select {
    fn choose(&self, fitnesses: Vec<f32>) -> Vec<Vec<usize>>;
}