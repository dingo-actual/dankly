use crate::generics::population::*;

struct GeneralPopulation<G: Clone + Copy, P: Clone + Copy> {
    solutions: Vec<Solution<G,P>>,
    age: f32,
}

impl Population<G, P> for GeneralPopulation {
    fn new() -> GeneralPopulation<G,P> {
        GeneralPopulation {
            solutions: Vec::new(),
            age: 0,
        }
    }
    fn get(&self, n: &usize) -> &Solution<G,P> {
        &self.solutions[n]
    }
    fn set(&mut self, n: &usize, soln: Solution<G,P>) -> () {
        self.solutions[n] = soln;
    }
    fn len(&self) -> usize {
        self.solutions.len();
    }
    fn get_age(&self) -> &f32 {
        &self.age
    }
    fn add(&mut self, soln: Solution<G,P>) -> () {
        self.solutions.push(soln);
    }
    fn rem(&mut self, n: &usize) {
        self.solutions.remove(n);
    }
    fn inc_age(&mut self) -> () {
        self.age = self.age + 1;
    }
    fn clear(&mut self) -> () {
        self.solutions = Vec::new();
    }
}