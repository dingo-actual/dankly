use crate::generics::population::Population;

fn call<A, B>(f: &dyn Fn(&A) -> B, args: &A) -> B {
    f(args)
}

pub struct Terminator<G: Clone, P: Clone> {
    max_age: Option<f32>,
    min_diversity: Option<f32>,
    max_iter: usize,
    fitness_atol: Option<f32>,
    fitness_rtol: Option<f32>,
    prev_fitness: f32,
    crnt_fitness: f32,
    crnt_iter: usize,
    crnt_age: f32,
    crnt_diversity: f32,
    div_measure: Box<dyn Fn(&Population<G,P>) -> f32>,
    fitness_agg: Box<dyn Fn(&Population<G,P>) -> f32>,
}

impl <G: Clone, P: Clone> Terminator<G,P> {
    fn check_diversity(&self) -> bool {
        match self.min_diversity {
            None => false,
            Some(m_div) => m_div >= self.crnt_diversity,
        }
    }
    fn check_age(&self) -> bool {
        match self.max_age {
            None => false,
            Some(m_age) => m_age <= self.crnt_age,
        }
    }
    fn check_iter(&self) -> bool {
        self.max_iter <= self.crnt_iter
    }
    fn check_fitness_abs(&self) -> bool {
        match self.fitness_atol {
            None => false,
            Some(atol) => {
                let prog_abs = (self.crnt_fitness - self.prev_fitness).abs();
                prog_abs <= atol
            },
        }
    }
    fn check_fitness_rel(&self) -> bool {
        match self.fitness_rtol {
            None => false,
            Some(rtol) => {
                let prog_rel = (self.crnt_fitness - self.prev_fitness).abs() / self.prev_fitness;
                prog_rel <= rtol
            },
        }
    }
    pub fn check(&self) -> bool {
        self.check_age() || self.check_diversity() || self.check_fitness_abs() || self.check_fitness_rel() || self.check_iter()
    }
    pub fn update(&mut self, pop: &Population<G, P>) -> () {
        self.crnt_age = *pop.get_age();
        self.crnt_iter = self.crnt_iter + 1;
        self.prev_fitness = self.crnt_fitness.clone();
        self.crnt_fitness = call(&self.fitness_agg, pop);
        self.crnt_diversity = call(&self.div_measure, pop);
    }
}