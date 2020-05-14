use crate::generics::terminator::Terminator;
use crate::generics::population::Population;

pub struct GeneralTerminator {
    progress_rtol: Option<f32>,
    progress_atol: Option<f32>,
    max_iter: usize,
    max_age: Option<f32>,
    min_diversity: Option<f32>,
    last_fitness: f32,
    crnt_fitness: f32,
    crnt_iter: usize,
    crnt_age: f32,
    crnt_diversity: f32,
}

impl GeneralTerminator {
    pub fn new(max_iter: usize, rtol: Option<f32>, atol: Option<f32>, max_age: Option<f32>, min_diversity: Option<f32>) -> GeneralTerminator {
        GeneralTerminator {
            progress_atol: atol,
            progress_rtol: rtol,
            max_iter: max_iter,
            max_age: max_age,
            min_diversity: min_diversity,
            last_fitness: f32::MIN,
            crnt_fitness: f32::MIN,
            crnt_iter: 0,
            crnt_age: 0,
        }
    }
}

impl Terminator for GeneralTerminator {
    fn check(&self) -> bool {
        let progress = self.crnt_fitness - self.last_fitness;
        if self.crnt_iter == self.max_iter {
            return True;
        }
        match self.progress_atol {
            Some(atol) => {
                if progress < atol {
                    return True;
                }
            },
            None => {},
        }
        match self.progress_rtol {
            Some(rtol) => {
                if progress / self.last_fitness < rtol {
                    return True;
                }
            },
            None => {},
        }
        match self.max_age {
            Some(age) => {
                if self.crnt_age >= age {
                    return True;
                } 
            },
            None => {},
        }
        match self.min_diversity {
            Some(diverse) => {
                if self.crnt_diversity <= diverse {
                    return True;
                }
            },
            None => {},
        }
        return False;
    }
    fn update(&mut self, pop: &Population<G,P>) -> () {
        self.crnt_iter = self.crnt_iter + 1;
        self.last_fitness = self.crnt_fitness.copy();
        self.crnt_fitness = pop.get_fitness();
        self.crnt_diversity = pop.get_diversity();
        self.crnt_age = pop.get_age();
    }
}