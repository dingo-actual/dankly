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
            return true;
        }
        match self.progress_atol {
            Some(atol) => {
                if progress < atol {
                    return true;
                }
            },
            None => {},
        }
        match self.progress_rtol {
            Some(rtol) => {
                if progress / self.last_fitness < rtol {
                    return true;
                }
            },
            None => {},
        }
        match self.max_age {
            Some(age) => {
                if self.crnt_age >= age {
                    return true;
                } 
            },
            None => {},
        }
        match self.min_diversity {
            Some(diverse) => {
                if self.crnt_diversity <= diverse {
                    return true;
                }
            },
            None => {},
        }
        false
    }
    fn update(&mut self, fitnesses: Vec<f32>, diversity: f32, age: f32) -> () {
        self.crnt_iter = self.crnt_iter + 1;
        self.crnt_diversity = diversity;
        self.crnt_age = pop.age;
        self.last_fitness = self.crnt_fitness.copy();
        
        //compute mean of fitnesses
        self.crnt_fitness = 0;
        let mut denom = 0;
        for fit in fitnesses.iter() {
            self.crnt_fitness = self.crnt_fitness + fit;
            denom = denom + 1;
        }
        self.crnt_fitness = self.crnt_fitness / denom;
    }
}