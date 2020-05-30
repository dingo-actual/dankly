use rand::thread_rng;
use rand_distr::{Normal, LogNormal, Binomial, Distribution};

#[derive(Clone)]
pub enum AdaptDist {
    Const,
    Normal,
    LogNormal,
    Binomial,
}


#[derive(Clone)]
pub struct SelfAdaptation {
    pub prob_distr: AdaptDist,
    pub strength_distr: AdaptDist,
    pub mutation_prob: f32,
    pub mutation_strength: f32,
    max_mutation_strength: Option<f32>,
    min_mutation_strength: Option<f32>,
}

impl SelfAdaptation {
    pub fn new(prob_distr: AdaptDist, strength_distr: AdaptDist, mutation_prob: f32, mutation_strength: f32, max_mutation_strength: Option<f32>, min_mutation_strength: Option<f32>) -> SelfAdaptation {
        SelfAdaptation {
            prob_distr: prob_distr,
            strength_distr: strength_distr,
            mutation_prob: mutation_prob,
            mutation_strength: mutation_strength,
            max_mutation_strength: max_mutation_strength,
            min_mutation_strength: min_mutation_strength,
        }
    }
    pub fn self_adapt(&mut self) -> () {
        let mut rng = thread_rng();
        match &self.prob_distr {
            AdaptDist::Const => (),
            AdaptDist::Normal => {
                let normal = Normal::new(0.0, 1.0).unwrap();
                self.mutation_prob = self.mutation_prob + normal.sample(&mut rng);
                if self.mutation_prob > 1.0 {
                    self.mutation_prob = 1.0 - 1e-10;
                } else if self.mutation_prob < 0.0 {
                    self.mutation_prob = 1e-10;
                }
            },
            AdaptDist::LogNormal => {
                let lognormal = LogNormal::new(0.0, 1.0).unwrap();
                self.mutation_prob = self.mutation_prob * lognormal.sample(&mut rng);
                if self.mutation_prob > 1.0 {
                    self.mutation_prob = 1.0 - 1e-10;
                }
            },
            AdaptDist::Binomial => {
                panic!("Invalid sampling distribution for mutation probability");
            }
        }
        match &self.strength_distr {
            AdaptDist::Const => (),
            AdaptDist::Normal => {
                let normal = Normal::new(0.0, 1.0).unwrap();
                self.mutation_strength = self.mutation_strength + normal.sample(&mut rng);
                if self.max_mutation_strength.is_some() && self.mutation_prob > self.max_mutation_strength.unwrap() {
                    self.mutation_prob = self.max_mutation_strength.unwrap() - 1e-10;
                } else if self.min_mutation_strength.is_some() && self.mutation_prob < self.min_mutation_strength.unwrap() {
                    self.mutation_prob = self.min_mutation_strength.unwrap() + 1e-10;
                }
            },
            AdaptDist::LogNormal => {
                let lognormal = LogNormal::new(0.0, 1.0).unwrap();
                self.mutation_strength = self.mutation_strength * lognormal.sample(&mut rng);
                if self.max_mutation_strength.is_some() && self.mutation_prob > self.max_mutation_strength.unwrap() {
                    self.mutation_prob = self.max_mutation_strength.unwrap() - 1e-10;
                } else if self.min_mutation_strength.is_some() && self.mutation_prob < self.min_mutation_strength.unwrap() {
                    self.mutation_prob = self.min_mutation_strength.unwrap() + 1e-10;
                }
            },
            AdaptDist::Binomial => {
                if self.max_mutation_strength.is_none() || self.min_mutation_strength.is_none() {
                    panic!("Binomial self-adaptation requires mutation strength bounds");
                }
                let max_str = self.max_mutation_strength.unwrap() as isize;
                let min_str = self.min_mutation_strength.unwrap() as isize;
                let range = (max_str - min_str) as u64;
                let binomial = Binomial::new(range, 0.5).unwrap();
                self.mutation_strength = (binomial.sample(&mut rng) as isize + min_str) as f32;
            }
        }
    }
}

#[derive(Clone)]
pub struct Genotype<G: Clone> {
    genes: G,
    self_adapt: SelfAdaptation,
}

impl <G: Clone> Genotype<G> {
    pub fn new(genes: G, self_adapt: SelfAdaptation) -> Genotype<G> {
        Genotype {
            genes: genes,
            self_adapt: self_adapt,
        }
    }
    pub fn get_genes(&self) -> &G {
        &self.genes
    }
    pub fn get_genes_mut(&mut self) -> &mut G {
        &mut self.genes
    }
    pub fn get_sa(&self) -> &SelfAdaptation {
        &self.self_adapt
    }
    pub fn get_sa_mut(&mut self) -> &mut SelfAdaptation {
        &mut self.self_adapt
    }
}

#[derive(Clone)]
pub struct Solution<G: Clone, P: Clone> {
    genotype: Genotype<G>,
    phenotype: P,
    fitness: f32,
}

impl <G: Clone, P: Clone> Solution<G,P> {
    pub fn get_gtype(&self) -> &Genotype<G> {
        &self.genotype
    }
    pub fn get_gtype_mut(&mut self) -> &mut Genotype<G> {
        &mut self.genotype
    }
    pub fn get_ptype(&self) -> &P {
        &self.phenotype
    }
    pub fn get_ptype_mut(&mut self) -> &mut P {
        &mut self.phenotype
    }
    pub fn get_fitness(&self) -> &f32 {
        &self.fitness
    }
    pub fn get_fitness_mut(&mut self) -> &mut f32 {
        &mut self.fitness
    }
}

pub struct Population<G: Clone, P: Clone> {
    solutions: Vec<Solution<G, P>>,
    age: f32,
}

impl <G: Clone, P: Clone> Population<G,P> {
    pub fn new() -> Population<G,P> {
        Population {
            solutions: Vec::new(),
            age: 0.0,
        }
    }
    pub fn get(&self, n: &usize) -> &Solution<G,P> {
        &self.solutions[*n]
    }
    pub fn get_mut(&mut self, n: &usize) -> &mut Solution<G,P> {
        &mut self.solutions[*n]
    }
    pub fn add(&mut self, soln: Solution<G,P>) -> () {
        self.solutions.push(soln);
    }
    pub fn get_age(&self) -> &f32 {
        &self.age
    }
    pub fn get_age_mut(&mut self) -> &mut f32 {
        &mut self.age
    }
    pub fn len(&self) -> usize {
        self.solutions.len()
    }
    pub fn clear(&mut self) -> () {
        self.solutions = Vec::new();
    }
    pub fn rem(&mut self, n: &usize) -> () {
        self.solutions.remove(*n);
    }
    pub fn get_fitnesses(&self) -> Vec<&f32> {
        let mut out = Vec::new();
        for n in 0..self.len() {
            out.push(self.get(&n).get_fitness());
        }
        out
    }
}