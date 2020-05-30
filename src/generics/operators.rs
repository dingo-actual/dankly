use crate::generics::population::{Solution, Population, Genotype};
use rand::thread_rng;
use rand_distr::{Uniform, Distribution};

pub struct Crossover<G: Clone, P: Clone> {
    chooser: Box<dyn Fn(&Vec<&f32>) -> Vec<Vec<usize>>>,
    crossover: Box<dyn Fn(&Vec<&Solution<G, P>>) -> Vec<Solution<G,P>>>,
}

impl <G: Clone, P: Clone> Crossover<G,P> {
    pub fn crossover_pop(&self, pop: &Population<G,P>) -> Vec<Solution<G,P>> {
        let fitnesses = pop.get_fitnesses();
        let parent_ixs = (&self.chooser)(&fitnesses);
        let mut out = Vec::new();
        for ixs in parent_ixs.iter() {
            let mut parents = Vec::new();
            for ix in ixs.iter() {
                parents.push(pop.get(ix));
            }
            let mut children = (&self.crossover)(&parents);
            while let Some(child) = children.pop() {
                out.push(child);
            }
        }
        out
    }
}

pub struct LocalSearch<G: Clone, P: Clone> {
    chooser: Box<dyn Fn(&Vec<&f32>) -> (Vec<usize>, Vec<Vec<usize>>)>,
    search: Box<dyn Fn(&mut Solution<G, P>, &Vec<&Solution<G,P>>) -> ()>,
}

impl <G: Clone, P: Clone> LocalSearch<G,P> {
    pub fn search_pop(&self, pop: &mut Population<G,P>) -> () {
        let fitnesses = (&*pop).get_fitnesses();
        let (pt_ixs, nhds) = (&self.chooser)(&fitnesses);
        for (pt_ix, nhd_ixs) in pt_ixs.iter().zip(nhds.iter()) {
            let mut nhd = Vec::new();
            for nhd_ix in nhd_ixs.iter() {
                let nhd_member = pop.get(nhd_ix).clone();
                nhd.push(nhd_member);
            }
            let pt = pop.get_mut(pt_ix);
            {
                (&self.search)(pt, &nhd.iter().collect());
            }
        }
    }
}

fn self_adapt<G: Clone>(gtype: &mut Genotype<G>) -> () {
    let adapt = gtype.get_sa_mut();
    adapt.self_adapt();
}

pub struct Mutation<G: Clone> {
    pub mutation_op: Box<dyn Fn(&mut Genotype<G>) -> ()>,
}

impl<G: Clone> Mutation<G> {
    pub fn mutate<P: Clone>(&self, soln: &mut Solution<G, P>) -> () {
        let gtype = soln.get_gtype_mut();
        let mut rng = thread_rng();
        let unif = Uniform::new(0.0, 1.0);
        let mut_samp = unif.sample(&mut rng);
        if mut_samp < gtype.get_sa().mutation_prob {
            (self.mutation_op)(gtype);
        }
        self_adapt(gtype);
    }
}

pub struct Selector {
    select_fn: Box<dyn Fn(&Vec<&f32>) -> Vec<usize>>,
}

impl Selector {
    pub fn select<G: Clone, P: Clone>(&self, pop: &Population<G,P>) -> Vec<usize> {
        let fitnesses = pop.get_fitnesses();
        (&self.select_fn)(&fitnesses)
    }
}