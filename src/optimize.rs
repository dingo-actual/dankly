use crate::generics::operators::*;
use crate::generics::population::*;
use crate::generics::terminator::Terminator;

pub fn optimize<G: Clone, P: Clone>(pop: &mut Population<G,P>, xover: &Crossover<G,P>, mutate: &Mutation<G>, local_search: &LocalSearch<G,P>, select_survivors: &Selector, fitness: Box<dyn Fn(&Solution<G,P>) -> f32>, term: &mut Terminator<G,P>) -> () {
    let mut terminate = false;
    while !terminate {
        let mut new_solutions = Vec::new();
        // create new solutions through crossover
        {
            let mut children = xover.crossover_pop(&*pop);
            while let Some(child) = children.pop() {
                new_solutions.push(child);
            }
        }
        // mutate new solutions
        for solution in new_solutions.iter_mut() {
            (&mutate.mutation_op)(solution.get_gtype_mut());
        }
        // add new solutions to population
        while let Some(mut solution) = new_solutions.pop() {
            let fit = fitness(&solution).clone();
            {
                let soln_fit = solution.get_fitness_mut();
                *soln_fit = fit;
            }
            pop.add(solution);
        }
        // perform local search on population
        {
            local_search.search_pop(pop);
        }
        // choose survivors
        {
            let survivor_ixs = select_survivors.select(&*pop);
            let mut survivors = Vec::new();
            for ix in survivor_ixs.iter() {
                survivors.push(pop.get(ix).clone());
            }
            {
                pop.clear();
            }
            while let Some(survivor) = survivors.pop() {
                pop.add(survivor);
            }
        }
        {
            // TODO: make age the median over solution ages -- needs solution age implemented
            let age = pop.get_age_mut();
            *age = *age + 1.0;
        }
        // update term conditions
        {
            term.update(&*pop);
        }
        terminate = term.check();
    }
}