use crate::generics::operators::*;
use crate::generics::population::*;
use crate::generics::terminator::Terminator;

fn optimize<G: Clone + Copy, P: Clone + Copy>(pop: &Population<G,P>, fitness: &Fitness<G,P>, diversity: &Diversity<G,P>, xover: &Crossover<G,P>, mutate: &Mutation<G,P>, local_search: &LocalSearch<G,P>, select_parents: SelectParents<G,P>, select_survivors: SelectSurvivors<G,P>, term: Terminator<G,P>) -> () {
    let mut terminate = False;
    let fitnesses = Vec::new();
    for n in 0..pop.len() {
        fitnesses.push(fitness.eval(pop.get(n)));
    }
    while !terminate {
        let mut new_solutions = Vec::new();
        // select parents
        let parent_ixs = select_parents.choose(&pop);
        // create new solutions through crossover
        for ixs in parent_ixs.iter() {
            let parents = Vec::new();
            for ix in ixs.iter() {
                parents.push(pop.get(ix));
            }
            let child_vec = xover.crossover(parents);
            while let Some(child) = child_vec.pop() {
                new_solutions.push(child);
            }
        }
        // mutate new solutions
        for solution in new_solutions.iter_mut() {
            mutate.mutate(solution);
        }
        // perform local search on mutated solutions
        for solution in new_solutions.iter_mut() {
            local_search.search(solution);
        }
        // add new solutions to population
        while let Some(solution) = new_solutions.pop() {
            fitnesses.push(fitness.eval(&solution));
            pop.add(solution);
        }
        // choose survivors
        let survivor_ixs = select_survivors.choose(fitnesses);
        let survivors = Vec::new();
        for ix in survivor_ixs.iter() {
            survivors.push(pop.get(ix).clone());
        }
        while let Some(survivor) = survivors.pop() {
            pop.add(survivor);
        }
        pop.inc_age();
        // update term conditions
        div = diversity.eval(&pop);
        term.update(fitnesses, div, pop.get_age());
        terminate = term.check();
    }
}