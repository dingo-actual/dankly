use crate::generics::operators::*;
use crate::generics::population::*;
use crate::generics::terminator::Terminator;

fn optimize<G: Clone + Copy, P: Clone + Copy>(pop: &Population<G,P>, fitness: &Fitness<G,P>, diversity: &Diversity<G,P>, xover: &Crossover<G,P>, mutate: &Mutation<G,P>, local_search: &LocalSearch<G,P>, select_survivors: Select, term: Terminator<G,P>) -> () {
    let mut terminate = False;
    let mut fitnesses = Vec::new();
    for n in 0..pop.len() {
        fitnesses.push(fitness.eval(pop.get(n)));
    }
    while !terminate {
        let mut new_solutions = Vec::new();
        // select parents
        let parent_ixs = xover.choose(&fitnesses);
        // create new solutions through crossover
        let children = xover.crossover_all(&pop, &parent_ixs);
        for child_vec in children.iter() {
            for child in child_vec.iter() {
                new_solutions.push(child);
            }
        }
        // mutate new solutions
        for solution in new_solutions.iter_mut() {
            mutate.mutate(solution);
        }
        // add new solutions to population
        while let Some(solution) = new_solutions.pop() {
            fitnesses.push(fitness.eval(&solution));
            pop.add(solution);
        }
        // perform local search on mutated solutions
        let search_ixs = local_search.choose(&fitnesses);
        local_search.search_all(&mut pop, &search_ixs);
        for ix_vec in search_ixs.iter() {
            fitnesses[ix_vec[0]] = fitness.eval(pop.get(&ix_vec[0]));
        }
        // choose survivors
        let survivor_ixs = select_survivors.choose(fitnesses);
        let survivors = Vec::new();
        for ix_vec in survivor_ixs.iter() {
            survivors.push(*pop.get(ix_vec[0]).copy());
        }
        pop.clear();
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