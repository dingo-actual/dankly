use crate::generics::operators::crossover::Crossover;
use crate::generics::operators::local_search::LocalSearch;
use crate::generics::operators::mutation::Mutation;
use crate::generics::operators::select_parents::SelectParents;
use crate::generics::operators::select_survivors::SelectSurvivors;
use crate::generics::population::population::Population;
use crate::generics::population::solution::Solution;
use crate::generics::terminator::Terminator;
use crate::util::uid::UidStack;

fn optimize<G,P>(pop: &Population<G,P>, xover: &Crossover<G,P>, mutate: &Mutation<G,P>, select_parents: SelectParents<G,P>, select_survivors: SelectSurvivors<G,P>, term: Terminator<G,P>) -> () {
    let mut terminate = False;
    let mut fitnesses = vec![0; pop.len()];
    while !terminate {
        // update fitnesses

        // select parents

        // create new solutions through crossover

        // mutate new solutions

        // perform local search on mutated solutions

        // update term conditions
    }
}