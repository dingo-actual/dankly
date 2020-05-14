use crate::generics::metrics::population::Population;

pub trait Terminator<G,P> {
    fn check(&self) -> bool;
    fn update(&mut self, pop: &Population<G,P>) -> ();
}