pub trait Solution<G: Genotype, P: Phenotype> {
    fn set_gtype(&mut self, gtype_new: &G) -> ();
    fn set_ptype(&mut self, ptype_new: &P) -> ();
    fn get_gtype(&self) -> &G;
    fn get_ptype(&self) -> &P;
    fn induce_gtype(&mut self) -> ();
    fn induce_ptype(&mut self) -> ();
}

pub trait Population<G: Genotype, P: Phenotype> {
    fn new(n: &usize) -> Self;
    fn get(&self, n: &usize) -> &Solution<G,P>;
    fn set(&mut self, n: &usize, soln: Solution<G,P>) -> ();
    fn len(&self) -> usize;
    fn get_age(&self) -> &f32;
    fn add(&mut self, soln: Solution<G,P>) -> ();
    fn rem(&mut self, n: &usize) -> ();
    fn inc_age(&mut self) -> ();
    fn clear(&mut self) -> ();
}

pub trait Genotype<G: Clone + Copy> {
    fn get(&self) -> &G;
    fn set(&mut self, other: &G) -> ();
}

pub trait Phenotype<P: Clone + Copy> {
    fn get(&self) -> &P;
    fn set(&mut self, other: &P) -> ();
}