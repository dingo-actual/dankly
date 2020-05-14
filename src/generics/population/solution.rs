pub trait Solution<G: Clone + Copy, P: Clone + Copy> {
    fn new() -> Self;
    fn set_gtype(&mut self, gtype_new: &G) -> ();
    fn set_ptype(&mut self, ptype_new: &P) -> ();
    fn get_gtype(&self) -> &G;
    fn get_ptype(&self) -> &P;
    fn induce_gtype(&mut self) -> ();
    fn induce_ptype(&mut self) -> ();
    fn update_fitness(&mut self) -> ();
    fn get_fitness(&self) -> f32;
}