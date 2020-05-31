use crate::generics::operators::Mutation;
use crate::generics::population::Genotype;
use crate::concrete::population::solution::*;
use rand::thread_rng;
use rand::seq::SliceRandom;
use rand_distr::{Normal, Uniform, Distribution};
use ndarray::Array1;

fn randint(lo: usize, hi: usize) -> usize {
    let mut rng = thread_rng();
    let unif = Uniform::new(lo, hi);
    unif.sample(&mut rng)
}

fn get_blocks(block_size: &usize, max_ix: &usize) -> ((usize, usize), (usize, usize)) {
    let buf = randint(0, max_ix - 2 * block_size + 1);
    let lo = randint(0, max_ix - 2 * block_size - buf + 1);
    let hi = lo + block_size + buf;
    ((lo, lo + block_size), (hi, hi + block_size))
}

// mutation for Permutation

fn transpose_mutation_factory(block_size: usize) -> Box<dyn Fn(&mut Permutation) -> ()>{
    fn transpose_mutation(perm: &mut Permutation, block_size: usize) -> () {
        let max_ix = perm.get_len().clone();
        let ((lo1, hi1), (lo2, hi2)) = get_blocks(&block_size, &max_ix);
        perm.swap_range((&lo1, &hi1), (&lo2, &hi2));
    }
    Box::new(move |perm: &mut Permutation| transpose_mutation(perm, block_size))
}

fn transpose_mutation_op(gtype: &mut Genotype<Permutation>) -> () {
    let block_size;
    {
        let sa = gtype.get_sa();
        block_size = sa.mutation_strength.clone() as usize;
    }
    let perm = gtype.get_genes_mut();
    let op = transpose_mutation_factory(block_size);
    op(perm);
}

pub fn get_transpose_mutator() -> Mutation<Permutation> {
    Mutation{mutation_op: Box::new(transpose_mutation_op)}
}

// mutation for RealVec

fn normal_vect(mu: f32, sigma: f32, len: usize) -> Array1<f32> {
    let mut rng = thread_rng();
    let normal = Normal::new(mu, sigma).unwrap();
    let mut out_v: Vec<f32> = Vec::with_capacity(len);
    for n in 0..len {
        out_v[n] = normal.sample(&mut rng);
    }
    Array1::from(out_v)
}

fn normal_mutation_factory(sigma: f32) -> Box<dyn Fn(&mut Array1<f32>) -> ()> {
    fn normal_mutation(vect: &mut Array1<f32>, sigma: f32) -> () {
        let offset = normal_vect(0.0, sigma, vect.len().clone());
        for n in  0..vect.len() {
            vect[n] = vect[n] + offset[n];
        }
    }
    Box::new(move |vect: &mut Array1<f32>| normal_mutation(vect, sigma))
}

fn normal_mutation_op(gtype: &mut Genotype<Array1<f32>>) -> () {
    let sigma: f32;
    {
        let sa = gtype.get_sa();
        sigma = sa.mutation_strength;
    }
    let vect = gtype.get_genes_mut();
    let op = normal_mutation_factory(sigma);
    op(vect);
}

pub fn get_normal_mutator() -> Mutation<Array1<f32>> {
    Mutation{mutation_op: Box::new(normal_mutation_op)}
}

// mutation for BoundedRealVec

fn bounded_normal_mutation_factory() -> Box<dyn Fn(&mut BoundedRealVec) -> ()> {
    fn normal_mutation(vect: &mut BoundedRealVec) -> () {
        let offset = BoundedRealVec::new_normal(&vect.len(), vect.lbound(), vect.ubound());
        *vect = vect.clone() + offset;
    }
    Box::new(move |vect: &mut BoundedRealVec| normal_mutation(vect))
}

fn bounded_normal_mutation_op(gtype: &mut Genotype<BoundedRealVec>) -> () {
    let vect = gtype.get_genes_mut();
    let op = bounded_normal_mutation_factory();
    op(vect);
}

pub fn get_bounded_normal_mutator() -> Mutation<BoundedRealVec> {
    Mutation{mutation_op: Box::new(bounded_normal_mutation_op)}
}

// mutation for Categorical sequence

fn category_mutation_factory(n: usize) -> Box<dyn Fn(&mut CategSeq) -> ()> {
    fn change_categories(catseq: &mut CategSeq, n: usize) -> () {
        if n >= catseq.len() {
            panic!("Cannot change more than all items in sequence")
        }
        let mut ixs: Vec<usize> = Vec::with_capacity(n);
        for ix in 0..catseq.len() {
            ixs.push(ix);
        }
        let mut rng = thread_rng();
        ixs.shuffle(&mut rng);
        for k in 0..n {
            let ix = ixs[k];
            let unif = Uniform::new(0, catseq.n_cat(ix));
            let new_cat = unif.sample(&mut rng);
            let cat = catseq.get_mut(ix);
            *cat = new_cat;
        }
    }
    Box::new(move |catseq: &mut CategSeq| change_categories(catseq, n))
}

fn category_mutation_op(gtype: &mut Genotype<CategSeq>) -> () {
    let k: usize;
    {
        let sa = gtype.get_sa();
        k = sa.mutation_strength as usize;
    }
    let catseq = gtype.get_genes_mut();
    let op = category_mutation_factory(k);
    op(catseq);
}

pub fn get_category_mutator() -> Mutation<CategSeq> {
    Mutation{mutation_op: Box::new(category_mutation_op)}
}