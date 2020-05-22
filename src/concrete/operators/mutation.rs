use crate::generics::operators::Mutation;
use crate::concrete::population::solution::*;
use rand::thread_rng;
use rand::distributions::{Bernoulli, Distribution};
use rand_distr::{Normal, Uniform, Distribution};

struct TransposeMutator {
    block_size: usize,
}

impl TransposeMutator {
    pub fn new(block_size: &usize) -> TransposeMutator {
        if block_size == 0 {
            panic!("Cannot have block size zero");
        }
        TransposeMutator {
            block_size: block_size,
        }
    }
    fn get_blocks(&self, max_ix: usize) -> ((usize, usize), (usize, usize)) {
        let mut rng = thread_rng();
        let unif = Uniform::new(0, 1);
        let buf = (unif.sample(&mut rng) * (max_len - 2 * self.block_size + 1)).floor() as usize;
        let lo = (unif.sample(&mut rng) * (max_len - 2 * self.block_size - buf + 1)).floor() as usize;
        let hi = lo + self.block_size + buf;
        ((lo, lo + self.block_size), (hi, hi + block_size))
    }
}

impl Mutation<Permutation, Permutation> for TransposeMutator {
    fn mutate_solution(&self, soln: &mut Solution<Permutation, Permutation>) -> () {
        {
            if self.block_size > soln.len / 2 {
                panic!("Cannot mutate with given block size");
            }
        }
        let len = soln.get_gtype().get_len().copy();
        let ((lo1, hi1), (lo2, hi2)) = self.get_blocks(&len);
        {
            let perm_g = soln.get_gtype_mut();
            perm_g.swap_range((&lo1, &hi1), (&lo2, &hi2));
        }
        {
            let perm_p = soln.get_ptype_mut();
            perm_p.swap_range((&lo1, &hi1), (&lo2, &hi2));
        }
    }
}