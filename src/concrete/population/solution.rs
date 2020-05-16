use crate::generics::population::Solution;
use rand::seq::SliceRandom;
use rand::thread_rng;

//permutation

#[derive(Clone, Copy)]
struct Permutation {
    seq: Vec<usize>,
    len: usize,
}

impl Permutation {
    fn new(len: &usize) -> Permutation {
        let mut rng = thread_rng();
        let mut seq = Vec::with_capacity(len);
        for n in 0..len {
            seq.push(n);
        }
        seq.shuffle(&mut rng);
        Permutation{
            seq: seq,
            len: len,
        }
    }
    fn swap(&mut self, ix1: &usize, ix2: &usize) -> Result<&usize, &str> {
        if self.check_bounds(ix1) && self.check_bounds(ix2) {
            if ix1 != ix2 {
                x1 = self.seq[*ix1].copy();
                self.seq[*ix1] = self.seq[*ix2].copy();
                self.seq[*ix2] = x1;
                return Ok(&x1);
            } else {
                return Ok(&0);
            }
        } else if self.check_bounds(ix1) {
            return Err("Index 2 out of bounds");
        } else if self.check_bounds(ix2) {
            return Err("Index 1 out of bounds");
        } else {
            return Err("Both indices out of bounds");
        }
    }
    fn swap_range(&mut self, bounds1: (&usize, &usize), bounds2: (&usize, &usize)) -> Result<usize, &str> {
        let (lo1, hi1) = bounds1;
        let (lo2, hi2) = bounds2;
        if lo1 >= hi1 {
            return Err(format!("Invalid range ({}, {})", lo1, hi1));
        }
        if lo2 >= hi2 {
            return Err(format!("Invalid range ({}, {})", lo2, hi2));
        }
        if lo1 >= lo2 || hi1 >= hi2 {
            return Err("Invalid ranges");
        }
        let ixs = [lo1, lo2, hi1, hi2];
        for ix in ixs.iter() {
            if !self.check_bounds(ix) {
                return Err(format!("Invalid index: {}", ix));
            }
        }
        if lo2 < hi1 {
            return Err("Overlapping ranges");
        }
        let mut out = Vec::with_capacity(self.len);
        for n in 0..lo1 {
            out.push(self.seq[n]);
        }
        for n in lo2..hi2 {
            out.push(self.seq[n]);
        }
        for n in hi1..lo2 {
            out.push(self.seq[n]);
        }
        for n in lo1..hi1 {
            out.push(self.seq[n]);
        }
        for n in hi2..self.len {
            out.push(self.seq[n]);
        }
        self.seq = out;
        Ok(&0);
    }
    fn check_bounds(&self, ix: &usize) -> bool {
        return ix >= 0 && ix < self.len;
    }
}

struct SolnIdentPermu {
    gtype: Permutation,
    ptype: Permutation,
    len: usize,
}

impl SolnIdentPermu {
    fn new(len: &usize) {
        let perm = Permutation::new();
        SolnIdentPermu {
            gtype: perm,
            ptype: perm.copy(),
            len: *len,
        }
    }
}

impl Solution<Permutation, Permutation> for SolnIdentPermu {
    fn set_gtype(&mut self, gtype_new: &Permutation) -> () {
        self.gtype = *gtype_new.copy();
    }
    fn set_ptype(&mut self, ptype_new: &Permutation) -> () {
        self.ptype = *ptype_new.copy();
    }
    fn get_gtype(&self) -> &Permutation {
        &self.gtype
    }
    fn get_ptype(&self) -> &Permutation {
        &self.ptype
    }
    fn induce_gtype(&mut self) -> () {
        self.gtype = self.ptype.copy();
    }
    fn induce_ptype(&mut self) -> () {
        self.ptype = self.gtype.copy();
    }
}

//real vector

//bounded real vector

//binary string

//categorical sequence