use std::ops::{Add, Sub, Mul, Div};
use rand::seq::SliceRandom;
use rand::thread_rng;
use rand_distr::{Normal, Uniform, Distribution};
use ndarray::Array1;

//permutation

#[derive(Clone)]
pub struct Permutation {
    seq: Array1<usize>,
    len: usize,
}

impl Permutation {
    pub fn new_random(len: &usize) -> Permutation {
        let mut rng = thread_rng();
        let mut seq = Vec::with_capacity(*len);
        for n in 0..*len {
            seq.push(n);
        }
        seq.shuffle(&mut rng);
        Permutation{
            seq: Array1::from(seq),
            len: *len,
        }
    }
    pub fn get_len(&self) -> &usize {
        &self.len
    }
    pub fn swap(&mut self, ix1: &usize, ix2: &usize) -> () {
        if self.check_bounds(ix1) && self.check_bounds(ix2) {
            if ix1 != ix2 {
                let x1 = self.seq[*ix1].clone();
                self.seq[*ix1] = self.seq[*ix2].clone();
                self.seq[*ix2] = x1;
            }
        } else if self.check_bounds(ix1) {
            panic!("Index 2 out of bounds");
        } else if self.check_bounds(ix2) {
            panic!("Index 1 out of bounds");
        } else {
            panic!("Both indices out of bounds");
        }
    }
    pub fn swap_range(&mut self, bounds1: (&usize, &usize), bounds2: (&usize, &usize)) -> () {
        let (lo1, hi1) = bounds1;
        let (lo2, hi2) = bounds2;
        if lo1 >= hi1 {
            panic!("Invalid range");
        }
        if lo2 >= hi2 {
            panic!("Invalid range");
        }
        if lo1 >= lo2 || hi1 >= hi2 {
            panic!("Invalid range");
        }
        let ixs = [lo1, lo2, hi1, hi2];
        for ix in ixs.iter() {
            if !self.check_bounds(ix) {
                panic!("Invalid index");
            }
        }
        if lo2 < hi1 {
            panic!("Overlapping ranges");
        }
        let mut out = Vec::with_capacity(self.len);
        for n in 0..*lo1 {
            out.push(self.seq[n]);
        }
        for n in *lo2..*hi2 {
            out.push(self.seq[n]);
        }
        for n in *hi1..*lo2 {
            out.push(self.seq[n]);
        }
        for n in *lo1..*hi1 {
            out.push(self.seq[n]);
        }
        for n in *hi2..self.len {
            out.push(self.seq[n]);
        }
        self.seq = Array1::from(out);
    }
    fn check_bounds(&self, ix: &usize) -> bool {
        return *ix < self.len;
    }
}

//bounded real vector

#[derive(Clone)]
pub struct BoundedRealVec {
    len: usize,
    seq: Array1<f32>,
    ubounds: Array1<f32>,
    lbounds: Array1<f32>,
}

impl Add<BoundedRealVec> for BoundedRealVec {
    type Output = BoundedRealVec;

    fn add(self, other: BoundedRealVec) -> BoundedRealVec {
        if !self.check_compatible(&other) {
            panic!("Tried to add incompatible bounded vectors");
        }
        let mut out = self.seq + other.seq;
        BoundedRealVec::clip(&mut out, &self.lbounds, &self.ubounds);
        BoundedRealVec{
            len: self.len, 
            seq: out,
            lbounds: self.lbounds.clone(),
            ubounds: self.ubounds.clone(),
        }
    }
}

impl Sub<BoundedRealVec> for BoundedRealVec {
    type Output = BoundedRealVec;

    fn sub(self, other: BoundedRealVec) -> BoundedRealVec {
        if !self.check_compatible(&other) {
            panic!("Tried to subtract incompatible bounded vectors");
        }
        let mut out = self.seq - other.seq;
        BoundedRealVec::clip(&mut out, &self.lbounds, &self.ubounds);
        BoundedRealVec{
            len: self.len, 
            seq: out,
            lbounds: self.lbounds.clone(),
            ubounds: self.ubounds.clone(),
        }
    }
}

impl Mul<BoundedRealVec> for BoundedRealVec {
    type Output = BoundedRealVec;

    fn mul(self, rhs: BoundedRealVec) -> BoundedRealVec {
        if !self.check_compatible(&rhs) {
            panic!("Tried to multiply incompatible bounded vectors");
        }
        let mut out = self.seq * rhs.seq;
        BoundedRealVec::clip(&mut out, &self.lbounds, &self.ubounds);
        BoundedRealVec{
            len: self.len, 
            seq: out,
            lbounds: self.lbounds.clone(),
            ubounds: self.ubounds.clone(),
        }
    }
}

impl Div<BoundedRealVec> for BoundedRealVec {
    type Output = BoundedRealVec;

    fn div(self, rhs: BoundedRealVec) -> BoundedRealVec {
        if !self.check_compatible(&rhs) {
            panic!("Tried to divide incompatible bounded vectors");
        }
        let mut out = self.seq / rhs.seq;
        BoundedRealVec::clip(&mut out, &self.lbounds, &self.ubounds);
        BoundedRealVec{
            len: self.len, 
            seq: out,
            lbounds: self.lbounds.clone(),
            ubounds: self.ubounds.clone(),
        }
    }
}


impl BoundedRealVec {
    fn check_compatible(&self, other: &BoundedRealVec) -> bool {
        if self.len != other.len {
            return false
        }
        for (lo1, lo2) in self.lbounds.iter().zip(other.lbounds.iter()) {
            if lo1 != lo2 {
                return false
            }
        }
        for (hi1, hi2) in self.ubounds.iter().zip(other.ubounds.iter()) {
            if hi1 != hi2 {
                return false
            }
        }
        true
    }
    fn check_bounds(len: &usize, lbounds: &Array1<f32>, ubounds: &Array1<f32>) -> bool {
        if ubounds.len() != *len || lbounds.len() != *len {
            return false
        }
        for (lo, hi) in lbounds.iter().zip(ubounds.iter()) {
            if lo > hi {
                return false
            }
        }
        true
    }
    fn within_bounds(seq: &Array1<f32>, lbounds: &Array1<f32>, ubounds: &Array1<f32>) -> bool {
        if seq.len() != lbounds.len() || seq.len() != ubounds.len() {
            return false;
        }
        for (n, x) in seq.iter().enumerate() {
            if *x < lbounds[n] || *x > ubounds[n] {
                return false;
            }
        }
        true
    }
    pub fn clip(seq: &mut Array1<f32>, lbounds: &Array1<f32>, ubounds: &Array1<f32>) -> () {
        for n in 0..seq.len() {
            if seq[n] < lbounds[n] {
                seq[n] = lbounds[n];
            } else if seq[n] > ubounds[n] {
                seq[n] = ubounds[n];
            }
        }
    }
    fn bounded_normal(lbound: &f32, ubound: &f32) -> f32 {
        let mut rng = thread_rng();
        let mu = (ubound + lbound) / 2.0;
        let sigma = (ubound - lbound) / 4.0;
        let normal = Normal::new(mu, sigma).unwrap();
        let mut samp = normal.sample(&mut rng);
        while samp < *lbound || samp > *ubound {
            samp = normal.sample(&mut rng);
        }
        samp
    }
    fn new_from(len: &usize, lbounds: &Array1<f32>, ubounds: &Array1<f32>, seq: Array1<f32>) -> BoundedRealVec {
        if !BoundedRealVec::check_bounds(len, lbounds, ubounds) {
            panic!("Invalid bounds");
        }
        if !BoundedRealVec::within_bounds(&seq, lbounds, ubounds) {
            panic!("Attempted creation with out of bounds values")
        }
        BoundedRealVec{
            len: *len, 
            seq: seq,
            ubounds: ubounds.clone(),
            lbounds: lbounds.clone(),
        }
    }
    pub fn new_lo(len: &usize, lbounds: &Array1<f32>, ubounds: &Array1<f32>) -> BoundedRealVec {
        BoundedRealVec::new_from(len, lbounds, ubounds, lbounds.clone())
    }
    pub fn new_hi(len: &usize, lbounds: &Array1<f32>, ubounds: &Array1<f32>) -> BoundedRealVec {
        BoundedRealVec::new_from(len, lbounds, ubounds, ubounds.clone())
    }
    pub fn new_normal(len: &usize, lbounds: &Array1<f32>, ubounds: &Array1<f32>) -> BoundedRealVec {
        let mut seq = Vec::with_capacity(*len);
        for n in 0..*len {
            seq.push(BoundedRealVec::bounded_normal(&lbounds[n], &ubounds[n]))
        }
        BoundedRealVec::new_from(len, lbounds, ubounds, Array1::from(seq))
    }
    pub fn new_uniform(len: &usize, lbounds: &Array1<f32>, ubounds: &Array1<f32>) -> BoundedRealVec {
        let mut rng = thread_rng();
        let unif = Uniform::new_inclusive(0.0, 1.0);
        let mut seq = Vec::with_capacity(*len);
        for n in 0..*len {
            let samp = unif.sample(&mut rng);
            let scaled_samp = lbounds[n] + samp * (ubounds[n] - lbounds[n]);
            seq.push(scaled_samp);
        }
        BoundedRealVec::new_from(len, lbounds, ubounds, Array1::from(seq))
    }
}

//categorical sequence

#[derive(Clone)]
pub struct CategSeq {
    seq: Array1<usize>,
    len: usize,
    n_categories: Array1<usize>,
}

impl CategSeq {
    fn check_categories(seq: &Array1<usize>, n_categories: &Array1<usize>) -> bool {
        for (ix, k) in seq.iter().enumerate() {
            if *k >= n_categories[ix] {
                return false;
            }
        }
        true
    }
    fn new_from(seq: Array1<usize>, len: &usize, n_categories: &Array1<usize>) -> CategSeq {
        if seq.len() != *len {
            panic!("Attempted creation from incorrectly sized vector");
        }
        if !CategSeq::check_categories(&seq, n_categories) {
            panic!("Attempted creation outside category range")
        }
        CategSeq {
            seq: seq,
            len: *len,
            n_categories: n_categories.clone(),
        }
    }
    pub fn new_random(len: &usize, n_categories: &Array1<usize>) -> CategSeq {
        let mut rng = thread_rng();
        let mut seq = Vec::with_capacity(*len);
        let unif = Uniform::new_inclusive(0.0, 1.0);
        for k in 0..*len {
            let samp = unif.sample(&mut rng);
            let samp_int = (samp * (n_categories[k] as f32)).floor() as usize;
            seq.push(samp_int);
        }
        CategSeq::new_from(Array1::from(seq), len, n_categories)
    }
}