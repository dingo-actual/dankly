use std::ops::{Add, Sub, Mul, Div};
use crate::generics::population::Solution;
use rand::seq::SliceRandom;
use rand::thread_rng;
use rand_distr::{Normal, Uniform, Distribution};

//permutation

#[derive(Clone, Copy)]
pub struct Permutation {
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

pub struct SolnIdentPermu {
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
#[derive(Clone, Copy)]
pub struct RealVec {
    len: usize,
    seq: Vec<f32>,
}

impl Add<RealVec> for RealVec {
    type Output = RealVec;

    fn add(self, other: RealVec) -> Option<RealVec> {
        if self.len != other.len {
            panic!("Tried to add incompatible shapes");
        }
        let mut out = Vec::with_capacity(self.len);
        for (x1, x2) in self.seq.iter().zip(other.seq.iter()) {
            out.push(x1 + x2);
        }
        RealVec{len: self.len, seq: out}
    }
}

impl Sub<RealVec> for RealVec {
    type Output = RealVec;

    fn sub(&self, other: &RealVec) -> RealVec {
        if self.len != other.len {
            panic!("Tried to subtract incompatible shapes");
        }
        let mut out = Vec::with_capacity(self.len);
        for (x1, x2) in self.seq.iter().zip(other.seq.iter()) {
            out.push(x1 - x2);
        }
        RealVec{len: self.len, seq: out}
    }
}

impl Mul<RealVec> for RealVec {
    type Output = RealVec;

    fn mul(&self, rhs: &RealVec) -> RealVec {
        if self.len != rhs.len {
            panic!("Tried to multiply incompatible shapes");
        }
        let mut out = Vec::with_capacity(self.len);
        for (x1, x2) in self.seq.iter().zip(rhs.seq.iter()) {
            out.push(x1 * x2);
        }
        Some(RealVec{len: self.len, seq: out})
    }
}

impl Div<RealVec> for RealVec {
    type Output = RealVec;

    fn div(&self, rhs: &RealVec) -> RealVec {
        if self.len != rhs.len {
            panic!("Tried to divide incompatible shapes");
        }
        if rhs.seq.iter().map(|x| x != 0).any() {
            panic!("Division by zero")
        }
        let mut out = Vec::with_capacity(self.len);
        for (x1, x2) in self.seq.iter().zip(rhs.seq.iter()) {
            out.push(x1 / x2);
        }
        RealVec{len: self.len, seq: out}
    }
}

impl RealVec {
    fn new(len: usize) -> RealVec {
        RealVec{len: len, seq: Vec::with_capacity(len)}
    }
    fn new_zeros(len: usize) -> RealVec {
        RealVec{len: len, seq: vec![0; len]}
    }
    fn new_ones(len: usize) -> RealVec {
        RealVec{len: len, seq: vec![1, len]}
    }
    fn new_normal(len: usize, mu: f32, sigma: f32) -> RealVec {
        let rng = thread_rng();
        let normal = Normal::new(mu, sigma).unwrap();
        let out: Vec<f32> = Normal.sample_iter(rng).take(len).collect();
        RealVec{len: len, seq: out}
    }
    fn new_uniform(len: usize, lo: f32, hi: f32) -> RealVec {
        let rng = thread_rng();
        let unif = Uniform::new_inclusive(lo, hi);
        let out: Vec<f32> = Uniform.sample_iter(rng).take(len).collect();
        RealVec{len: len, seq: out}
    }
    fn set(&mut self, other: &RealVec) -> () {
        if self.len != other.len {
            panic!("Tried to set with incompatible shapes");
        }
        for n in 0..self.shape {
            self.seq[n] = other.seq[n].copy();
        }
    }
}
pub struct SolnIdentRealVec {
    gtype: RealVec,
    ptype: RealVec,
    len: usize,
}

impl SolnIdentRealVec {
    pub fn new_ones(len: usize) -> SolnIdentRealVec {
        let ones = RealVec::new_ones(len);
        SolnIdentRealVec {
            gtype: ones.copy(),
            ptype: ones.copy(),
            len: len,
        }
    }
    pub fn new_zeros(len: usize) -> SolnIdentRealVec {
        let zeros = RealVec::new_zeros(len);
        SolnIdentRealVec {
            gtype: zeros.copy(),
            ptype: zeros.copy(),
            len: len,
        }
    }
    pub fn new_normal(len: usize, mu: f32, sigma: f32) -> SolnIdentRealVec {
        let out = RealVec::new_normal(len, mu, sigma);
        SolnIdentRealVec {
            gtype: out.copy(),
            ptype: out.copy(),
            len: len,
        }
    }
    pub fn new_uniform(len: usize, lo: f32, hi: f32) -> SolnIdentRealVec {
        let out = RealVec::new_uniform(len, lo, hi);
        SolnIdentRealVec {
            gtype: out.copy(),
            ptype: out.copy(),
            len: len,
        }
    }
}

impl Solution<RealVec, RealVec> for SolnIdentRealVec {
    fn set_gtype(&mut self, gtype_new: &RealVec) -> () {
        self.gtype.set(gtype_new);
    }
    fn set_ptype(&mut self, ptype_new: &RealVec) -> () {
        self.ptype.set(ptype_new);
    }
    fn get_gtype(&self) -> &RealVec {
        &self.gtype
    }
    fn get_ptype(&self) -> &RealVec {
        &self.ptype
    }
    fn induce_gtype(&mut self) -> () {
        self.gtype = self.ptype.copy();
    }
    fn induce_ptype(&mut self) -> () {
        self.ptype = self.gtype.copy();
    }
}

//bounded real vector

#[derive(Clone, Copy)]
pub struct BoundedRealVec {
    len: usize,
    seq: Vec<f32>,
    ubounds: Vec<f32>,
    lbounds: Vec<f32>,
}

pub struct SolnIdentBoundedRealVec {
    gtype: BoundedRealVec,
    ptype: BoundedRealVec,
    len: usize,
    ubounds: Vec<f32>,
    lbounds: Vec<f32>,
}

impl Add<BoundedRealVec> for BoundedRealVec {
    type Output = BoundedRealVec;

    fn add(&self, other: &BoundedRealVec) -> BoundedRealVec {
        if !self.check_compatible(other) {
            panic!("Tried to add incompatible bounded vectors")
        }
        let mut out = Vec::with_capacity(self.len);
        for (x1, x2) in self.seq.iter().zip(other.seq.iter()) {
            out.push(x1 + x2);
        }
        self.clip(&mut out, &self.lbounds, &self.hbounds);
        BoundedRealVec{
            len: self.len, 
            seq: out,
            lbounds: *self.lbounds.copy(),
            hbounds: *self.hbounds.copy(),
        }
    }
}

impl Sub<BoundedRealVec> for BoundedRealVec {
    type Output = BoundedRealVec;

    fn sub(&self, other: &BoundedRealVec) -> BoundedRealVec {
        if !self.check_compatible(other) {
            panic!("Tried to subtract incompatible bounded vectors")
        }
        let mut out = Vec::with_capacity(self.len);
        for (x1, x2) in self.seq.iter().zip(other.seq.iter()) {
            out.push(x1 - x2);
        }
        self.clip(&mut out, &self.lbounds, &self.hbounds);
        BoundedRealVec{
            len: self.len, 
            seq: out,
            lbounds: *self.lbounds.copy(),
            hbounds: *self.hbounds.copy(),
        }
    }
}

impl Mul<BoundedRealVec> for BoundedRealVec {
    type Output = BoundedRealVec;

    fn mul(&self, rhs: &BoundedRealVec) -> BoundedRealVec {
        if !self.check_compatible(rhs) {
            panic!("Tried to multiply incompatible bounded vectors")
        }
        let mut out = Vec::with_capacity(self.len);
        for (x1, x2) in self.seq.iter().zip(rhs.seq.iter()) {
            out.push(x1 * x2);
        }
        self.clip(&mut out, &self.lbounds, &self.hbounds);
        BoundedRealVec{
            len: self.len, 
            seq: out,
            lbounds: *self.lbounds.copy(),
            hbounds: *self.hbounds.copy(),
        }
    }
}

impl Div<BoundedRealVec> for BoundedRealVec {
    type Output = BoundedRealVec;

    fn div(&self, rhs: &BoundedRealVec) -> BoundedRealVec {
        if !self.check_compatible(rhs) {
            panic!("Tried to divide incompatible bounded vectors")
        }
        if rhs.seq.iter().map(|x| x != 0).any() {
            panic!("Division by zero")
        }
        let mut out = Vec::with_capacity(self.len);
        for (x1, x2) in self.seq.iter().zip(rhs.seq.iter()) {
            out.push(x1 * x2);
        }
        self.clip(&mut out, &self.lbounds, &self.hbounds);
        BoundedRealVec{
            len: self.len, 
            seq: out,
            lbounds: *self.lbounds.copy(),
            hbounds: *self.hbounds.copy(),
        }
    }
}


impl BoundedRealVec {
    fn check_compatible(&self, other: &BoundedRealVec) -> bool {
        if self.len != other.len {
            return False
        }
        for (lo1, lo2) in self.lbounds.iter().zip(other.lbounds.iter()) {
            if lo1 != lo2 {
                return False
            }
        }
        for (hi1, hi2) in self.hbounds.iter().zip(other.hbounds.iter()) {
            if hi1 != hi2 {
                return False
            }
        }
        True
    }
    fn check_bounds(len: usize, lbounds: &Vec<f32>, ubounds: &Vec<f32>) -> bool {
        if ubounds.len() != len || lbounds.len() != len {
            return False
        }
        for (lo, hi) in lbounds.iter().zip(ubounds.iter()) {
            if lo > hi {
                return False
            }
        }
        True
    }
    fn clip(seq: &mut Vec<f32>, lbounds: &Vec<f32>, ubounds: &Vec<f32>) -> () {
        for n in 0..seq.len() {
            if seq[n] < lbounds[n] {
                seq[n] = lbounds[n];
            } else if seq[n] > ubounds[n] {
                seq[n] = ubounds[n];
            }
        }
    }
    fn bounded_normal(lbound: &f32, ubound: &f32) -> f32 {
        let rng = thread_rng();
        let mu = (ubound + lbound) / 2;
        let sigma = (hbound - lbound) / 4;
        let normal = Normal::new(mu, sigma).unwrap();
        let mut samp = normal.sample(&mut rng);
        while samp < lbound || samp > ubound {
            samp = normal.sample(&mut rng);
        }
        samp;

    }
    pub fn new(len: usize, lbounds: &Vec<f32>, ubounds: &Vec<f32>) -> RealVec {
        if !self.check_bounds(len, lbounds, ubounds) {
            panic!("Invalid bounds");
        }
        BoundedRealVec{
            len: len, 
            seq: Vec::with_capacity(len),
            ubounds: *ubounds.copy(),
            lbounds: *lbounds.copy(),
        }
    }
    pub fn new_lo(len: usize, lbounds: &Vec<f32>, ubounds: &Vec<f32>) -> BoundedRealVec {
        if !self.check_bounds(len, lbounds, ubounds) {
            panic!("Invalid bounds");
        }
        BoundedRealVec{
            len: len, 
            seq: *lbounds.copy(), 
            ubounds: *ubounds.copy(),
            lbounds: *lbounds.copy(),
        }
    }
    pub fn new_hi(len: usize, lbounds: &Vec<f32>, ubounds: &Vec<f32>) -> BoundedRealVec {
        if !self.check_bounds(len, lbounds, ubounds) {
            panic!("Invalid bounds");
        }
        BoundedRealVec{
            len: len, 
            seq: *hbounds.copy(), 
            ubounds: *ubounds.copy(),
            lbounds: *lbounds.copy(),
        }
    }
    pub fn new_normal(len: usize, lbounds: &Vec<f32>, ubounds: &Vec<f32>) -> BoundedRealVec {
        if ! self.check_bounds(len, lbounds, ubounds) {
            panic!("Invalid bounds");
        }
        let mut out = Vec::with_capacity(len);
        for n in 0..len {
            out.push(self.bounded_normal(&lbounds[n], &ubounds[n]))
        }
        BoundedRealVec{
            len: len, 
            seq: out, 
            ubounds: *ubounds.copy(),
            lbounds: *lbounds.copy(),
        }
    }
    pub fn new_uniform(len: usize, lbounds: &Vec<f32>, ubounds: &Vec<f32>) -> RealVec {
        if !self.check_bounds(len, lbounds, ubounds) {
            panic!("Invalid bounds");
        }
        let rng = thread_rng();
        let unif = Uniform::new_inclusive(0, 1);
        let mut out = Vec::with_capacity(len);
        for n in 0..len {
            let samp = unif.sample(&mut rng);
            let scaled_samp = lbounds[n] + samp * (ubounds[n] - lbounds[n]);
            out.push(scaled_samp);
        }
        BoundedRealVec{
            len: len,
            seq: out,
            ubounds: *ubounds.copy(),
            lbounds: *lbounds.copy(),
        }
    }
    pub fn set(&mut self, other: &BoundedRealVec) -> () {
        if self.check_compatible(other) {
            panic!("Set with incompatible shapes");
        }
        for n in 0..self.shape {
            self.seq[n] = other.seq[n].copy();
        }
        self.clip(&mut self.seq, &self.lbounds, &self.ubounds);
    }
}

impl Solution<BoundedRealVec, BoundedRealVec> for SolnIdentBoundedRealVec {
    fn set_gtype(&mut self, gtype_new: &BoundedRealVec) -> () {
        self.gtype.set(gtype_new);
    }
    fn set_ptype(&mut self, ptype_new: &BoundedRealVec) -> () {
        self.ptype.set(ptype_new);
    }
    fn get_gtype(&self) -> &BoundedRealVec {
        &self.gtype
    }
    fn get_ptype(&self) -> &BoundedRealVec {
        &self.ptype
    }
    fn induce_gtype(&mut self) -> () {
        self.gtype = self.ptype.copy();
    }
    fn induce_ptype(&mut self) -> () {
        self.ptype = self.gtype.copy();
    }
}

//binary string

//categorical sequence