use std::ops::{Add, Sub, Mul, Div, BitAnd, BitOr, BitXor, Not};
use crate::generics::population::{Solution, Genotype, Phenotype};
use rand::seq::SliceRandom;
use rand::thread_rng;
use rand::distributions::{Bernoulli, Distribution};
use rand_distr::{Normal, Uniform, Distribution};

//permutation

#[derive(Clone)]
pub struct Permutation {
    seq: Vec<usize>,
    len: usize,
}

pub struct SolnIdentPermu {
    gtype: Permutation,
    ptype: Permutation,
    len: usize,
}

impl Genotype<Vec<usize>> for Permutation {
    fn get(&self) -> &Vec<usize> {
        &self.seq
    }
    fn set(&self, other: &Vec<usize>) {
        if other.len() != self.len {
            panic!("Attempted set with mismatched shapes")
        }
        for n in 0..self.len {
            self.seq[n] = other[n].clone();
        }
    }
}

impl Phenotype<Vec<usize>> for Permutation {
    fn get(&self) -> &Vec<usize> {
        &self.seq
    }
    fn set(&self, other: &Vec<usize>) {
        if other.len() != self.len {
            panic!("Attempted set with mismatched shapes")
        }
        for n in 0..self.len {
            self.seq[n] = other[n].clone();
        }
    }
}

//TODO: change Result return type to panic instead

impl Permutation {
    fn new_random(len: &usize) -> Permutation {
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
    fn get_len(&self) {
        self.len
    }
    fn swap(&mut self, ix1: &usize, ix2: &usize) -> Result<&usize, &str> {
        if self.check_bounds(ix1) && self.check_bounds(ix2) {
            if ix1 != ix2 {
                x1 = self.seq[*ix1].clone();
                self.seq[*ix1] = self.seq[*ix2].clone();
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

impl SolnIdentPermu {
    pub fn new_random(len: &usize) -> SolnIdentPermu {
        let perm = Permutation::new_random(len);
        SolnIdentPermu {
            gtype: perm,
            ptype: perm.clone(),
            len: *len,
        }
    }
}

impl Solution<Permutation, Permutation> for SolnIdentPermu {
    fn set_gtype(&mut self, gtype_new: &Permutation) -> () {
        self.gtype.set(gtype_new.get());
    }
    fn set_ptype(&mut self, ptype_new: &Permutation) -> () {
        self.ptype.set(ptype_new.get());
    }
    fn get_gtype(&self) -> &Permutation {
        &self.gtype
    }
    fn get_ptype(&self) -> &Permutation {
        &self.ptype
    }
    fn get_gtype_mut(&mut self) -> &mut Permutation {
        &mut self.gtype
    }
    fn get_ptype_mut(&mut self) -> &mut Permutation {
        &mut self.ptype
    }
    fn induce_gtype(&mut self) -> () {
        self.gtype.set(self.ptype.get().clone());
    }
    fn induce_ptype(&mut self) -> () {
        self.ptype.set(self.gtype.get().clone());
    }
}


//TODO: fix Solution implementations
//TODO: copy -> clone

//real vector
#[derive(Clone, Copy)]
pub struct RealVec {
    len: usize,
    seq: Vec<f32>,
}

pub struct SolnIdentRealVec {
    gtype: RealVec,
    ptype: RealVec,
    len: usize,
}

impl Genotype for RealVec {
    fn get(&self) -> &Vec<f32> {
        &self.seq
    }
    fn set(&self, other: &Vec<f32>) {
        if other.len() != self.len {
            panic!("Attempted set with mismatched shapes")
        }
        for n in 0..self.len {
            self.seq[n] = other[n].copy();
        }
    }
}

impl Phenotype for RealVec {
    fn get(&self) -> &Vec<f32> {
        &self.seq
    }
    fn set(&self, other: &Vec<f32>) {
        if other.len() != self.len {
            panic!("Attempted set with mismatched shapes")
        }
        for n in 0..self.len {
            self.seq[n] = other[n].copy();
        }
    }
}

impl Add<&RealVec> for RealVec {
    type Output = RealVec;

    fn add(self, other: &RealVec) -> RealVec {
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

impl Sub<&RealVec> for RealVec {
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

impl Mul<&RealVec> for RealVec {
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

impl Div<&RealVec> for RealVec {
    type Output = RealVec;

    fn div(&self, rhs: &RealVec) -> RealVec {
        if self.len != rhs.len {
            panic!("Tried to divide incompatible shapes");
        }
        if rhs.seq.iter().map(|x| x == 0).any() {
            panic!("Division by zero");
        }
        let mut out = Vec::with_capacity(self.len);
        for (x1, x2) in self.seq.iter().zip(rhs.seq.iter()) {
            out.push(x1 / x2);
        }
        RealVec{len: self.len, seq: out}
    }
}

impl RealVec {
    fn new_from(len: &usize, seq: Vec<f32>) -> RealVec {
        if len != seq.len() {
            panic!("Incorrect initial sequence length");
        }
        RealVec{len: len, seq: seq}
    }
    pub fn new_zeros(len: &usize) -> RealVec {
        let seq = vec![0; len];
        RealVec::new_from(len, seq)
    }
    pub fn new_ones(len: &usize) -> RealVec {
        let seq = vec![1; len];
        RealVec::new_from(len, seq)
    }
    pub fn new_normal(len: &usize, mu: &f32, sigma: &f32) -> RealVec {
        let rng = thread_rng();
        let normal = Normal::new(mu, sigma).unwrap();
        let seq: Vec<f32> = Normal.sample_iter(rng).take(len).collect();
        RealVec::new_from(len, seq)
    }
    pub fn new_uniform(len: usize, lo: f32, hi: f32) -> RealVec {
        let rng = thread_rng();
        let unif = Uniform::new_inclusive(lo, hi);
        let seq: Vec<f32> = Uniform.sample_iter(rng).take(len).collect();
        RealVec::new_from(len, seq)
    }
}

impl SolnIdentRealVec {
    fn new_from(len: &usize, seq: RealVec) {
        SolnIdentRealVec {
            gtype: seq.copy(),
            ptype: seq.copy(),
            len: len,
        }
    }
    pub fn new_ones(len: usize) -> SolnIdentRealVec {
        let seq = RealVec::new_ones(len);
        SolnIdentRealVec::new_from(len, seq)
    }
    pub fn new_zeros(len: usize) -> SolnIdentRealVec {
        let seq = RealVec::new_zeros(len);
        SolnIdentRealVec::new_from(len, seq)
    }
    pub fn new_normal(len: usize, mu: f32, sigma: f32) -> SolnIdentRealVec {
        let seq = RealVec::new_normal(len, mu, sigma);
        SolnIdentRealVec::new_from(len, seq)
    }
    pub fn new_uniform(len: usize, lo: f32, hi: f32) -> SolnIdentRealVec {
        let seq = RealVec::new_uniform(len, lo, hi);
        SolnIdentRealVec::new_from(len, seq)
    }
}

impl Solution for SolnIdentRealVec {
    fn set_gtype(&mut self, gtype_new: &Vec<f32>) -> () {
        self.gtype.set(gtype_new);
    }
    fn set_ptype(&mut self, ptype_new: &Vec<f32>) -> () {
        self.ptype.set(ptype_new);
    }
    fn get_gtype(&self) -> &Vec<f32> {
        &self.gtype.get()
    }
    fn get_ptype(&self) -> &Vec<f32> {
        &self.ptype.get()
    }
    fn induce_gtype(&mut self) -> () {
        self.gtype.set(self.ptype.get().copy());
    }
    fn induce_ptype(&mut self) -> () {
        self.ptype.set(self.gtype.get().copy());
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

impl Genotype for BoundedRealVec {
    fn get(&self) -> &Vec<f32> {
        &self.seq
    }
    fn set(&self, other: &Vec<f32>) {
        if other.len() != self.len {
            panic!("Attempted set with mismatched shapes")
        }
        if !BoundedRealVec::within_bounds(other, &self.lbounds, &self.ubounds) {
            panic!{"Attempted set outside of bounds"}
        }
        for n in 0..self.len {
            self.seq[n] = other[n].copy();
        }
    }
}

impl Phenotype for BoundedRealVec {
    fn get(&self) -> &Vec<f32> {
        &self.seq
    }
    fn set(&self, other: &Vec<f32>) {
        if other.len() != self.len {
            panic!("Attempted set with mismatched shapes")
        }
        if !BoundedRealVec::within_bounds(other, &self.lbounds, &self.ubounds) {
            panic!{"Attempted set outside of bounds"}
        }
        for n in 0..self.len {
            self.seq[n] = other[n].copy();
        }
    }
}

impl Add<&BoundedRealVec> for BoundedRealVec {
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

impl Sub<&BoundedRealVec> for BoundedRealVec {
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

impl Mul<&BoundedRealVec> for BoundedRealVec {
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

impl Div<&BoundedRealVec> for BoundedRealVec {
    type Output = BoundedRealVec;

    fn div(&self, rhs: &BoundedRealVec) -> BoundedRealVec {
        if !self.check_compatible(rhs) {
            panic!("Tried to divide incompatible bounded vectors")
        }
        if rhs.seq.iter().map(|x| x == 0).any() {
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
            return false
        }
        for (lo1, lo2) in self.lbounds.iter().zip(other.lbounds.iter()) {
            if lo1 != lo2 {
                return false
            }
        }
        for (hi1, hi2) in self.hbounds.iter().zip(other.hbounds.iter()) {
            if hi1 != hi2 {
                return false
            }
        }
        true
    }
    fn check_bounds(len: &usize, lbounds: &Vec<f32>, ubounds: &Vec<f32>) -> bool {
        if ubounds.len() != len || lbounds.len() != len {
            return false
        }
        for (lo, hi) in lbounds.iter().zip(ubounds.iter()) {
            if lo > hi {
                return false
            }
        }
        true
    }
    fn within_bounds(seq: &Vec<f32>, lbounds: &Vec<f32>, ubounds: &Vec<f32>) -> bool {
        if seq.len() != lbounds.len() || seq.len() != ubounds.len() {
            return false;
        }
        for (n, x) in seq.iter().enumerate() {
            if x < lbounds[n] || x > ubounds[n] {
                return false;
            }
        }
        true
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
    fn new_from(len: &usize, lbounds: &Vec<f32>, ubounds: &Vec<f32>, seq: Vec<f32>) -> BoundedRealVec {
        if !BoundedRealVec::check_bounds(len, lbounds, ubounds) {
            panic!("Invalid bounds");
        }
        if !BoundedRealVec::within_bounds(&seq, lbounds, ubounds) {
            panic!("Attempted creation with out of bounds values")
        }
        BoundedRealVec{
            len: len, 
            seq: seq,
            ubounds: *ubounds.copy(),
            lbounds: *lbounds.copy(),
        }
    }
    pub fn new_lo(len: &usize, lbounds: &Vec<f32>, ubounds: &Vec<f32>) -> BoundedRealVec {
        BoundedRealVec::new_from(len, lbounds, ubounds, *lbounds.copy())
    }
    pub fn new_hi(len: usize, lbounds: &Vec<f32>, ubounds: &Vec<f32>) -> BoundedRealVec {
        BoundedRealVec::new_from(len, lbounds, ubounds, *ubounds.copy())
    }
    pub fn new_normal(len: usize, lbounds: &Vec<f32>, ubounds: &Vec<f32>) -> BoundedRealVec {
        let mut seq = Vec::with_capacity(len);
        for n in 0..len {
            seq.push(self.bounded_normal(&lbounds[n], &ubounds[n]))
        }
        BoundedRealVec::new_from(len, lbounds, ubounds, seq)
    }
    pub fn new_uniform(len: usize, lbounds: &Vec<f32>, ubounds: &Vec<f32>) -> RealVec {
        let rng = thread_rng();
        let unif = Uniform::new_inclusive(0, 1);
        let mut seq = Vec::with_capacity(len);
        for n in 0..len {
            let samp = unif.sample(&mut rng);
            let scaled_samp = lbounds[n] + samp * (ubounds[n] - lbounds[n]);
            seq.push(scaled_samp);
        }
        BoundedRealVec::new_from(len, lbounds, ubounds, seq)
    }
}

impl SolnIdentBoundedRealVec {
    fn new_from(len: &usize, lbounds: &Vec<f32>, ubounds: &Vec<f32>, realvec: &BoundedRealVec) -> SolnIdentBoundedRealVec {
        SolnIdentBoundedRealVec {
            gtype: *realvec.copy(),
            ptype: *realvec.copy(),
            len: len,
            ubounds: *ubounds.copy(),
            lbounds: *lbounds.copy(),
        }
    }
    pub fn new_lo(len: &usize, lbounds: &Vec<f32>, ubounds: &Vec<f32>) -> SolnIdentBoundedRealVec {
        let realvec = BoundedRealVec::new_lo(len, lbounds, ubounds);
        SolnIdentBoundedRealVec::new_from(len, lbounds, ubounds, &realvec)
    }
    pub fn new_hi(len: &usize, lbounds: &Vec<f32>, ubounds: &Vec<f32>) -> SolnIdentBoundedRealVec {
        let realvec = BoundedRealVec::new_hi(len, lbounds, ubounds);
        SolnIdentBoundedRealVec::new_from(len, lbounds, ubounds, &realvec)
    }
    pub fn new_normal(len: &usize, lbounds: &Vec<f32>, ubounds: &Vec<f32>) -> SolnIdentBoundedRealVec {
        let realvec = BoundedRealVec::new_normal(len, lbounds, ubounds);
        SolnIdentBoundedRealVec::new_from(len, lbounds, ubounds, &realvec)
    }
    pub fn new_uniform(len: &usize, lbounds: &Vec<f32>, ubounds: &Vec<f32>) -> SolnIdentBoundedRealVec {
        let realvec = BoundedRealVec::new_uniform(len, lbounds, ubounds);
        SolnIdentBoundedRealVec::new_from(len, lbounds, ubounds, &realvec)
    }
}

impl Solution for SolnIdentBoundedRealVec {
    fn set_gtype(&mut self, gtype_new: &Vec<f32>) -> () {
        self.gtype.set(gtype_new);
    }
    fn set_ptype(&mut self, ptype_new: &Vec<f32>) -> () {
        self.ptype.set(ptype_new);
    }
    fn get_gtype(&self) -> &Vec<f32> {
        &self.gtype.get()
    }
    fn get_ptype(&self) -> &Vec<f32> {
        &self.ptype.get()
    }
    fn induce_gtype(&mut self) -> () {
        self.gtype.set(self.ptype.get().copy());
    }
    fn induce_ptype(&mut self) -> () {
        self.ptype.set(self.gtype.get().copy());
    }
}

//binary string

#[derive(Clone, Copy)]
struct BinarySeq {
    len: usize,
    seq: Vec<bool>,
}

pub struct SolnIdentBinarySeq {
    gtype: BinarySeq,
    ptype: BinarySeq,
    len: usize,
}

impl BinarySeq {
    fn check_compatible(&self, other: &BinarySeq) {
        return self.len == other.len;
    }
    fn new_from(len: &usize, seq: Vec<bool>) {
        BinarySeq {
            len: len,
            seq: seq,
        }
    }
    pub fn new_true(len: &usize) -> BinarySeq {
        let seq = vec![true, len];
        BinarySeq::new_from(len, seq)
    }
    pub fn new_false(len: &usize) -> BinarySeq {
        let seq = vec![false, len];
        BinarySeq::new_from(len, seq)
    }
    pub fn new_random(len: &usize, p: &f32) {
        if p < 0 || p > 1 {
            panic!("Invalid probability")
        }
        let rng = rand::thread_rng();
        let bern = Bernoulli::new(p).unwrap();
        let seq = bern.sample_iter(rng).take(len).collect();
        BinarySeq::new_from(len, seq)
    }
    pub fn parity(&self) -> bool {
        self.seq.iter().fold(false, |acc, x| acc ^ x)
    }
    pub fn count(&self) -> usize {
        self.seq.iter().filter(|&x| x).count()
    }
}

impl BitAnd for BinarySeq {
    type Output = Self;

    fn bitand(&self, rhs: &BinarySeq) -> Self::Output {
        if !self.check_compatible(rhs) {
            panic!("Attempted AND on incompatible shapes")
        }
        let mut out = Vec::with_capacity(self.len);
        for (x1, x2) in self.seq.iter().zip(other.seq.iter()) {
            out.push(x1 && x2);
        }
        BinarySeq {
            len: self.len,
            seq: out,
        }
    }
}

impl BitOr for BinarySeq {
    type Output = Self;

    fn bitor(&self, rhs: &BinarySeq) -> Self::Output {
        if !self.check_compatible(rhs) {
            panic!("Attempted OR on incompatible shapes")
        }
        let mut out = Vec::with_capacity(self.len);
        for (x1, x2) in self.seq.iter().zip(other.seq.iter()) {
            out.push(x1 || x2);
        }
        BinarySeq {
            len: self.len,
            seq: out,
        }
    }
}

impl BitXor for BinarySeq {
    type Output = Self;

    fn bitxor(&self, rhs: &BinarySeq) -> Self::Output {
        if !self.check_compatible(rhs) {
            panic!("Attempted XOR on incompatible shapes")
        }
        let mut out = Vec::with_capacity(self.len);
        for (x1, x2) in self.seq.iter().zip(other.seq.iter()) {
            out.push(x1 ^ x2);
        }
        BinarySeq {
            len: self.len,
            seq: out,
        }
    }
}

impl SolnIdentBinarySeq {
    fn new_from(len: &usize, seq: BinarySeq) {
        SolnIdentBinarySeq {
            len: len,
            gtype: seq.copy(),
            ptype: seq.copy(),
        }
    }
    pub fn new_true(len: &usize) -> SolnIdentBinarySeq {
        let seq = BinarySeq::new_true(len);
        SolnIdentBinarySeq::new_from(len, seq)
    }
    pub fn new_false(len: &usize) -> SolnIdentBinarySeq {
        let seq = BinarySeq::new_false(len);
        SolnIdentBinarySeq::new_from(len, seq)
    }
    pub fn new_random(len: &usize, p: &f32) -> SolnIdentBinarySeq {
        let seq = BinarySeq::new_random(len, p);
        SolnIdentBinarySeq::new_from(len, seq)
    }
}

impl Solution for SolnIdentRealVec {
    fn set_gtype(&mut self, gtype_new: &Vec<bool>) -> () {
        self.gtype.set(gtype_new);
    }
    fn set_ptype(&mut self, ptype_new: &Vec<bool>) -> () {
        self.ptype.set(ptype_new);
    }
    fn get_gtype(&self) -> &Vec<bool> {
        &self.gtype.get()
    }
    fn get_ptype(&self) -> &Vec<bool> {
        &self.ptype.get()
    }
    fn induce_gtype(&mut self) -> () {
        self.gtype.set(self.ptype.get().copy());
    }
    fn induce_ptype(&mut self) -> () {
        self.ptype.set(self.gtype.get().copy());
    }
}

//categorical sequence

pub struct CategSeq {
    seq: Vec<usize>,
    len: usize,
    n_categories: Vec<usize>,
}

pub struct SolnIdentCategSeq {
    gtype: CategSeq,
    ptype: CategSeq,
    len: usize,
    n_categories: Vec<usize>,
}

impl Genotype for CategSeq {
    fn get(&self) -> &Vec<usize> {
        &self.seq
    }
    fn set(&self, other: &Vec<usize>) {
        if other.len() != self.len {
            panic!("Attempted set with mismatched shapes");
        }
        if !CategSeq::check_categories(other, &self.n_categories) {
            panic!("Attempted set outside category range");
        }
        for n in 0..self.len {
            self.seq[n] = other[n].copy();
        }
    }
}

impl Phenotype for CategSeq {
    fn get(&self) -> &Vec<usize> {
        &self.seq
    }
    fn set(&self, other: &Vec<usize>) {
        if other.len() != self.len {
            panic!("Attempted set with mismatched shapes");
        }
        if !CategSeq::check_categories(other, &self.n_categories) {
            panic!("Attempted set outside category range");
        }
        for n in 0..self.len {
            self.seq[n] = other[n].copy();
        }
    }
}

impl CategSeq {
    fn check_categories(seq: &Vec<usize>, n_categories: &Vec<usize>) -> bool {
        for k in seq.iter() {
            if k < 0 || k >= n_categories[k] {
                return false;
            }
        }
        true
    }
    fn new_from(seq: Vec<usize>, len: &usize, n_categories: &Vec<usize>) -> CategSeq {
        if seq.len() != len {
            panic!("Attempted creation from incorrectly sized vector");
        }
        if !CategSeq::check_categories(&seq, n_categories) {
            panic!("Attempted creation outside category range")
        }
        CategSeq {
            seq: seq,
            len: len,
            n_categories: n_categories,
        }
    }
    pub fn new_random(len: &usize, n_categories: &Vec<usize>) {
        let mut rng = thread_rng();
        let mut seq = Vec::with_capacity(len);
        let unif = Uniform::new(0, 1);
        for k in 0..len {
            let samp = unif.sample(&mut rng);
            let samp_int = (samp * n_categories[k]).floor() as usize;
            seq.push(samp_int);
        }
        CategSeq::new_from(seq, len, n_categories)
    }
}

impl SolnIdentCategSeq {
    fn new_from(seq: CategSeq, len: &usize, n_categories: &Vec<usize>) {
        SolnIdentCategSeq {
            gtype: seq.copy(),
            ptype: seq.copy(),
            len: len,
            n_categories: n_categories.copy(),
        }
    }
    pub fn new_random(len: &usize, n_categories: &Vec<usize>) {
        let seq = CategSeq::new_random(len, n_categories);
        SolnIdentCategSeq::new_from(seq, len, n_categories)
    }
}

impl Solution for SolnIdentCategSeq {
    fn set_gtype(&mut self, gtype_new: &Vec<usize>) -> () {
        self.gtype.set(gtype_new);
    }
    fn set_ptype(&mut self, ptype_new: &Vec<usize>) -> () {
        self.ptype.set(ptype_new);
    }
    fn get_gtype(&self) -> &Vec<usize> {
        &self.gtype.get()
    }
    fn get_ptype(&self) -> &Vec<usize> {
        &self.ptype.get()
    }
    fn induce_gtype(&mut self) -> () {
        self.gtype.set(self.ptype.get().copy());
    }
    fn induce_ptype(&mut self) -> () {
        self.ptype.set(self.gtype.get().copy());
    }
}