use num::Integer;
use rand::seq::SliceRandom;
use rand::thread_rng;
use rand_distr::{Uniform, Distribution};

// util functions

fn binsrch_probs(cum_prob: &Vec<f32>, target: &f32) -> usize {
    if *target < 0.0 || *target > 1.0 {
        panic!("Target for binary search must be between 0 and 1");
    }
    let mut lo = 0;
    let mut hi = cum_prob.len() - 1;
    let mut mid = (hi + lo) / 2;
    while !(cum_prob[mid] < *target && cum_prob[mid + 1] > *target) {
        if cum_prob[mid] > *target {
            hi = mid;
        } else if cum_prob[mid + 1] < *target {
            lo = mid;
        }
        mid = (hi + lo) / 2;
    }
    mid
}

fn choose_byprob(probs: &Vec<f32>, n: &usize, m: &usize) -> Vec<Vec<usize>> {
    let rng = thread_rng();
    let unif = Uniform::new(0.0, 1.0);
    let mut cum_prob: Vec<f32> = Vec::with_capacity(probs.len() + 1);
    let mut out: Vec<Vec<usize>> = Vec::with_capacity(*n);
    cum_prob.push(0.0);
    for (ix, x) in probs.iter().enumerate() {
        cum_prob.push(x + cum_prob[ix]);
    }
    for _ in 0..*n {
        let mut batch = Vec::with_capacity(*m);
        let samples: Vec<f32> = unif.sample_iter(rng).take(*m).collect();
        for samp in samples.iter() {
            batch.push(binsrch_probs(&cum_prob, &samp));
        }
        out.push(batch);
    }
    out
}

// choose function for survival, crossover and local search

pub fn choose_rand_byfitness(fitneseses: &Vec<f32>, n_batches: &usize, m_per_batch: &usize) -> Vec<Vec<usize>> {
    if *n_batches == 0 || *m_per_batch == 0 {
        panic!("Cannot choose zero members of population");
    }
    let mut sum_fit = 0.0;
    for x in fitneseses.iter() {
        sum_fit += *x;
    }
    let mut probs: Vec<f32> = Vec::with_capacity(fitneseses.len());
    for x in fitneseses.iter() {
        probs.push(*x / sum_fit);
    }
    choose_byprob(&probs, n_batches, m_per_batch)
}

// choose function for survival, crossover and local search

pub fn choose_rand_byloss(losses: &Vec<f32>, n_batches: &usize, m_per_batch: &usize) -> Vec<Vec<usize>> {
    if *n_batches == 0 || *m_per_batch == 0 {
        panic!("Cannot choose zero members of population");
    }
    if losses.iter().any(|x| *x <= 0.0) {
        panic!("Losses must by strictly positive");
    }
    let losses_inv = losses.iter().map(|x| 1.0 / *x).collect();
    choose_rand_byfitness(&losses_inv, n_batches, m_per_batch)
}

fn tournament(fitnesses: &Vec<f32>, ixs: &mut Vec<usize>, m_per_batch: &usize) {
    let (pairs, rem) = ixs.len().div_mod_floor(&2);
    if pairs + rem < *m_per_batch {
        ixs.sort_by(|a, b| fitnesses[*b].partial_cmp(&fitnesses[*a]).unwrap());
        while ixs.len() > *m_per_batch {
            let _ = ixs.pop();
        }
        return;
    } 
    let mut remove_ixs = Vec::with_capacity(pairs);
    let unif = Uniform::new(0.0, 1.0);
    let mut rng = thread_rng();
    for p in 0..pairs {
        let ix1 = 2 * p;
        let ix2 = 2 * p + 1;
        let prob_ix1 = fitnesses[ixs[ix1]] / (fitnesses[ixs[ix1]] + fitnesses[ixs[ix2]]);
        let samp = unif.sample(&mut rng);
        if samp < prob_ix1 {
            remove_ixs.push(ix2);
        } else {
            remove_ixs.push(ix1);
        }
    }
    remove_ixs.sort();
    remove_ixs.reverse();
    for ix in remove_ixs {
        ixs.remove(ix);
    }
    tournament(fitnesses, ixs, m_per_batch)
}

// choose function for survival, crossover and local search

pub fn choose_tournament_byfitness(fitneseses: &Vec<f32>, n_batches: &usize, m_per_batch: &usize) -> Vec<Vec<usize>> {
    if *n_batches == 0 || *m_per_batch == 0 {
        panic!("Cannot choose zero members of population");
    }
    if fitneseses.len() < *m_per_batch {
        panic!("Cannot choose more than entire population via tournament");
    }
    let mut out = Vec::with_capacity(*n_batches);
    for _ in 0..*n_batches {
        let mut ixs: Vec<usize> = (0..fitneseses.len()).collect();
        ixs.shuffle(&mut thread_rng());
        tournament(fitneseses, &mut ixs, m_per_batch);
        out.push(ixs);
    }
    out
}

// choose function for survival, crossover and local search

pub fn choose_tournament_byloss(losses: &Vec<f32>, n_batches: &usize, m_per_batch: &usize) -> Vec<Vec<usize>> {
    if losses.iter().any(|x| *x <= 0.0) {
        panic!("Losses must by strictly positive");
    }
    let losses_inv = losses.iter().map(|x| 1.0 / *x).collect();
    choose_tournament_byfitness(&losses_inv, n_batches, m_per_batch)
}

// choose function for survival

pub fn choose_topk_byfitness(fitneseses: &Vec<f32>, k: &usize) -> Vec<Vec<usize>> {
    if *k > fitneseses.len() {
        panic!("Cannot choose more than entire population via top k");
    }
    let mut ixs: Vec<usize> = (0..fitneseses.len()).collect();
    ixs.sort_by(|a, b| fitneseses[*b].partial_cmp(&fitneseses[*a]).unwrap());
    while ixs.len() > *k {
        let _ = ixs.pop();
    }
    vec![ixs]
}

// choose function for survival

pub fn choose_topk_byloss(losses: &Vec<f32>, k: &usize) -> Vec<Vec<usize>> {
    if *k > losses.len() {
        panic!("Cannot choose more than entire population via top k");
    }
    let mut ixs: Vec<usize> = (0..losses.len()).collect();
    ixs.sort_by(|a, b| losses[*a].partial_cmp(&losses[*b]).unwrap());
    while ixs.len() > *k {
        let _ = ixs.pop();
    }
    vec![ixs]
}

// choose function factory for survival

pub fn elitism_factory(m_per_batch: &usize, k_elite: &'static usize, other_chooser: Box<dyn Fn(&Vec<f32>, &usize, &usize) -> Vec<Vec<usize>>>) -> Box<dyn Fn(&Vec<f32>, &usize, &usize) -> Vec<Vec<usize>>>{
    if k_elite > m_per_batch {
        panic!("Cannot have more elite solutions than total number of solutions");
    }
    Box::new(move |f: &Vec<f32>, _n: &usize, m: &usize| choose_interrupt(f, m, k_elite, &other_chooser))
}

fn choose_interrupt(fitnesses: &Vec<f32>, m_per_batch: &usize, k_elite: &usize, other_chooser: &Box<dyn Fn(&Vec<f32>, &usize, &usize) -> Vec<Vec<usize>>>) -> Vec<Vec<usize>> {
    let mut top_k = choose_topk_byfitness(fitnesses, k_elite)[0].clone();
    let mut other_choice = other_chooser(fitnesses, &1, &(*m_per_batch - k_elite))[0].clone();
    let mut out: Vec<Vec<usize>> = Vec::new();
    top_k.append(&mut other_choice);
    out.push(top_k);
    out
}

// TODO: implement Selector creation functions