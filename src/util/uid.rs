use rand::Rng;

pub struct UidStack {
    uids: Vec<&String>,
    id_len: usize,
    size: usize,
    alphabet: String,
}

impl UidStack {
    pub fn new(id_len: &usize, size: &usize, alph: &String) -> UidStack {
        UidStack{
            uids: Vec::new(),
            id_len: *id_len,
            size: *n,
            alphabet: *alph,
        }
    }
    pub fn pop(&mut self) -> String {
        if self.uids.len() == 0 {
            self.fill();
        }
        *self.uids.pop()
    }
    fn fill(&mut self) {
        let mut n_fill = self.size - self.uids.len();
        while n_fill > 0 {
            // generate a new uid at random and add it to the stack
            self.uids.push(&self.gen());
        }
    }
    fn gen(&self) -> String {
        let mut rng = rand::thread_rng();
        let mut out = String::new();
        for _ in 0..self.id_len {
            let ix: usize = rng.gen_range(0, self.alphabet.len());
            out = out.add(&self.alphabet[ix]);
        }
        out
    }
}