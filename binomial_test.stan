data {
  int<lower=0> N;
  int n_seq[N];
  int delta_seq[N];
}s
parameters {
  real<lower = 0.0, upper = 1.0> theta[N];
}
model {
  theta ~ beta(1,1);
  delta_seq ~ binomial(n_seq, theta);
}

