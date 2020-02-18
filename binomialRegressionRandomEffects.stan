data {
    int<lower=1> Ntotal; // number of observations
    int<lower=1> Nlevels1; // number of random intercepts for factor 1 = nlevels of factor 1 
    int<lower=1> NfactorMap1[Ntotal]; // mapping variable to map each observation to a particular level of factor 1
    int y[Ntotal]; // response variable binomial distributed
    int Ntrials[Ntotal]; // total trials for binomial variable
}
parameters {
  // parameters to estimate in the model
    real intercept; // constant term
    //real<lower=0.01> sigmaFactor1; // population level sd for factor 1
    vector[Nlevels1] nCoefFactor1; // number of random jitters, for levels of factor 1
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  mu = nCoefFactor1[NfactorMap1];
  mu = inv_logit(mu);
}
model {
  // using generic priors to start with
  //sigmaFactor1 ~ exponential(1);
  intercept ~ normal(0, 2);
  nCoefFactor1 ~ normal(intercept, 1);
  // likelihood function
  y ~ binomial(Ntrials, mu);
}
