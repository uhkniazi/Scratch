data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Nclusters1; // number of levels for group 1 for random intercepts
  int<lower=1> Nclusters2; // number of levels for group 2 for random intercepts
  int<lower=1> Nclusters3; // number of levels for group 3 for random intercepts
  int<lower=1> NScaleBatches1; // number of batches of scale terms 
  int<lower=1, upper=Nclusters1> NgroupMap1[Ntotal]; // mapping variable to map each observation to group 1 
  int<lower=1, upper=Nclusters2> NgroupMap2[Ntotal]; // mapping variable to map each observation to group 2 
  int<lower=1, upper=Nclusters3> NgroupMap3[Ntotal]; // mapping variable to map each observation to group 3 
  int<lower=1, upper=NScaleBatches1> NBatchMap1[Nclusters1]; // expanding vector to map each variance term to the relevant 
                                                    // set of coefficients from rGroupsJitter1
  int y[Ntotal]; // response variable
  // additional parameters
  real gammaShape; // hyperparameters for the gamma distribution 
  real gammaRate;
}
// transformed data {
  // }
parameters {
  // parameters to estimate in the model
  real betas; // regression parameters
  real<lower=0.01> sigmaRan1[NScaleBatches1]; // random effect standard deviations for sub-batches in group 1
  real<lower=0.01> sigmaRan2; // random effect standard deviation for group 2
  real<lower=0.01> sigmaRan3; // random effect standard deviation for group 3
  vector[Nclusters1] rGroupsJitter1; // number of random jitters for each level of cluster/group 1
  vector[Nclusters2] rGroupsJitter2; // number of random jitters for each level of cluster/group 2
  vector[Nclusters3] rGroupsJitter3; // number of random jitters for each level of cluster/group 3
  real<lower=1, upper=100> iSize; // size parameter for the nb distribution
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  mu = exp(betas + rGroupsJitter1[NgroupMap1] + rGroupsJitter2[NgroupMap2]*sigmaRan2 + rGroupsJitter3[NgroupMap3]*sigmaRan3);
}
model {
  real sigmaRan1_expanded[Nclusters1]; 
  sigmaRan1 ~ gamma(gammaShape, gammaRate);
  sigmaRan2 ~ gamma(gammaShape, gammaRate);
  sigmaRan3 ~ gamma(gammaShape, gammaRate);
  // random effects sample
  sigmaRan1_expanded = sigmaRan1[NBatchMap1];
  rGroupsJitter1 ~ normal(0, sigmaRan1_expanded);
  rGroupsJitter2 ~ normal(0, 1);
  rGroupsJitter3 ~ normal(0, 1);
  // likelihood function
  y ~ neg_binomial_2(mu, iSize);
}
