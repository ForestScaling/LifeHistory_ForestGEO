data {
//number of stems in the data must be integer >=1
int <lower = 1> n_stem; 

// number of species in the dataset, integer 
int <lower = 1, upper = n_stem> n_spp; 

// vector of integers used to denote species
int <lower = 1, upper = n_spp> species[n_stem]; 

// Binary survival 1-true, 0-false
int <lower = 0, upper = 1> bindat[n_stem]; 

// time interval between observations
vector<lower = 0>[n_stem] t_interv; 

}

transformed data{
// Force time to be on 5yr, value is arbitrary (makes the math easier for the computer)
vector[n_stem] t_5; 
t_5 = t_interv/5; 

}

parameters{
// sequence of odds ratios of survival for each species
real odds_rat[n_spp]; 

// Average odds ratio pooling across all of the species
real mean_odds; 

// SD of odds ratio 
real <lower = 0> sd_odds; 
}


transformed parameters{
// Defines the probabilities of survival
real <lower = 0, upper = 1> p[n_spp];

// Transforms survival into a binomial probability, I don't know of a clean closed form for accounting
// for different durations between observations without modifying p directly.
p = inv_logit(odds_rat);

}

model{
// define probability of survival, accountign for difference in length between obs. E.g. if obs
// for species A are at different time intervals we need to make sure that the "true" probability
// reflects the difference
real p_rlz;

// loops through all the stems
for(i in 1:n_stem){

// stem specific "fix" for different time duration between observations
p_rlz = pow(p[species[i]], t_5[i]);

// update posterior
target += bernoulli_lpmf(bindat[i] | p_rlz);

}

//Close the loop

// prior on mean_odds is "non-informative uniform" on -Inf +Inf, 
// prior on sd_odds is "uniform" >0
// I think the machine artificially keeps these from being improper, if necessary we can always change these to some 
// really diffuse prior
odds_rat ~ normal(mean_odds, sd_odds);

}




