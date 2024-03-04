data {
// Notice how we don't have real growth values in here;
int <lower = 1> n_stem;
int <lower = 1, upper = n_stem> n_spp;
int <lower = 1, upper = n_spp> species[n_stem];
vector<lower = 0>[n_stem] t_interv;

}

transformed data{
vector[n_stem] t_5;
t_5 = t_interv/5;
}

generated quantities{
// parameters from the underlying beta distr.
real <lower = 0> alpha;
real <lower = 0> beta;
// survival simulated from prior only
vector[n_spp] surv_p;
// Randomly draw survival probabilities for all the spp

//beta and alpha draws
alpha = abs(normal_rng(0, 25000));
beta =  abs(normal_rng(0, 25000));



for(j in 1:n_spp){
surv_p[j] = beta_rng(alpha, beta);
}

vector[n_stem] surv_bindat_pred;
// realized survival
real p_rlz;
for(i in 1:n_stem){
p_rlz = pow(surv_p[species[i]], t_5[i]);
surv_bindat_pred[i] = bernoulli_rng(p_rlz);
}
}
