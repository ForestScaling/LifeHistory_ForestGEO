data {
// Number of stems, integer
int <lower = 1> n_stem;

// Integer of the number of species being estimated.
int <lower = 1, upper = n_stem> n_spp;

// Integers for assigning stems to the species, each stem has an integert identifying which species it belongs to
int <lower = 1, upper = n_spp> species[n_stem];

// Log growth values for each of the stems. Already calculated on an annual basis
real loggrowth_vals[n_stem];
}

parameters {
real z_sig[n_spp];

// mean growth for the species
real growparam_sp_CL[n_spp];


// mean and SD of average values for species growth 
real mu_growparam_CL;
real <lower = 0> sigma_growparam_CL;

// mean and SD of the SD of the growth values for the different species (lognormal)
real logmu_growsd_CL;
real <lower = 0> logsigma_growsd_CL;

}

transformed parameters{
real <lower = 0> growsd_sp_CL[n_spp];

for(j in 1:n_spp){
growsd_sp_CL[j] = exp(logmu_growsd_CL + logsigma_growsd_CL*z_sig[j]);
}
}




model {
// loops through the stems, likelihood of the observations for the species specific average and SD in growth



for(i in 1:n_stem){

target += normal_lpdf(loggrowth_vals[i] | growparam_sp_CL[species[i]], growsd_sp_CL[species[i]]);

}


// priors
// Heirachical drawing of mean values, drawn from a normal
growparam_sp_CL ~ normal(mu_growparam_CL, sigma_growparam_CL);

sigma_growparam_CL ~ exponential(0.1);
mu_growparam_CL ~ normal(0,10);

// Drawing values of SD heierarchically from a lognormal dist (forces SD to be >0)

z_sig ~ normal(0,1);

logsigma_growsd_CL ~ exponential(0.1);
logmu_growsd_CL ~ normal(0,10);

// priors on the hypoer parameters not specified (technically improper, but STAN doesn't 
// seem to mind.
// Constraints on variance to be >0 

}

generated quantities {
// define the posterior predictive "observations", simulate data from the posterior distribution.
// we should expect to see that the true data is qualitatively similar to the posterior
real g_rep[n_stem];
vector[n_stem] log_lik;

// run through all the stems, simulate a value of growth, we can then check things like Kurtosis, ske,
// etc. to make sure that our data is similar in shape to that which is implied by the model.
// loop through the stems

for(y in 1:n_stem){

//simulated  value
g_rep[y] = normal_rng(growparam_sp_CL[species[y]] , growsd_sp_CL[species[y]]);

// log likelihood
log_lik[y] = normal_lpdf(loggrowth_vals[y] | growparam_sp_CL[species[y]], growsd_sp_CL[species[y]]);
}

// End of loop

}













