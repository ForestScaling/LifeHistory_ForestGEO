data {
int <lower = 1> n_stem;
int <lower = 1, upper = n_stem> n_spp;
int <lower = 1, upper = n_spp> species[n_stem];
real <lower = 0> growth2[n_stem];


}

parameters {
real sp_growth[n_spp];
real <lower = 0> sp_sigma[n_spp];
real mean_growth_mu;
real <lower = 0> mean_growth_sd;
real mean_sigma_mu;
real <lower = 0> mean_sigma_sd; 
}

// Likelihood
model {
for(i in 1:n_stem){
target += normal_lpdf(growth2[i]| sp_growth[species[i]], sp_sigma[species[i]]);
}

// partial pooling of species average growth rate
target += normal_lpdf(sp_growth | mean_growth_mu, mean_growth_sd);

// partial pooling of species variability in growth rate
target += lognormal_lpdf(sp_sigma | mean_sigma_mu, mean_sigma_sd);

// Priors on these are currently unspecified

}




