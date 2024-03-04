data {
int <lower = 1> n_stem;
int <lower = 1, upper = n_stem> n_spp;
int <lower = 1, upper = n_spp> species[n_stem];
real growth_vals[n_stem];


}

parameters {

real growparam_sp_CL[n_spp];
real <lower = 0> growsd_sp_CL[n_spp];
real mu_growparam_CL;
real <lower = 0> sigma_growparam_CL;
real logmu_growsd_CL;
real <lower = 0> logsigma_growsd_CL;
}

model {
for(i in 1:n_stem){
target += lognormal_lpdf(growth_vals[i] | growparam_sp_CL[species[i]], growsd_sp_CL[species[i]]);
}

//priors
target += normal_lpdf(growparam_sp_CL | mu_growparam_CL, sigma_growparam_CL);
target += lognormal_lpdf(growsd_sp_CL | logmu_growsd_CL, logsigma_growsd_CL);

}

generated quantities {
// define the posterior predictive "observations"
real g_rep[n_stem];

for(y in 1:n_stem){
g_rep[y] = lognormal_rng(growparam_sp_CL[species[y]] , growsd_sp_CL[species[y]]);

}

}
