data {
int <lower = 1> n_obs;
int <lower = 1,  upper = n_obs> n_OTU;
int <lower = 1, upper = n_OTU> OTU[n_obs];
real <lower = 0, upper = 1> propn[n_obs];
}

parameters {
// individual species survival 
real <lower = 0> alpha_sp_surv_CL[n_OTU];
real <lower = 0> beta_sp_surv_CL[n_OTU];

// Priors for the parameters of the beta distribution
real <lower = 0> alpha_avg_CL;
real <lower = 0> alpha_sdv_CL;

real <lower = 0> beta_avg_CL;
real <lower = 0> beta_sdv_CL;
}

model {
for(i in 1:n_obs){
target += beta_lpdf(propn[i] | alpha_sp_surv_CL[OTU[i]], beta_sp_surv_CL[OTU[i]]);
}

/// Half normal since alpha nad beta are constrained to be >0
target += lognormal_lpdf(alpha_sp_surv_CL | alpha_avg_CL, alpha_sdv_CL);
target += lognormal_lpdf(beta_sp_surv_CL | beta_avg_CL, beta_sdv_CL);



}


