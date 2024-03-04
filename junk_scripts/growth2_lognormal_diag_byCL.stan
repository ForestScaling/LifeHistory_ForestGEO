data {
int <lower = 1> n_stem;
int <lower = 1, upper = n_stem> n_lvl;
int <lower = 1, upper = n_lvl> levels[n_stem];
real <lower = 0> growth2[n_stem];

}

parameters {
real lognorm_mu_CL[n_lvl];
real <lower = 0> lognorm_sigma_CL[n_lvl];

}

model {
for(i in 1:n_stem){
target += lognormal_lpdf(growth2[i] | lognorm_mu_CL[levels[i]], lognorm_sigma_CL[levels[i]]);
}



}

generated quantities {
// define the posterior predictive "observations"
real <lower = 0> g2_rep[n_stem];

for(y in 1:n_stem){
g2_rep[y] = lognormal_rng(lognorm_mu_CL[levels[y]] , lognorm_sigma_CL[levels[y]]);

}

}
