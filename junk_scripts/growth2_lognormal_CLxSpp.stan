data {
int <lower = 1> n_stem;
int <lower = 1, upper = n_stem> n_spp;
int <lower = 1, upper = n_spp> species[n_stem];
real <lower = 0> growth2[n_stem];


}

parameters {
real lognorm_mu_CL[n_spp];
real <lower = 0> lognorm_sigma_CL[n_spp];

}

model {
for(i in 1:n_stem){
target += lognormal_lpdf(growth2[i] | lognorm_mu_CL[species[i]], lognorm_sigma_CL[species[i]]);
}



}

generated quantities {
// define the posterior predictive "observations"
real <lower = 0> g2_rep[n_stem];

for(y in 1:n_stem){
g2_rep[y] = lognormal_rng(lognorm_mu_CL[species[y]] , lognorm_sigma_CL[species[y]]);

}

}
