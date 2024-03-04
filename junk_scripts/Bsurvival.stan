data {
int <lower = 1> n_stem;
int <lower = 1, upper = n_stem> n_spp;
int <lower = 1, upper = n_spp> species[n_stem];
int <lower = 0, upper = 1> surv_bindat[n_stem];
vector<lower = 0>[n_stem] t_interv;

}

transformed data{

vector[n_stem] t_5;
t_5 = t_interv/5;

}


parameters {
real <lower = 0, upper = 1> surv_p[n_spp];
real <lower = 0> alpha;
real <lower = 0> beta; 

}

model{
real p_rlz;

for(i in 1:n_stem){
p_rlz = pow(surv_p[species[i]], t_5[i]);
target += bernoulli_lpmf(surv_bindat[i] | p_rlz);
}

surv_p ~ beta(alpha, beta);
alpha ~ normal(0, 100);
beta ~ normal(0, 100);

}

