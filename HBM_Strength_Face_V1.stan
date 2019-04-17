// Fitting Single subjects' 'recall errors across 

data {
int<lower=1> ncond;    // number of trials
int<lower=1> ntrials;    // number of trials per condition per sub
int<lower=1> nsub;      // number of subjects for each conditions
real<lower= -pi(),upper = pi()>  allerrors [ntrials];    //Errors
vector [ncond] RadiusLevel; // all z hand
}

parameters {
real <lower=0.001> Scale_kappa1;
real <lower=0.001> Shape_kappa1;
real <lower=0.001, upper = 100> kappa1 [nsub];

real <lower=0.001, upper = 30> alpha_pm;
real <lower=0.001, upper = 30> beta_pm;
real <lower=0.001, upper = 1> pm [nsub];

real <lower = -pi()/6, upper = pi()/6> GrpMu1;
real <lower = -pi()/6, upper = pi()/6> GrpMu2;

real <lower = -pi()/6, upper = pi()/6>  Prior_GrpMu1;
real <lower = -pi()/6, upper = pi()/6>  Prior_GrpMu2;

real <lower=0.0001> SD_mu1;
real <lower = -pi()/6, upper = pi()/6> mu1 [nsub];
real <lower = mu1 [nsub], upper = pi()/6> mu2 [nsub];
}


transformed parameters {

}


model {
    Scale_kappa1 ~ gamma(2,2);
    Shape_kappa1 ~ gamma(2,2);
    alpha_pm ~ gamma(1,1);
    beta_pm ~ gamma(1,1);
    SD_mu1 ~ gamma(0.5,0.5); 

    for(isub in 1:nsub){
    
    target += normal_lpdf(mu1[isub] | GrpMu1, SD_mu1);
    target += normal_lpdf(mu2[isub] | GrpMu2, SD_mu1);

    pm[isub] ~ beta(alpha_pm, beta_pm);
    kappa1[isub] ~ gamma(Shape_kappa1, Scale_kappa1);    

    for(i in 1:trials){         
           target += log_sum_exp(log(pm[isub])+log_sum_exp(log(0.5)+von_mises_log(errors[i,isub],mu1[isub],kappa1[isub]), log(0.5)+von_mises_log(errors[i,isub],mu2[isub],kappa1[isub])), log((1-pm[isub]) / (2*pi())));

}
}
}

generated quantities {
real GrpMuDiff;
real Prior_GrpMuDiff;

GrpMuDiff = GrpMu2 - GrpMu1;    
Prior_GrpMuDiff = Prior_GrpMu2 - Prior_GrpMu1;
}
