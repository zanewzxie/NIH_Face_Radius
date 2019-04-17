
data {
int<lower=1> nhand;    // number of handcondition
int<lower=1> nsub;      // number of subjects for each conditions
int<lower=1> ncolor;    // number of set size 
int AllSub [ncolor*nhand*nsub]; // all sub
int AllHand [ncolor*nhand*nsub]; // all hand
int AllColor [ncolor*nhand*nsub]; // all hand
vector [ncolor*nhand*nsub] y; // all y
vector [nhand] handlevel; // all z hand
vector [ncolor] colorlevel; // all z color
}


parameters {
real <lower = 0.0001>  Sigma ;
vector <lower = 0.0001> [nsub] gamma;
real <lower = 0.0001>  gamma_mu;  
real <lower = 0.0001>  gamma_sd; 

real <lower = 0.0001> fa_mu;                     
real <lower = 0.0001> fb_mu;    
real <lower = 0.0001> sa_mu;
real <lower = 0.0001> sb_mu; 

vector <lower = 0.0001> [nsub] fa;
vector <lower = 0.0001> [nsub] fb;
vector <lower = 0.0001> [nsub] sa;
vector <lower = 0.0001> [nsub] sb;

real <lower = 0.0001> fa_sd;                     
real <lower = 0.0001> fb_sd;    
real <lower = 0.0001> sa_sd;
real <lower = 0.0001> sb_sd; 
}

transformed parameters {
matrix [nhand,nsub] UF ;   // Utility associated with each hand level for each subject
matrix [ncolor,nsub] US ;   // Utility associated with each color level for each subject
matrix [nhand,nsub] EF ;   // Effort associated with each hand level for each subject
matrix [ncolor,nsub] ES ;   // Effort associated with each color level for each subject
vector  [ncolor*nhand*nsub] Predict; // predict probability

for (i in 1:ncolor*nhand*nsub) {
EF[AllHand[i],AllSub[i]] = (fb[AllSub[i]]+pow(handlevel[AllHand[i]],fa[AllSub[i]]));
UF[AllHand[i],AllSub[i]] =  -log(EF[AllHand[i],AllSub[i]]); // for 5 hand levels

ES[AllColor[i],AllSub[i]] = (sb[AllSub[i]]+pow(colorlevel[AllColor[i]],sa[AllSub[i]])); 
US[AllColor[i],AllSub[i]] =  -log(ES[AllColor[i],AllSub[i]]); // for 5 color levels

Predict[i] = Phi((UF[AllHand[i],AllSub[i]] - US[AllColor[i],AllSub[i]])/gamma[AllSub[i]]);
}
}

model {
    fa_mu ~ normal(1,10)T[0.0001,];
    fb_mu ~ normal(0,100)T[0.0001,];
    sa_mu ~ normal(1,10)T[0.0001,];
    sb_mu ~ normal(0,100)T[0.0001,];

    fa_sd ~ cauchy(0,20);
    fb_sd ~ cauchy(0,20);
    sa_sd ~ cauchy(0,20);
    sb_sd ~ cauchy(0,20);
    gamma_mu ~ cauchy(0,20);
    gamma_sd ~ cauchy(0,20);
    Sigma ~  cauchy(0,20);

for (isub in 1:nsub) {
       fa[isub] ~ normal(fa_mu,fa_sd)T[0.0001,];
       sa[isub] ~ normal(sa_mu,sa_sd)T[0.0001,];
       fb[isub] ~ normal(fb_mu,fb_sd)T[0.0001,];
       sb[isub] ~ normal(sb_mu,sb_sd)T[0.0001,];

       gamma[isub] ~ normal(gamma_mu,gamma_sd)T[0.0001,];
}

    for (i in 1:ncolor*nhand*nsub) {
     y[i] ~ normal(Predict[i],Sigma);
}
}


generated quantities {
   vector[ncolor*nhand*nsub] log_lik;

   for (n in 1:ncolor*nhand*nsub){
      log_lik[n] <- normal_log(y[n], Predict[n], Sigma);
   }

}

