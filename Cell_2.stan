//
//  Stan model for simple ODE model for cell concentrations
//

data {
  int N; // total observations
  int J; // number of mice
  int Nt;  // length of simulated times
  int<lower=1,upper=Nt> time[N];
  int<lower=1,upper=J> id[N];
  vector[N] Tc;  // uninfected cell concentration 
  vector[N] Ic;  // infected cell concentration
  vector[N] Vc;  // viral concentration
  real s;        // volume of sample
  real C;        // reference volume for concentration
}
transformed data {
  real sd_factor;
  sd_factor = sqrt(C/s);
}
parameters {

  // starting values

  vector<lower=0>[J] T1; // uninfected cells 
  vector<lower=0>[J] I1; // infected cells
  vector<lower=0>[J] V1; // viruses

  // parameters

  real <lower=0> lambda; // autonomous cell generation rate 
  real <lower=0> dT;     // death rate of uninfected cells
  real <lower=0> dI;     // death rate of infected cells
  real <lower=0> beta;   // infection rate 
  real <lower=0> p;      // production rate of viruses per infected cell
  real <lower=0> c;      // clearance rate of viruses
}
model {
  int t;
  vector[N] Texp;  // expected values
  vector[N] Iexp;
  vector[N] Vexp;
  vector[J] tcurr; // incremental values 
  vector[J] icurr;
  vector[J] vcurr;
  vector[J] tnext; 
  vector[J] inext;
  vector[J] vnext;
  
  // priors 
  
  beta ~ normal(0,.1);
  lambda ~ normal(0,.1);
  dT ~ normal(0,.1);
  dI ~ normal(0,.1);
  p ~ normal(0,.1);
  c ~ normal(0,.1);
  
  // model
  
  tcurr = T1;
  icurr = I1;
  vcurr = V1;
  t = 1;
  for(n in 1:N) {
    while( t < time[n]){
      tnext = tcurr + lambda - dT*tcurr - beta * tcurr .* vcurr;
      inext = icurr + beta * tcurr .* vcurr - dI * icurr;
      vnext = vcurr + p*icurr - c*vcurr;
      tcurr = tnext;
      icurr = inext;
      vcurr = vnext;
      t = t + 1;
    }
    Texp[n] = tcurr[id[n]];
    Iexp[n] = icurr[id[n]];
    Vexp[n] = vcurr[id[n]];
  }

  // likelihood

  Tc ~ normal(Texp, sd_factor * Texp);
  Ic ~ normal(Iexp, sd_factor * Iexp);
  Vc ~ normal(Vexp, sd_factor * Vexp);
}
