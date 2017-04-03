//
//  Stan model for simple ODE model for cell counts
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
  // structural parameters
  // 
  // starting values
  row_vector<lower=0>[J] T1; // uninfected cells 
  row_vector<lower=0>[J] I1; // infected cells
  row_vector<lower=0>[J] V1; // viruses
  real <lower=0> lambda; // autonomous cell generation rate 
  real <lower=0> dT;     // death rate of uninfected cells
  real <lower=0> dI;     // death rate of infected cells
  real <lower=0> beta;   // infection rate 
  real <lower=0> p;      // production rate of viruses per infected cell
  real <lower=0> c;      // clearance rate of viruses
}
transformed parameters {
  //
  // Expected cell count in sample treated as
  // deterministic because it's happening
  // in a huge population and random variation
  // is small relative to sampling variation
  // 
  matrix[Nt,J] T;
  matrix[Nt,J] I;
  matrix[Nt,J] V;
  T[1,] = T1;
  I[1,] = I1;
  V[1,] = V1;
  for(t in 1:(Nt-1)) {
    row_vector[J] tcurr;
    row_vector[J] icurr;
    row_vector[J] vcurr;
    tcurr = T[t,];
    icurr = I[t,];
    vcurr = V[t,];
    T[t+1,] = tcurr + lambda - dT*tcurr - beta * tcurr .* vcurr;
    I[t+1,] = icurr + beta * tcurr .* vcurr - dI * icurr;
    V[t+1,] = vcurr + p*icurr - c*vcurr;
  }
}
model {
  // priors
  vector[N] Texp;
  vector[N] Iexp;
  vector[N] Vexp;
  beta ~ cauchy(0,2.5);
  lambda ~ cauchy(0,2.5);
  dT ~ cauchy(0,2.5);
  dI ~ cauchy(0,2.5);
  p ~ cauchy(0,2.5);
  c ~ cauchy(0,2.5);
  // likelihood
  for(n in 1:N){
    Texp[n] = T[time[n],id[n]];  
    Iexp[n] = I[time[n],id[n]];  
    Vexp[n] = V[time[n],id[n]];
   // if(n<10) print("n = ",n," T = ", Texp[n], " I = ",Iexp[n]," V = ",Vexp[n]);
   // if(n%100 == 0 ) print("n = ",n," T = ", Texp[n], " I = ",Iexp[n]," V = ",Vexp[n]);
  }
  Tc ~ normal(Texp, sd_factor * Texp);
  Ic ~ normal(Iexp, sd_factor * Iexp);
  Vc ~ normal(Vexp, sd_factor * Vexp);
}
