//
// This Stan program defines an epidemiological model.
// Recorded cumulative case counts on each day are considered to be
// normally distributed (with standard deviation sigma) realizations 
// of the sum of infected and recovered patients at the corresponding 
// timestep of a simple SIR model (approximated with Euler's method)
// with initial infected I_init, initial susceptible S_init, and 
// parameters beta and gamma.
//

// The input data are N, R_init, cases, and pop.
data {
  int<lower=0> N; // Number of observations
  real<lower=0> R_init; // Initial recovered
  real cases[N]; // case measurements
  real<lower=0> pop; // county population
}

// The parameters accepted by the model. Our model
// accepts the parameters I_init, S_init, gamma, beta, and sigma.
parameters {
  real<lower=0> I_init; // Initial infected
  real<lower=0> S_init; // Initial susceptible
  real<lower=0> gamma; // parameter of interest
  real<lower=0> beta; // parameter of interest
  real<lower=0> sigma; // Regression noise
}

transformed parameters {
      real<lower=0> I_vals[N]; // Theoretical infections
      real<lower=0> S_vals[N]; // Theoretical susceptible
      real<lower=0> R_vals[N]; // Theoretical recovered
      real<lower=0> case_vals[N]; // Theoretical infectious+recovered
      
      // Use the initial values of I, S, and R and Euler's method to estimate
      // the values of I, S, and R at successive timesteps. Note that more
      // complex and accurate numerical methods were considered but were too
      // time-consuming and did not notably improve performance.
      
      I_vals[1] = I_init;
      S_vals[1] = S_init;
      R_vals[1] = R_init;
      case_vals[1] = I_init + R_init;
      
      for(i in 2:N) {
        I_vals[i] = I_vals[i-1] + ((beta*S_vals[i-1]*I_vals[i-1]) - (gamma*I_vals[i-1]));
        S_vals[i] = S_vals[i-1] + (-1*beta*S_vals[i-1]*I_vals[i-1]);
        R_vals[i] = R_vals[i-1] + (gamma*I_vals[i-1]);
        case_vals[i] = I_vals[i] + R_vals[i];
      }
}

// The model to be estimated.
model {
    cases ~ normal(case_vals, sigma);
    S_init ~ normal(0.75*pop, sqrt(pop));
    I_init ~ lognormal(log(1.1*cases[1]), log(1.1*cases[1]));
    gamma ~ gamma(0.01, 1);
    beta ~ gamma((0.05/pop), 1);
    sigma ~ gamma(2,2);
}
