data{
  int n;
  vector[n] Day;
  vector[n] Proportion;
  array[n] int Species;
  int n_Species;
  array[n] int Sediment;
  int n_Sediment;
  vector[n] Load;
}

parameters{
  // Decay constant k
  real alpha_mu; // Global intercept
  real<lower=0> alpha_sigma_sp; // Interspecific variation
  real<lower=0> alpha_sigma_sed; // Sediment variation
  vector[n_Species] alpha_z_sp; // z-scores
  vector[n_Sediment] alpha_z_sed;
  real beta; // Slope over carbon load
  
  // "Refractory" proportion r
  real gamma_mu;
  real<lower=0> gamma_sigma_sp;
  real<lower=0> gamma_sigma_sed; 
  vector[n_Species] gamma_z_sp; 
  vector[n_Sediment] gamma_z_sed;
  real delta;
  
  // Likelihood precision
  real<lower=0> nu;
}

transformed parameters{
  // Convert z-scores
  vector[n_Species] alpha_sp = alpha_z_sp * alpha_sigma_sp + 0;
  vector[n_Sediment] alpha_sed = alpha_z_sed * alpha_sigma_sed + 0;
  vector[n_Species] gamma_sp = gamma_z_sp * gamma_sigma_sp + 0;
  vector[n_Sediment] gamma_sed = gamma_z_sed * gamma_sigma_sed + 0;
}

model{
  // Priors
  // Decay constant k
  alpha_mu ~ normal( log(0.1) , 1 );
  alpha_sigma_sp ~ normal( 0 , 0.3 )T[0,];
  alpha_sigma_sed ~ normal( 0 , 0.3 )T[0,];
  alpha_z_sp ~ normal( 0 , 1 );
  alpha_z_sed ~ normal( 0 , 1 );
  beta ~ normal( 0 , 0.1 );
  
  // "Refractory" proportion r
  gamma_mu ~ normal( logit(0.1) , 1 );
  gamma_sigma_sp ~ normal( 0 , 0.3 )T[0,];
  gamma_sigma_sed ~ normal( 0 , 0.3 )T[0,];
  gamma_z_sp ~ normal( 0 , 1 );
  gamma_z_sed ~ normal( 0 , 1 );
  delta ~ normal( 0 , 0.1 );
  
  // Likelihood precision
  nu ~ gamma( square(30) / square(20) , 30 / square(20) );
  
  // Model
  vector[n] k = exp(
    alpha_mu + alpha_sp[Species] + alpha_sed[Sediment] + beta * Load
  );
  vector[n] r = inv_logit(
    gamma_mu + gamma_sp[Species] + gamma_sed[Sediment] + delta * Load
  );
  vector[n] mu = (1 - r) .* exp( -k .* Day ) + r;

  // Beta likelihood
  Proportion ~ beta( mu * nu , (1 - mu) * nu );
}