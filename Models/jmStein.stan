// Joint model — analytical (closed-form) version.
// Stein-Fojo biomarker: f(t) = f_0 * (exp(-k_d t) + exp(k_g t) - 1)
// Weibull survival:    S(t) = exp(-t^gamma * exp(eta))
//
// Stan >= 2.33. Same standata as jmODEklog.stan; unused NONMEM fields
// are simply ignored.

data {

  int<lower = 1> nId;
  int<lower = 1> nObsPD;

  // Time and index of ctDNA observations
  vector<lower = 0>[nObsPD] time_ctdna;      // observation times, one per obs
  array[nObsPD] int<lower = 1, upper = nId> subject_of_obs;

  // Observed biomarker
  vector<lower = 0>[nObsPD] PDObs;

  // Survival information
  vector<lower = 0, upper = 1>[nId] censoring_status;
  vector<lower = 0>[nId] survival_time;
  vector<lower = 0, upper = 1>[nId] treatment;

  // BLQ indicator (unused in current likelihood; kept for consistency with
  // ODE model standata — remove if you want a cleaner interface)
  vector<lower = 0>[nObsPD] ISBLQ_PD;

  // Priors
  real KGPrior;
  real KDPrior;
  real BASEPrior;
  real mu0Prior;
  real shapePrior;
  real priorSigmaPD;

  real KGPriorCV;
  real KDPriorCV;
  real BASEPriorCV;
  real mu0PriorCV;
  real shapePriorCV;
}

transformed data {
  int<lower = 1> nRandom = 3;
  vector[nObsPD] logPDObs = log(PDObs);
}

parameters {
  real<lower = 0> KDhat;
  real<lower = 0> KGhat;
  real<lower = 0> BASEhat;

  real beta_kg;
  real beta_kd;
  real beta_base;
  real beta_x;
  real mu0;
  real<lower = 0> shape;

  real<lower = 0> sigmaPD;

  matrix[nRandom, nId] etaStd;
  cholesky_factor_corr[nRandom] L;
  vector<lower = 0>[nRandom] omega;
}

transformed parameters {
  vector<lower = 0>[nRandom] thetahat = to_vector({KDhat, KGhat, BASEhat});
  matrix<lower = 0>[nId, nRandom] theta;
  vector[nId] eta;

  theta = (rep_matrix(thetahat, nId) .*
           exp(diag_pre_multiply(omega, L * etaStd)))';

  for (j in 1:nId) {
    eta[j] = mu0
           + beta_x    * treatment[j]
           + beta_kg   * log(theta[j, 2])   // KG
           + beta_kd   * log(theta[j, 1])   // KD
           + beta_base * log(theta[j, 3]);  // BASE
  }
}

model {
  // ── Priors ───────────────────────────────────────────────────
  KDhat   ~ lognormal(log(KDPrior),   KDPriorCV);   // KDPriorCV bumped to 0.3
  KGhat   ~ lognormal(log(KGPrior),   KGPriorCV);   // KGPriorCV bumped to 0.3
  BASEhat ~ lognormal(log(BASEPrior), BASEPriorCV); // BASEPriorCV bumped to 0.3
  
  beta_kg   ~ std_normal();
  beta_kd   ~ std_normal();
  beta_base ~ std_normal();
  beta_x    ~ std_normal();
  
  mu0   ~ normal(mu0Prior, mu0PriorCV);
  shape ~ lognormal(log(shapePrior), shapePriorCV); // was normal, now lognormal
  
  sigmaPD ~ normal(0, 0.5);   // was exponential, now half-normal
  
  omega ~ normal(0, 1);     // was normal(0, 0.5), truncated to (0,1); now half-normal(0, ∞)
  L     ~ lkj_corr_cholesky(1);
  
  to_vector(etaStd) ~ normal(0, 1);

  // ── ctDNA likelihood (analytical Stein-Fojo) ─────────────────
  {
    vector[nObsPD] pd_pred;
    for (n in 1:nObsPD) {
      int j = subject_of_obs[n];
      real kd_j   = theta[j, 1];
      real kg_j   = theta[j, 2];
      real base_j = theta[j, 3];
      real t      = time_ctdna[n];
      pd_pred[n] = fmax(machine_precision(),
                        base_j * (exp(-kd_j * t) + exp(kg_j * t) - 1.0));
    }
    target += normal_lpdf(logPDObs | log(pd_pred), sigmaPD);
  }

  // ── Survival likelihood (analytical Weibull) ─────────────────
  for (j in 1:nId) {
    real log_S = -pow(survival_time[j], shape) * exp(eta[j]);
    target += log_S;
    if (censoring_status[j] == 1) {
      // Event observed: add log-hazard
      target += log(shape) + (shape - 1) * log(survival_time[j]) + eta[j];
    }
  }
}

generated quantities {

  // Log-likelihood: per observation
  // First nObsPD entries: one per ctDNA observation.
  // Next  nId    entries: one per subject's survival record.
  vector[nObsPD + nId] log_lik;

  // Posterior predictive
  vector[nObsPD] logPD_rep;
  vector<lower = 0>[nId] event_time_rep;

  corr_matrix[nRandom] rho = multiply_lower_tri_self_transpose(L);

  {
    // ctDNA log-lik and PPC replicates
    for (n in 1:nObsPD) {
      int j = subject_of_obs[n];
      real kd_j   = theta[j, 1];
      real kg_j   = theta[j, 2];
      real base_j = theta[j, 3];
      real t      = time_ctdna[n];
      real pd_pred = fmax(machine_precision(),
                          base_j * (exp(-kd_j * t) + exp(kg_j * t) - 1.0));
      real mu_n = log(pd_pred);
      log_lik[n]   = normal_lpdf(logPDObs[n] | mu_n, sigmaPD);
      logPD_rep[n] = normal_rng(mu_n, sigmaPD);
    }

    // Survival log-lik
    for (j in 1:nId) {
      real ll = -pow(survival_time[j], shape) * exp(eta[j]);
      if (censoring_status[j] == 1) {
        ll += log(shape) + (shape - 1) * log(survival_time[j]) + eta[j];
      }
      log_lik[nObsPD + j] = ll;
    }

    // Posterior predictive event times
    // T ~ Weibull(shape, scale) with scale = exp(-eta / shape)
    for (j in 1:nId) {
      real scale_j = exp(-eta[j] / shape);
      event_time_rep[j] = weibull_rng(shape, scale_j);
    }
  }
}