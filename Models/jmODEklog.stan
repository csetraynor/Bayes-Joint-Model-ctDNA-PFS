// Joint model

functions{
  
  vector strictly_positive(vector x) {
    int N = rows(x);
    vector[N] res;
    for (n in 1:N)
      res[n] = fmax(x[n], machine_precision());
    return res;
  }
  
  // define ODE system for Joint model
  vector ModelJMODE(real t, vector x, real[] parms, real[] rdummy, int[] idummy) {
      
    // Parameters

    real kg       = parms[1];
    real kd       = parms[2];
    real baseline = parms[3];
    
    real eta   = parms[4];
    real shape   = parms[5];
    
    // Return object (derivative)
    vector[3] y; // 3 element per compartment of
    
    // the model
    real A1 = x[1];
    real A2 = fmax(machine_precision(), x[2] + baseline); // model initial conditions
    real p11 = fmax(machine_precision(), x[3] + 1.0);  // Survival Probability
    
    // PD component of the ODE system
    y[1] = A2;
    y[2] = (kg - kd) * A2 + (kg*kd) * A1;
    // survival equation (Weibull)
    y[3] = - p11 * shape * pow (fmax(machine_precision(), t) , shape - 1.0) * eta;
    
    return y;
  }
}

data  {
  
  int<lower = 1> nId;
  int<lower = 1> nt;
  int<lower = 1> nObsPD; // number of observations
  int<lower = 1> iObsPD[nObsPD]; // index of observation
  int<lower = 1> iObsSURV[nId]; // index of survival observation
  
  // NONMEM data
  int<lower = 1> cmt[nt];
  int<lower = 1> start[nId];
  int<lower = 1> end[nId];
  
  int evid[nt];
  real amt[nt];
  real time[nt];
  real<lower = 0> rate[nt];
  real<lower = 0> ii[nt];
  int<lower = 0> addl[nt];
  
  vector<lower = 0>[nObsPD] PDObs; // observed ADC concentration (dependent variable)
  
  // survival information
  vector<lower = 0, upper = 1>[nId] censoring_status;
  vector<lower = 0>[nId] survival_time; 
  
  vector<lower = 0, upper = 1>[nId] treatment; // treatment indicator
  
  // BLQ
  vector<lower = 0>[nObsPD] ISBLQ_PD; 
  
  
  // priors
  real KGPrior ;
  real KDPrior ;
  real BASEPrior ;
  
  real beta_g_Prior;
  real beta_x_Prior;
  real mu0Prior;
  real shapePrior;
  
  real priorSigmaPD;
  
  real KGPriorCV ;
  real KDPriorCV ;
  real BASEPriorCV ;
  
  real beta_g_PriorCV;
  real beta_x_PriorCV;
  real mu0PriorCV;
  real shapePriorCV;
  
}
  
transformed data {
  int<lower = 0> ss[nt] = rep_array(0, nt);
  int<lower = 1> nCmt = 3;
  real biovar[nt, nCmt];
  real tlag[nt, nCmt];
  vector[nObsPD] logPDObs = log(PDObs);
  
  int<lower = 1> nParam = 5;
  int<lower = 1> nRandom = 3;
  int len[nId];
  
  int nPred = 1;

  for (j in 1:nId) {
    len[j] = end[j] - start[j] + 1;
  }
  
  for (j in 1:nt) {
    for (i in 1:nCmt) {
      biovar[j, i] = 1;
      tlag[ j, i] = 0;
      }
  }
  
}

parameters  {
  
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
  
  // covariate effects
  
  matrix[nRandom, nId] etaStd;
  cholesky_factor_corr[nRandom] L;
  vector<lower = 0, upper = 1>[nRandom] omega;
}

transformed parameters{
  vector<lower = 0>[nRandom] thetahat =
    to_vector({KDhat, KGhat, BASEhat});
  matrix<lower = 0>[nId, nRandom] theta;
  theta = (rep_matrix(thetahat, nId) .* 
	   exp(diag_pre_multiply(omega, L * etaStd)))';

}

model{
  vector[nt] PD;
  vector[nt] SURV;
  
  vector[nObsPD] PDhatObs; // predicted PD concentration 
  vector[nId] SurvhatObs; // predicted PD concentration
  vector[nId] PDatEventHatObs;
  
  matrix[nCmt, nt] x;
  real parms[nId, nParam];
  real KG[nId];
  real KD[nId];
  real BASE[nId];
  real eta[nId];
  
  // PK model is allometrically scaled using West et al. assumptions https://www.science.org/doi/10.1126/science.276.5309.122
  for(j in 1:nId){ 
    
    KD[j]    = theta[j, 1] ;
    KG[j]    = theta[j, 2] ;
    BASE[j]  = theta[j, 3] ;
    
    eta[j]   = mu0 + beta_x * treatment[j] + beta_kg * log(KG[j]) + beta_kd * log(KD[j]) + beta_base * log(BASE[j]);
    
    parms[j, ] = {KG[j], KD[j], BASE[j] , exp(eta[j]), shape};
  }
  
  // print("params : ", parms);
  
  x = pmx_solve_group_rk45(ModelJMODE, nCmt, len,
                         time, amt, rate, ii, evid, cmt, addl, ss,
                         parms,
                         1e-6, 1e-6, 1E5);

  for(j in 1:nId) {
      PD[start[j]:end[j]] = x[2, start[j]:end[j]]' + BASE[j] ;
      SURV[start[j]:end[j]] = x[3, start[j]:end[j]]' + 1.0 ;
    }
    
  PDhatObs = strictly_positive(PD[iObsPD]);
  SurvhatObs = strictly_positive(SURV[iObsSURV]);
  PDatEventHatObs = strictly_positive(PD[iObsSURV]);
  
  // informative prior
  KDhat ~ lognormal(log(KDPrior), KDPriorCV);
  KGhat ~ lognormal(log(KGPrior), KGPriorCV);
  BASEhat ~ lognormal(log(BASEPrior), BASEPriorCV);
  
  beta_kg  ~ std_normal();
  beta_kd  ~ std_normal();
  beta_base  ~ std_normal();
  beta_x ~ std_normal();
  
  mu0 ~ normal(mu0Prior, mu0PriorCV);
  shape ~ normal(shapePrior, shapePriorCV);
  
  // SIGMA
  sigmaPD ~ exponential(priorSigmaPD);
  
  // OMEGA
  omega ~ normal(0, 0.5);
  L ~ lkj_corr_cholesky(1);
  
  // Inter-individual variability
  to_vector(etaStd) ~ normal(0, 1);
  
  // observed data likelihood
  target += normal_lpdf(logPDObs | log(PDhatObs), sigmaPD); 
  
  // survival likelihood
  for(j in 1:nId) {
    if(censoring_status[j] == 1){
      target += log(SurvhatObs[j]);
      target += log(shape) + (shape - 1) * log(survival_time[j]) + eta[j];
    } else {
      target += log(SurvhatObs[j]);
    }
  }
}

generated quantities {
  
  vector[nId] KG_IPRED;
  vector[nId] KD_IPRED;
  vector[nId] BASE_IPRED;
  
  // Variables for IIV  
  matrix[nRandom, nPred] etaStdPred;
  matrix<lower = 0>[nPred, nRandom] thetaPredM;
  corr_matrix[nRandom] rho;

  rho = L * L';
  
  for(i in 1:nPred) {
    for(j in 1:nRandom) {
      etaStdPred[j, i] = normal_rng(0, 1);
    }
  }
    
  thetaPredM = (rep_matrix(thetahat, 1) .* 
                exp(diag_pre_multiply(omega, L * etaStdPred)))';
                
   for (j in 1:nId) {
     
    KD_IPRED[j]    = theta[j, 1] ;
    KG_IPRED[j]    = theta[j, 2] ;
    BASE_IPRED[j]  = theta[j, 3] ;
    
  }
  
}
