 $PARAM @annotated
  
  kd      : 0.418   : Death constant (1/time) // THETA(1)
  kg      : 0.838   : Growth constant (1/time) // THETA(2)
  a10     : 0       : Initial conditions // THETA(3) 
  a20     : 1       : Initial conditions // THETA(4) 
  
  x       :  1      : baseline covariate // boolean 0 or 1
  beta_x  :  0.1    : covariate effect baseline
  beta_kg  : 0.1    : covariate effect of kg
  beta_kd  : 0.1    : covariate effect of kg
  beta_base  : 0.1    : covariate effect of kg
  mu      :  0.1    : intercept
  

  shape   :  1     : shape parameter

$CMT
  
  A1 
  A2 // Biomarker
  
  p11 // survival probability

$MAIN
  
 A1_0 = a10;
 A2_0 = a20;
 p11_0 = 1.0;
 
 double eta_x = beta_x * x + beta_kg * log(kg) + beta_kd * log(kd) + beta_base * log(a20);
 
 double eta = exp(mu) * exp(eta_x);

$ODE
  
  dxdt_A1 = A2;
  dxdt_A2 = (kg - kd) * A2 + (kg*kd) * A1;
  

  
  if(SOLVERTIME >= 1E-5) {
    
    dxdt_p11 =  - p11 * shape * pow (SOLVERTIME, shape - 1) * eta;
  
  } else {
    
    dxdt_p11 = - p11 * eta; // approximate via exponential model
  
  }
  
  
$TABLE
    
    double CHANGE = (kg - kd) * A2 + (kg*kd) * A1;
    double CTDNA = A2*10;
  
$CAPTURE  
  CHANGE CTDNA