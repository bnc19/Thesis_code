data{
   int<lower = 1> C; // No. current immune status
   int<lower = 1> B; // No. baseline serostatus
   int<lower = 1> K; // No. serotypes
   int<lower = 1> D; // No. vcd data
   int<lower = 1> V; // No. trial arms 
   int<lower = 1> J; // No. age groups
   int<lower = 1> T; // No. model time points (0:53)
   int<lower = 1> R; // No. outcomes 
   int<lower = 0> time[T];  // time 
   real<lower = 0> mu[B,K]; // titres at T0
   int<lower = 0> HI;       // period of heterotypic immunity
   array[J] int<lower=0> pop_J;  // baseline pop
   array[B,V,J,D] int<lower=0> pop;   // pop
   // FLAGS
   int<lower = 0, upper = 1> include_eps;      // include enhanced secondary hosp? (T/F)
   int<lower = 0, upper = 3> include_beta;     // include age-specific nc50?: 1= change age groups 1 & 2 for both outcomes / 2 = only age grp 1, sep for each outcome / 3 = only age grp 1 for both outcomes 
 int<lower = 0, upper = 2> mono_lc_SN;       // 0 for serotype specific lc SN / 1 for single / 2 for offset from MO
   int<lower = 0, upper = 2> mono_lc_MU;       // 0 for serotype specific lc MO / 1 for single / 2 for offset from MO
   int<lower = 0, upper = 1> rho_K;            // 0 for mono rho, 1 for serotype rho 
   int<lower = 0, upper = 1> L_K;              // 0 for mono L, 1 for serotype K 
   int<lower = 0, upper = 2> w_CK;             // 0 for mono w, 1 for serostatus w, 2 for serotype w 
   int<lower = 0, upper = 2> alpha_CK;         // 0 for mono alpha, 1 for serostatus alpha, 2 for serotype alpha 
   int<lower = 0, upper = 1> tau_K;            // 0 for mono tau, 1 for serotype tau
   int<lower = 0, upper = 1> MU_test_SN;       // multitypic test SN at baseline (T/F)
   int<lower = 0, upper = 1> MU_symp;          // 0 for only 1' and 2' symp infections, 1 for post-sec symp inf   
   int<lower = 0, upper = 1> enhancement;      // 0 for no vac enhancement of SN, 1 for enhancement 
   int<lower = -1, upper=0> lower_bound_L;     // enhancement lower bound 
  // PARAMETERS 
  real<lower = 0> hs[B];                       // short term decay 
  real<lower = 12> hl;                         // long term decay 
  real<lower = -12> ts[B];                     // time switch decay rate
  real<lower = 0> lambda_D[K,D];               // FOI  
  real<lower = 0, upper = 1> p[J];             // serotype-specific prob exposure  
  real<lower = 1> epsilon ;                    // prob hosp 2' rel to 1'/3'/4'
  real<lower = 0, upper = 1> gamma;            // prob symp 2' 
  real<lower = 0> rho[K];                      // 1/risk symp 1' rel to 2'
  real<lower = 0, upper = 1> phi;              // prob symp 3'/4' rel to 1'
  real<lower = 0, upper = 1> delta[K];         // prob hosp
  real<lower = lower_bound_L> L[K];            // enhancement parameter 
  real<lower = 0> tau[K];                      // multiply L for hosp
  real<lower = 0> w[K];                        // shape parameter 
  real<lower = 0> lc[C,K];                     // 50% symp protection
  real<lower = 0> alpha[K];                    // reduction in titre for protection against hosp
  real<lower = 0> omega;                       // reduction in titre for protection against MU compared to MO 
  real<lower = 0> kappa;                       // increase in titre for protection against SN compared to MO 
  real beta[J-1];                              // increase in lc50 for protection in younger 
  real<lower = 0, upper = 1> sens ;            // baseline test sensitivity 
  real<lower = 0, upper = 1> spec ;            // baseline test specificity 
}

model {}

generated quantities{
  real p_KJ[K,J] ; 
  real pSP[J]; // prob SP (baseline for likelihood)
  real N_SP[J]; // number SP
  real<lower = 0> C_BVKJRD[B,V,K,J,R,D];
  real<lower = 0> n[B,K,T]; // titres
  real pi_1[B];   // decay rate 1 
  real pi_2;      // decay rate 2
  real<lower = 0> RR_symp[C,V,K,J,T];             
  real<lower = 0> RR_hosp[C,V,K,J,T];
  vector[K] rhoT; // keep code simple by always referring to vector length K
  vector[K] gammaT;
  vector[K] phiT;
  real TpSN[J];              // true SN
  real TpMO[K,J];            // true monotypic  
  real TpMU2[(K+2),J];       // true multitypic
  real TpMU3[K,J];           // true multitypic  
  real TpMU4[J];             // true multitypic  
  real pSN[B,V,J,(T+1)];     // prob SN time plus baseline 
  real pMO[B,V,K,J,(T+1)];   // prob monotypic time plus baseline 
  real pMU2[B,V,(K+2),J,(T+1)]; // prob multitypic time plus baseline 
  real pMU3[B,V,K,J,(T+1)];  // prob multitypic time plus baseline 
  real n_C[C,K,T];           // allow for multitypic titres 
  real L_C[C,K,R];           // MO and MU L = 1 
  real nc50[C,K,J,R];        // nc50
  real lambda[K,T];          // FOI 
  real Inc1[B,V,K,J,T];
  real Inc2[B,V,K,J,T];
  real Inc3[B,V,K,J,T];
  real Inc4[B,V,K,J,T];
  real D1[B,V,K,J,T];
  real D2[B,V,K,J,T];
  real D34[B,V,K,J,T];
  real H1[B,V,K,J,T];
  real H2[B,V,K,J,T];
  real H34[B,V,K,J,T];
  real Di[B,V,K,J,D];
  real H[B,V,K,J,D];
  real pop_VJD[V,J,D];
 
// BIPHASIC TITRES 
 
for(b in 1:B) pi_1[b] = -log(2) / hs[b];
pi_2 = -log(2) / hl ; 

for(b in 1:B)
 for(k in 1:K)
  for(t in 1:T)
   n[b,k,t] = mu[b,k] * (exp(pi_1[b] * time[t] + pi_2 * ts[b]) + exp(pi_2 * time[t] + pi_1[b] * ts[b])) / (exp(pi_1[b] * ts[b]) + exp(pi_2 * ts[b])) ; 

 for(k in 1:K)
  for(t in 1:T){
    n_C[1,k,t] =  n[1,k,t];
    n_C[2,k,t] =  n[2,k,t]; // MO and MU have SP titres 
    n_C[3,k,t] =  n[2,k,t];
  }
  
// MODEL SELECTION 
if(rho_K==1) {
  for(k in 1:K) {
    rhoT[k] = rho[k];
    gammaT[k] = gamma;
    phiT[k] = phi;
  }
} else {
  for(k in 1:K) {
    rhoT[k] = rho[1];
    gammaT[k] = gamma; 
    phiT[k] = phi;
  }
}
// rho is reduction in primary symp compared to sec
// gamma is prob sec symp but gammaT is prob primary (divided by rho)
    for(k in 1:K) gammaT[k] /= rhoT[k]; 
      
// INITIAL CONDITIONS
// this defines the true serostatus populations 

for(k in 1:K)
 for(j in 1:J)
p_KJ[k,j] = p[j]; 

 for(j in 1:J){
    TpSN[j] = (1-p_KJ[1,j]) * (1-p_KJ[2,j]) * (1-p_KJ[3,j]) * (1-p_KJ[4,j]);
    // 1
    TpMO[1,j] = p_KJ[1,j]   * (1-p_KJ[2,j]) * (1-p_KJ[3,j]) * (1-p_KJ[4,j]);
    // 2
    TpMO[2,j] = p_KJ[2,j]   * (1-p_KJ[1,j]) * (1-p_KJ[3,j]) * (1-p_KJ[4,j]);
    // 3
    TpMO[3,j] = p_KJ[3,j]   * (1-p_KJ[1,j]) * (1-p_KJ[2,j]) * (1-p_KJ[4,j]);
    // 4
    TpMO[4,j] = p_KJ[4,j]   * (1-p_KJ[1,j]) * (1-p_KJ[2,j]) * (1-p_KJ[3,j]);
    // 1,2
    TpMU2[1,j] = p_KJ[1,j] * p_KJ[2,j] * (1-p_KJ[3,j]) * (1-p_KJ[4,j]);
    // 1,3
    TpMU2[2,j] = p_KJ[1,j] * p_KJ[3,j] * (1-p_KJ[2,j]) * (1-p_KJ[4,j]);
    // 1,4
    TpMU2[3,j] = p_KJ[1,j] * p_KJ[4,j] * (1-p_KJ[2,j]) * (1-p_KJ[3,j]);
    // 2,3
    TpMU2[4,j] = p_KJ[2,j] * p_KJ[3,j] * (1-p_KJ[1,j]) * (1-p_KJ[4,j]);
    // 2,4
    TpMU2[5,j] = p_KJ[2,j] * p_KJ[4,j] * (1-p_KJ[1,j]) * (1-p_KJ[3,j]);
    // 3,4 
    TpMU2[6,j] = p_KJ[3,j] * p_KJ[4,j] * (1-p_KJ[1,j]) * (1-p_KJ[2,j]);
    // not 1 
    TpMU3[1,j] = p_KJ[2,j] * p_KJ[3,j] * p_KJ[4,j] * (1-p_KJ[1,j]) ;
    // not 2
    TpMU3[2,j] = p_KJ[1,j] * p_KJ[3,j] * p_KJ[4,j] * (1-p_KJ[2,j]) ;
    // not 3
    TpMU3[3,j] = p_KJ[1,j] * p_KJ[2,j] * p_KJ[4,j] * (1-p_KJ[3,j]) ;
    // not 4 
    TpMU3[4,j] = p_KJ[1,j] * p_KJ[2,j] * p_KJ[3,j] * (1-p_KJ[4,j]) ;
    // all 
    TpMU4[j] = p_KJ[1,j] * p_KJ[2,j] * p_KJ[3,j] * p_KJ[4,j];
 }

// this accounts for imperfect test performance 
for(v in 1:V)
 for(j in 1:J)
   for(k in 1:K){
    pSN[1,v,j,1] = spec * TpSN[j] ;
    pSN[2,v,j,1] = (1-spec) * TpSN[j] ;
    pMO[1,v,k,j,1] = (1- sens) * TpMO[k,j] ;
    pMO[2,v,k,j,1] = sens * TpMO[k,j] ;
   }

if(MU_test_SN == 0){ // assume multitypic all classified SP
  for(v in 1:V)
   for(j in 1:J) {
    for(k in 1:6){  // there are 6 MU_2 combinations  
        pMU2[1,v,k,j,1] = 0;
        pMU2[2,v,k,j,1] = TpMU2[k,j] ; 
    }
    for(k in 1:K){  // there are 4 MU_3 combinations 
        pMU3[1,v,k,j,1] = 0; 
        pMU3[2,v,k,j,1] = TpMU3[k,j] ;
    } }
     } else { // allow mulitypic misclassification 
    for(v in 1:V)
     for(j in 1:J){
      for(k in 1:6){  // there are 6 MU_2 combinations  
          pMU2[1,v,k,j,1] = (1- sens) * TpMU2[k,j] ;
          pMU2[2,v,k,j,1] = sens * TpMU2[k,j] ; 
      }
      for(k in 1:K){  // there are 4 MU_3 combinations 
          pMU3[1,v,k,j,1] = (1- sens) * TpMU3[k,j] ; 
          pMU3[2,v,k,j,1] = sens * TpMU3[k,j] ;
      }}}

// probabilities of testing seropositive 
if(MU_test_SN == 0){  
 for(j in 1:J)  pSP[j] = (1-spec) * TpSN[j] + sens * sum(TpMO[ ,j]) + sum(TpMU2[ ,j]) + sum(TpMU3[ ,j]) + TpMU4[j];
} else{
 for(j in 1:J)  pSP[j] = (1-spec) * TpSN[j] + sens * (sum(TpMO[ ,j]) + sum(TpMU2[ ,j]) + sum(TpMU3[ ,j]) + TpMU4[j]) ;
}

 for(j in 1:J) N_SP[j] = pSP[j] * pop_J[j] ; 

// VACCINE RISK RATIO 

// outcome and age-group titre offsets 
real e_beta[2] = {exp(beta[1]), exp(beta[2])};
real e_alpha[C,K]; 

for(k in 1:K)
 for(c in 1:C){
   if(alpha_CK == 0){ 
     e_alpha[c,k] = exp(-alpha[1]);  // mono outcome offset 
     } else if(alpha_CK == 1) {
       e_alpha[c,k] = exp(-alpha[c]);  // serostatus outcome offset 
       } else {
         e_alpha[c,k] = exp(-alpha[k]);  // serotype outcome offset 
         }}
        
if(mono_lc_SN == 1){ // SN, oldest, symp 
  for(k in 1:K) nc50[1,k,3,1] = exp(lc[1,1]); 
   } else if (mono_lc_SN == 0) {
  for(k in 1:K) nc50[1,k,3,1] = exp(lc[1,k]); 
} else if(mono_lc_SN == 2){
  for(k in 1:K)  nc50[1,k,3,1] = exp(lc[2,k]) * exp(kappa); 
}

if(mono_lc_MU == 1){ // MU, oldest, symp
 for(k in 1:K)  nc50[3,k,3,1] = exp(lc[3,1]); 
   } else if(mono_lc_MU == 0) {
 for(k in 1:K) nc50[3,k,3,1] = exp(lc[3,k]); 
} else if (mono_lc_MU == 2){
   for(k in 1:K)  nc50[3,k,3,1] = exp(lc[2,k]) * exp(-omega); 
}


// MO, oldest, symp
for(k in 1:K) nc50[2,k,3,1] = exp(lc[2,k]); 

// age group offsets 
if(include_beta == 1){
  for(c in 1:C)
   for(k in 1:K){
   nc50[c,k,1,1] =   nc50[c,k,3,1] * e_beta[1]; // youngest
   nc50[c,k,2,1] =   nc50[c,k,3,1] * e_beta[2]; // middle 
 } 
} else if(include_beta == 3) {
  for(c in 1:C)
   for(k in 1:K){
   nc50[c,k,1,1] =   nc50[c,k,3,1] * e_beta[1]; // youngest
   nc50[c,k,2,1] =   nc50[c,k,3,1]; // middle 
   }
 } else {
  for(c in 1:C)
   for(k in 1:K){
   nc50[c,k,1,1] = nc50[c,k,3,1];
   nc50[c,k,2,1] = nc50[c,k,3,1];   
  }
}

// hosp 
for(c in 1:C)
 for(k in 1:K) 
  for(j in 1:J)
    nc50[c,k,j,2] =   nc50[c,k,j,1] * e_alpha[c,k];

if(include_beta == 2){
  for(c in 1:C)
   for(k in 1:K){
     nc50[c,k,1,1] =   nc50[c,k,1,1] * e_beta[1]; 
     nc50[c,k,1,2] =   nc50[c,k,1,2] * e_beta[2]; 
   }
}

// Enhancement 
for(k in 1:K){
  if(enhancement == 1){ 
   if(L_K == 0){
   L_C[1,k,1] = 1 + L[1]; // symp 
   L_C[1,k,2] = 1 + L[1] * ((tau_K==0)?tau[1]:tau[k]);
  } else if (L_K == 1){
   L_C[1,k,1] = 1 + L[k]; // symp 
   L_C[1,k,2] = 1 + L[k] * ((tau_K==0)?tau[1]:tau[k]);
 } 
 } else if (enhancement == 0){ // no SN enhancement 
    for(r in 1:R)  L_C[1,k,r] = 1 ; 
 } 

 for(r in 1:R) { // no enhancement if seropositive 
        L_C[2,k,r] = 1; // MO
        L_C[3,k,r] = 1; // MU
      }
   }
 
// Risk ratios - w can be mono, serostatus or serotype dependent  
for(c in 1:C)
 for(k in 1:K)
  for(j in 1:J)
   for(t in 1:T){
     if(w_CK == 0){
     RR_symp[c,2,k,j,t] = L_C[c,k,1] / (1 +  (n_C[c,k,t] / nc50[c,k,j,1])^w[1]) ;
     RR_hosp[c,2,k,j,t] = L_C[c,k,2] / (1 +  (n_C[c,k,t] / nc50[c,k,j,2])^w[1]) ;
     } else if(w_CK == 1){
     RR_symp[c,2,k,j,t] = L_C[c,k,1] / (1 +  (n_C[c,k,t] / nc50[c,k,j,1])^w[c]) ;
     RR_hosp[c,2,k,j,t] = L_C[c,k,2] / (1 +  (n_C[c,k,t] / nc50[c,k,j,2])^w[c]) ;  
     } else{
     RR_symp[c,2,k,j,t] = L_C[c,k,1] / (1 +  (n_C[c,k,t] / nc50[c,k,j,1])^w[k]) ;
     RR_hosp[c,2,k,j,t] = L_C[c,k,2] / (1 +  (n_C[c,k,t] / nc50[c,k,j,2])^w[k]) ;  
     }
     RR_symp[c,1,k,j,t] =  1 ;
     RR_hosp[c,1,k,j,t] =  1 ;
    }

// FOI  
for(k in 1:K){
  for(t in 1:12)  lambda[k,t] = exp(-lambda_D[k,1]);
  for(t in 13:18) lambda[k,t] = exp(-lambda_D[k,2]);
  for(t in 19:24) lambda[k,t] = exp(-lambda_D[k,3]);
  for(t in 25:36) lambda[k,t] = exp(-lambda_D[k,4]);
  for(t in 37:48) lambda[k,t] = exp(-lambda_D[k,5]);
  for(t in 49:54) lambda[k,t] = exp(-lambda_D[k,6]);
}

// SURVIVAL MODEL - prob of surviving each time point without infection 
for(b in 1:B)
 for(v in 1:V)
  for(j in 1:J)
   for(t in 1:T) pSN[b,v,j,(t+1)] = lambda[1,t]*lambda[2,t]*lambda[3,t]*lambda[4,t] * pSN[b,v,j,t];
   
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(j in 1:J) {
    for(t in 1:T) pMO[b,v,k,j,(t+1)] = lambda[1,t]*lambda[2,t]*lambda[3,t]*lambda[4,t]/lambda[k,t]* pMO[b,v,k,j,t] ; 
    for(t in (HI+1):T) pMO[b,v,k,j,(t+1)] += (1 - gammaT[k] * RR_symp[1,v,k,j,(t-HI)]) * (1- lambda[k,(t-HI)]) * pSN[b,v,j,(t-HI)] ;
    }  
   
for(b in 1:B)
 for(v in 1:V)
   for(j in 1:J) {
      for(t in 1:T){
          // MU_12
          pMU2[b,v,1,j,(t+1)] = lambda[3,t] * lambda[4,t] *  pMU2[b,v,1,j,t] ; 
           // MU_13
          pMU2[b,v,2,j,(t+1)] = lambda[2,t] * lambda[4,t] *  pMU2[b,v,2,j,t] ; 
           // MU_14
          pMU2[b,v,3,j,(t+1)] = lambda[2,t] * lambda[3,t] *  pMU2[b,v,3,j,t] ; 
          // MU_23
          pMU2[b,v,4,j,(t+1)] = lambda[1,t] * lambda[4,t] *  pMU2[b,v,4,j,t] ; 
          // MU_24
          pMU2[b,v,5,j,(t+1)] = lambda[1,t] * lambda[3,t] *  pMU2[b,v,5,j,t] ; 
          // MU_34
          pMU2[b,v,6,j,(t+1)] = lambda[1,t] * lambda[2,t] *  pMU2[b,v,6,j,t] ; 
      }
      for(t in (HI+1):T) {
          // MU_12
          pMU2[b,v,1,j,(t+1)] +=  (1 - gammaT[1] * rhoT[1] * RR_symp[2,v,1,j,(t-HI)]) * (1-lambda[1,(t-HI)]) * pMO[b,v,2,j,(t-HI)] + (1 - gammaT[2] * rhoT[2] * RR_symp[2,v,2,j,(t-HI)]) * (1-lambda[2,(t-HI)]) * pMO[b,v,1,j,(t-HI)] ;
          // MU_13
          pMU2[b,v,2,j,(t+1)] +=  (1 - gammaT[1] * rhoT[1] * RR_symp[2,v,1,j,(t-HI)]) * (1-lambda[1,(t-HI)]) * pMO[b,v,3,j,(t-HI)] + (1 - gammaT[3] * rhoT[3] * RR_symp[2,v,3,j,(t-HI)]) * (1-lambda[3,(t-HI)]) * pMO[b,v,1,j,(t-HI)] ;
          // MU_14
          pMU2[b,v,3,j,(t+1)] +=  (1 - gammaT[1] * rhoT[1] * RR_symp[2,v,1,j,(t-HI)]) * (1-lambda[1,(t-HI)]) * pMO[b,v,4,j,(t-HI)] + (1 - gammaT[4] * rhoT[4] * RR_symp[2,v,4,j,(t-HI)]) * (1-lambda[4,(t-HI)]) * pMO[b,v,1,j,(t-HI)] ;
          // MU_23
          pMU2[b,v,4,j,(t+1)] +=  (1 - gammaT[2] * rhoT[2] * RR_symp[2,v,2,j,(t-HI)]) * (1-lambda[2,(t-HI)]) * pMO[b,v,3,j,(t-HI)] + (1 - gammaT[3] * rhoT[3] * RR_symp[2,v,3,j,(t-HI)]) * (1-lambda[3,(t-HI)]) * pMO[b,v,2,j,(t-HI)] ;
          // MU_24
          pMU2[b,v,5,j,(t+1)] +=  (1 - gammaT[2] * rhoT[2] * RR_symp[2,v,2,j,(t-HI)]) * (1-lambda[2,(t-HI)]) * pMO[b,v,4,j,(t-HI)] + (1 - gammaT[4] * rhoT[4] * RR_symp[2,v,4,j,(t-HI)]) * (1-lambda[4,(t-HI)]) * pMO[b,v,2,j,(t-HI)] ;
          // MU_34
          pMU2[b,v,6,j,(t+1)] +=  (1 - gammaT[3] * rhoT[3] * RR_symp[2,v,3,j,(t-HI)]) * (1-lambda[3,(t-HI)]) * pMO[b,v,4,j,(t-HI)] + (1 - gammaT[4] * rhoT[4] * RR_symp[2,v,4,j,(t-HI)]) * (1-lambda[4,(t-HI)]) * pMO[b,v,3,j,(t-HI)] ;
        }}

for(b in 1:B)
 for(v in 1:V)
   for(j in 1:J){
     for(t in 1:T){ 
        // MU3 -1 
        pMU3[b,v,1,j,(t+1)] = lambda[1,t]  *  pMU3[b,v,1,j,t] ; 
        // MU3 -2 
        pMU3[b,v,2,j,(t+1)] = lambda[2,t]  *  pMU3[b,v,2,j,t] ; 
        // MU3 -3
        pMU3[b,v,3,j,(t+1)] = lambda[3,t]  *  pMU3[b,v,3,j,t] ; 
        // MU3 -4
        pMU3[b,v,4,j,(t+1)] = lambda[4,t]  *  pMU3[b,v,4,j,t] ; 
     }
     for(t in (HI+1):T) {
        // MU3 -1 
        pMU3[b,v,1,j,(t+1)] += (1 - phiT[4] * gammaT[4] *  RR_symp[3,v,4,j,(t-HI)]) * (1-lambda[4,(t-HI)]) * pMU2[b,v,4,j,(t-HI)] + (1 - phiT[3] * gammaT[3] *  RR_symp[3,v,3,j,(t-HI)]) * (1-lambda[3,(t-HI)]) * pMU2[b,v,5,j,(t-HI)] + (1 - phiT[2] * gammaT[2] *  RR_symp[3,v,2,j,(t-HI)]) * (1-lambda[2,(t-HI)]) * pMU2[b,v,6,j,(t-HI)] ; 
        // MU3 -2 
        pMU3[b,v,2,j,(t+1)] += (1 - phiT[4] * gammaT[4] *  RR_symp[3,v,4,j,(t-HI)]) * (1-lambda[4,(t-HI)]) * pMU2[b,v,2,j,(t-HI)] + (1 - phiT[3] * gammaT[3] *  RR_symp[3,v,3,j,(t-HI)]) * (1-lambda[3,(t-HI)]) * pMU2[b,v,3,j,(t-HI)] + (1 - phiT[1] * gammaT[1] *  RR_symp[3,v,1,j,(t-HI)]) * (1-lambda[1,(t-HI)]) * pMU2[b,v,6,j,(t-HI)] ; 
        // MU3 -3
        pMU3[b,v,3,j,(t+1)] += (1 - phiT[4] * gammaT[4] *  RR_symp[3,v,4,j,(t-HI)]) * (1-lambda[4,(t-HI)]) * pMU2[b,v,1,j,(t-HI)] + (1 - phiT[2] * gammaT[2] *  RR_symp[3,v,2,j,(t-HI)]) * (1-lambda[2,(t-HI)]) * pMU2[b,v,3,j,(t-HI)] + (1 - phiT[1] * gammaT[1] *  RR_symp[3,v,1,j,(t-HI)]) * (1-lambda[1,(t-HI)]) * pMU2[b,v,5,j,(t-HI)] ;  
        // MU3 -4
        pMU3[b,v,4,j,(t+1)] += (1 - phiT[3] * gammaT[3] *  RR_symp[3,v,3,j,(t-HI)]) * (1-lambda[3,(t-HI)]) * pMU2[b,v,1,j,(t-HI)] + (1 - phiT[2] * gammaT[2] *  RR_symp[3,v,2,j,(t-HI)]) * (1-lambda[2,(t-HI)])  * pMU2[b,v,2,j,(t-HI)] + (1 - phiT[1] * gammaT[1] *  RR_symp[3,v,1,j,(t-HI)]) * (1-lambda[1,(t-HI)]) * pMU2[b,v,4,j,(t-HI)] ;  
      }
    }  

// calculate infection incidence, symptomatic, hospitalised 
for(b in 1:B)
  for(v in 1:V)
    for(j in 1:J) {
      for(t in 1:T){
         Inc3[b,v,1,j,t] = (1 - lambda[1,t]) * (pMU2[b,v,4,j,t] + pMU2[b,v,5,j,t] + pMU2[b,v,6,j,t]);  
         Inc3[b,v,2,j,t] = (1 - lambda[2,t]) * (pMU2[b,v,2,j,t] + pMU2[b,v,3,j,t] + pMU2[b,v,6,j,t]);  
         Inc3[b,v,3,j,t] = (1 - lambda[3,t]) * (pMU2[b,v,1,j,t] + pMU2[b,v,3,j,t] + pMU2[b,v,5,j,t]);  
         Inc3[b,v,4,j,t] = (1 - lambda[4,t]) * (pMU2[b,v,1,j,t] + pMU2[b,v,2,j,t] + pMU2[b,v,4,j,t]);
         
         for(k in 1:K){
           Inc1[b,v,k,j,t] = (1 - lambda[k,t]) * pSN[b,v,j,t]; 
           Inc2[b,v,k,j,t] = (1 - lambda[k,t]) * (sum(pMO[b,v, ,j,t]) - pMO[b,v,k,j,t]);
           Inc4[b,v,k,j,t] = (1 - lambda[k,t]) * pMU3[b,v,k,j,t];  
           D1[b,v,k,j,t]  = Inc1[b,v,k,j,t] * gammaT[k] * RR_symp[1,v,k,j,t] ;
           D2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k] * RR_symp[2,v,k,j,t];
           D34[b,v,k,j,t] = (Inc3[b,v,k,j,t] + Inc4[b,v,k,j,t]) * gammaT[k] * phiT[k] * RR_symp[3,v,k,j,t];  
           H1[b,v,k,j,t]  = Inc1[b,v,k,j,t] * gammaT[k] * delta[k] * RR_hosp[1,v,k,j,t];
           if(include_eps == 0){
           H2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k]  * delta[k]  * RR_hosp[2,v,k,j,t];
           } else{
           H2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k]  * delta[k] * epsilon * RR_hosp[2,v,k,j,t]; 
           }
           H34[b,v,k,j,t] = (Inc3[b,v,k,j,t] + Inc4[b,v,k,j,t]) * gammaT[k] * phiT[k] * delta[k] * RR_hosp[3,v,k,j,t] ;   
         }}}

// Aggregate to match published time points e.g. 1-12 months 
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(j in 1:J){
   Di[b,v,k,j,1] = sum(D1[b,v,k,j,1:12])  + sum(D2[b,v,k,j,1:12])   ;
   H[b,v,k,j,1]  = sum(H1[b,v,k,j,1:12])  + sum(H2[b,v,k,j,1:12])   ;
   Di[b,v,k,j,2] = sum(D1[b,v,k,j,13:18]) + sum(D2[b,v,k,j,13:18])  ;
   H[b,v,k,j,2]  = sum(H1[b,v,k,j,13:18]) + sum(H2[b,v,k,j,13:18])  ;
   Di[b,v,k,j,3] = sum(D1[b,v,k,j,19:24]) + sum(D2[b,v,k,j,19:24])  ;
   H[b,v,k,j,3]  = sum(H1[b,v,k,j,19:24]) + sum(H2[b,v,k,j,19:24])  ;
   Di[b,v,k,j,4] = sum(D1[b,v,k,j,25:36]) + sum(D2[b,v,k,j,25:36])  ;
   H[b,v,k,j,4]  = sum(H1[b,v,k,j,25:36]) + sum(H2[b,v,k,j,25:36])  ;
   Di[b,v,k,j,5] = sum(D1[b,v,k,j,37:48]) + sum(D2[b,v,k,j,37:48])  ;
   H[b,v,k,j,5]  = sum(H1[b,v,k,j,37:48]) + sum(H2[b,v,k,j,37:48])  ;
   Di[b,v,k,j,6] = sum(D1[b,v,k,j,49:54]) + sum(D2[b,v,k,j,49:54])  ;
   H[b,v,k,j,6]  = sum(H1[b,v,k,j,49:54]) + sum(H2[b,v,k,j,49:54])  ;
   }

if(MU_symp == 1){ // if post-secondary can be symptomatic 
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(j in 1:J){
   Di[b,v,k,j,1] += sum(D34[b,v,k,j,1:12])  ;
   H[b,v,k,j,1]  += sum(H34[b,v,k,j,1:12])  ;
   Di[b,v,k,j,2] += sum(D34[b,v,k,j,13:18])  ;
   H[b,v,k,j,2]  += sum(H34[b,v,k,j,13:18])  ;
   Di[b,v,k,j,3] += sum(D34[b,v,k,j,19:24])  ;
   H[b,v,k,j,3]  += sum(H34[b,v,k,j,19:24])  ;
   Di[b,v,k,j,4] += sum(D34[b,v,k,j,25:36])  ;
   H[b,v,k,j,4]  += sum(H34[b,v,k,j,25:36])  ;
   Di[b,v,k,j,5] += sum(D34[b,v,k,j,37:48])  ;
   H[b,v,k,j,5]  += sum(H34[b,v,k,j,37:48])  ;
   Di[b,v,k,j,6] += sum(D34[b,v,k,j,49:54])  ;
   H[b,v,k,j,6]  += sum(H34[b,v,k,j,49:54])  ;
   }}
   
for(v in 1:V)
 for(j in 1:J)
  for(d in 1:D)
    pop_VJD[v,j,d] = sum(pop[ ,v,j,d]) ; // pop not by serostatus

// Cases 
for(b in 1:B) 
 for(v in 1:V) 
  for(k in 1:K)
   for(j in 1:J) 
    for(d in 1:D){
          C_BVKJRD[b,v,k,j,1,d] = Di[b,v,k,j,d] * pop_VJD[v,j,d]; 
          C_BVKJRD[b,v,k,j,2,d] = H[b,v,k,j,d]  * pop_VJD[v,j,d]; 
   }
}