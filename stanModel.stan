data {
  int<lower=1> n_days;
  int<lower=1> n_pre;
  int<lower=1> n_variant;
  int<lower=1> t_boostYoung;
  int<lower=1> t_boostYoungClose;

  int M;
  vector[M] N;
  real NTot;
  matrix[M,M] C;

  vector[M] IHR;
  vector[M] pICU;
  vector[M] pDICU;
  vector[M] pDH;

  vector[M] initialDistCases;
  vector[M] f_init;

  int cases[M, n_days];

  int hosp_cases[M, n_days];
  int ICUAgg[n_days];
  int DeathAgg[n_days];

  int n_seq[n_variant];
  int delta_seq[n_variant];

  matrix[M, n_days+n_pre] nu1;
  matrix[M, n_days+n_pre] nu2;
  real v1[M];
  real v2[M];

  real gA1;
  real gA2;
  real gB1;
  real gB2;
  real sA1;
  real sA2;
  real sB1;
  real sB2;

  real var1;
  real var2;

  real sigmaLatentA;
  real muA;
  real sigmaLatentB;
  real muB;

  int stepsDif;

  real oneDoseAlpha;
  real twoDoseAlpha;
  real oneDoseDelta;
  real twoDoseDelta;

  real sIHRAlpha;
  real sIHRDelta;
}
transformed data {
  real muDH = 0.5;
  real probDH = 0.001;
  real dt = 1.0/stepsDif;

}
parameters {
  matrix<lower=0.01>[M, n_days] beta;
  real<lower=0.0> alpha_inv;
  real<lower = 0.0> fB;
  real<lower = 0.0> phi_inv;

  real I0Tot_b;
  real I0B_b;
  real<lower =0.0> muHosp_inv;
  real<lower=0.0> muDICU_inv;
  real<lower=0.001, upper = 1.0> probICU;
  real<lower=0.001, upper = 1.0> probDICU;
  vector<lower = 0.001, upper = 1.0>[M] fD;
}

transformed parameters{

  real alpha = 1.0/alpha_inv;
  real muHosp = 1.0/muHosp_inv;
  real muDICU = 1.0/muDICU_inv;

  real I0Tot = exp(I0Tot_b);
  real I0B = exp(I0B_b);
  real phi = 1.0/phi_inv;

  vector[M] I0A = I0Tot * initialDistCases ./ fD / sum(initialDistCases ./ fD);

  real incA[n_days+n_pre] = rep_array(0.0, n_days+n_pre);
  real incB[n_days+n_pre] = rep_array(0.0, n_days+n_pre);
  real fractionDMod[n_days+n_pre] = rep_array(0.0, n_days+n_pre);

  matrix[M, n_days+n_pre] incTotAge = rep_matrix(0.0, M, n_days+n_pre);
  matrix[M, n_days+n_pre] incTotAgeA = rep_matrix(0.0, M, n_days+n_pre);
  matrix[M, n_days+n_pre] incTotAgeB = rep_matrix(0.0, M, n_days+n_pre);
  matrix[M, n_days+n_pre] incTotAgeNoVacc = rep_matrix(0.0, M, n_days+n_pre);
  matrix[M, n_days+n_pre] incTotAgeV1 = rep_matrix(0.0, M, n_days+n_pre);
  matrix[M, n_days+n_pre] incTotAgeV2 = rep_matrix(0.0, M, n_days+n_pre);

  matrix[M, n_days+n_pre] SN_t = rep_matrix(0.0, M, n_days+n_pre);
  matrix[M, n_days+n_pre] SV1_t = rep_matrix(0.0, M, n_days+n_pre);
  matrix[M, n_days+n_pre] SV2_t = rep_matrix(0.0, M, n_days+n_pre);

  matrix[M, n_days+n_pre] IA_t = rep_matrix(0.0, M, n_days+n_pre);
  matrix[M, n_days+n_pre] IAV1_t = rep_matrix(0.0, M, n_days+n_pre);
  matrix[M, n_days+n_pre] IAV2_t = rep_matrix(0.0, M, n_days+n_pre);

  matrix[M, n_days+n_pre] IB_t = rep_matrix(0.0, M, n_days+n_pre);
  matrix[M, n_days+n_pre] IBV1_t = rep_matrix(0.0, M, n_days+n_pre);
  matrix[M, n_days+n_pre] IBV2_t = rep_matrix(0.0, M, n_days+n_pre);

  matrix[M, n_days+n_pre] Itot_t = rep_matrix(0.0, M, n_days+n_pre);

  real incTot[n_days+n_pre] = rep_array(0.0, n_days+n_pre);
  real incTotD[n_days+n_pre] = rep_array(0.0, n_days+n_pre);
  real incTotVacc1[n_days+n_pre] = rep_array(0.0, n_days+n_pre);
  real incTotVacc2[n_days+n_pre] = rep_array(0.0, n_days+n_pre);
  real incTotNoVacc[n_days+n_pre] = rep_array(0.0, n_days+n_pre);
  real incTotFraccVac[n_days+n_pre] = rep_array(0.0, n_days+n_pre);

  matrix[M, n_days+n_pre] HospAge = rep_matrix(0.0, M, n_days+n_pre);
  real HospTot[n_days+n_pre] = rep_array(0.0, n_days+n_pre);
  real ICUTot[n_days+n_pre] = rep_array(0.0, n_days+n_pre);
  real DeathTot[n_days+n_pre] = rep_array(0.0, n_days+n_pre);

  real incTotVacc1Mod[n_days+n_pre] = rep_array(0.0, n_days+n_pre);
  real incTotVacc2Mod[n_days+n_pre] = rep_array(0.0, n_days+n_pre);
  real incTotNoVaccMod[n_days+n_pre] = rep_array(0.0, n_days+n_pre);


  //define variables
  vector[M] SN; vector[M] SNInter;
  vector[M] SV1; vector[M] SV1Inter;
  vector[M] SV2; vector[M] SV2Inter;
  vector[M] IA; vector[M] IAInter;
  vector[M] IAV1; vector[M] IAV1Inter;
  vector[M] IAV2; vector[M] IAV2Inter;
  vector[M] IB; vector[M] IBInter;
  vector[M] IBV1; vector[M] IBV1Inter;
  vector[M] IBV2; vector[M] IBV2Inter;

  vector[M] EA1; vector[M] EA1Inter;
  vector[M] EA1V1; vector[M] EA1V1Inter;
  vector[M] EA1V2; vector[M] EA1V2Inter;
  vector[M] EB1; vector[M] EB1Inter;
  vector[M] EB1V1; vector[M] EB1V1Inter;
  vector[M] EB1V2; vector[M] EB1V2Inter;

  vector[M] EA2; vector[M] EA2Inter;
  vector[M] EA2V1; vector[M] EA2V1Inter;
  vector[M] EA2V2; vector[M] EA2V2Inter;
  vector[M] EB2; vector[M] EB2Inter;
  vector[M] EB2V1; vector[M] EB2V1Inter;
  vector[M] EB2V2; vector[M] EB2V2Inter;

  vector[M] preHosp; vector[M] preHospInter;
  vector[M] hosp_compICU; vector[M] hosp_compICUInter;
  vector[M] hosp_compDeath; vector[M] hosp_compDeathInter;
  vector[M] ICU_compDeath; vector[M] ICU_compDeathInter;
  vector[M] R; vector[M] RInter;
  vector[M] RNoVacc;
  vector[M] RV1;

  vector[M] b;

  vector[M] A; vector[M] B;
  real ITot = sum(I0A);

  //fix initial condition
  for(i in 1:M){
    b[i] = beta[i,1];

    SN[i] = N[i] - I0A[i] - I0B*I0A[i]/ITot - v1[i] - v2[i] - N[i]*f_init[i]*(1-(v1[i]+v2[i]) / N[i]);
    SNInter[i] = SN[i];
    SV1[i] = v1[i] - v1[i]*f_init[i];
    SV1Inter[i] = SV1[i];
    SV2[i] = v2[i] - v2[i]*f_init[i];
    SV2Inter[i] = SV2[i];

    IA[i] = I0A[i]*1/(1+ sigmaLatentA / muA);
    IAInter[i] = IA[i];
    IAV1[i] = 0;
    IAV1Inter[i] = IAV1[i];
    IAV2[i] = 0;
    IAV2Inter[i] = IAV2[i];

    EA1[i] = I0A[i]*1/(1+sigmaLatentA/muA)/2;
    EA1Inter[i] = EA1[i];
    EA1V1[i] = 0;
    EA1V1Inter[i] = 0;
    EA1V2[i] = 0;
    EA1V2Inter[i] = 0;

    EA2[i] = I0A[i]*1/(1+sigmaLatentA/muA)/2;
    EA2Inter[i] = EA2[i];
    EA2V1[i] = 0;
    EA2V1Inter[i] = 0;
    EA2V2[i] = 0;
    EA2V2Inter[i] = 0;

    IB[i] = I0B*I0A[i]/ITot*1/(1+muB/sigmaLatentB);
    IBInter[i] = IB[i];
    IBV1[i] = 0;
    IBV1Inter[i] = IBV1[i];
    IBV2[i] = 0;
    IBV2Inter[i] = IBV2[i];

    EB1[i] = I0B*I0A[i]/ITot*1/(1+sigmaLatentB/muB)/2;
    EB1Inter[i] = EB1[i];
    EB1V1[i] = 0;
    EB1V1Inter[i] = 0;
    EB1V2[i] = 0;
    EB1V2Inter[i] = 0;

    EB2[i] = I0B*I0A[i]/ITot*1/(1+sigmaLatentB/muB)/2;
    EB2Inter[i] = EB2[i];
    EB2V1[i] = 0;
    EB2V1Inter[i] = 0;
    EB2V2[i] = 0;
    EB2V2Inter[i] = 0;

    R[i] = N[i]*f_init[i];
    RInter[i] = R[i];
    RNoVacc[i] =  N[i]*f_init[i]*(1-(v1[i]+ v2[i]) /N[i]);
    RV1[i] = f_init[i]*v1[i];

    preHosp[i] = 0;
    preHospInter[i] = 0;

    hosp_compICU[i] = 0.0;
    hosp_compICUInter[i] = hosp_compICU[i];
    hosp_compDeath[i] = 0.0;
    hosp_compDeathInter[i] = hosp_compDeath[i];
    ICU_compDeath[i] = 0.0;
    ICU_compDeathInter[i] = ICU_compDeath[i];
  }

  // loop over days
  for(t in 1:(n_days+n_pre)){
    if(t > n_pre){
      for(i in 1:M){
        b[i] = beta[i,t-n_pre];
      }
    }
    else{
      for(i in 1:M){
        b[i] = beta[i,1];
      }
    }

    for(k in 1:stepsDif){
      // calculate A & B
      A = (C*(IA+gA1*IAV1+gA2*IAV2))./N;
      B = (C*fB*(IB+gB1*IBV1+gB2*IBV2))./N;

      SNInter  = SN + (- col(nu1,t) .* (SN) ./ (SN+ RNoVacc) - SN .* b .*(A + B))*dt;
      SV1Inter  = SV1 +  (col(nu1,t) .* (SN) ./ (SN + RNoVacc) - col(nu2,t) .* (SV1) ./ (SV1 + RV1) - SV1 .*b .*(sA1*A + sB1*B))*dt;
      SV2Inter  = SV2 + (col(nu2,t) .* (SV1) ./ (SV1 + RV1)  - SV2 .*b .*(sA2*A + sB2*B))*dt;

      EA1Inter = EA1 +  (SN .*b .*A - 2*sigmaLatentA*EA1)*dt;
      EA2Inter = EA2 +  (EA1 - EA2)*2*sigmaLatentA*dt;
      IAInter = IA + (-IA*muA + 2*sigmaLatentA*EA2)*dt;

      EA1V1Inter = EA1V1 + (SV1 .*b*sA1 .*A - 2*sigmaLatentA*EA1V1)*dt;
      EA2V1Inter = EA2V1 +  (EA1V1 - EA2V1)*2*sigmaLatentA*dt;
      IAV1Inter = IAV1 + (-IAV1*muA + 2*sigmaLatentA*EA2V1)*dt;

      EA1V2Inter = EA1V2 + (SV2 .*b*sA2 .*A - 2*sigmaLatentA*EA1V2)*dt;
      EA2V2Inter = EA2V2 +  (EA1V2 - EA2V2)*2*sigmaLatentA*dt;
      IAV2Inter = IAV2 + (-IAV2*muA +2*sigmaLatentA*EA2V2)*dt;

      EB1Inter = EB1 +  (SN .*b .*B - 2*sigmaLatentB*EB1)*dt;
      EB2Inter = EB2 +  (EB1 - EB2)*2*sigmaLatentB*dt;
      IBInter = IB + (-IB*muB + 2*sigmaLatentB*EB2)*dt;

      EB1V1Inter = EB1V1 + (SV1 .*b*sB1 .*B - 2*sigmaLatentB*EB1V1)*dt;
      EB2V1Inter = EB2V1 +  (EB1V1 - EB2V1)*2*sigmaLatentB*dt;
      IBV1Inter = IBV1 + (-IBV1*muB + 2*sigmaLatentB*EB2V1)*dt;

      EB1V2Inter = EB1V2 + (SV2 .*b*sB2 .*B - 2*sigmaLatentB*EB1V2)*dt;
      EB2V2Inter = EB2V2 +  (EB1V2 - EB2V2)*2*sigmaLatentB*dt;
      IBV2Inter = IBV2 + (-IBV2*muB +2*sigmaLatentB*EB2V2)*dt;

      RInter  = R + muA*(IA+ IAV1 + IAV2)*dt + muB*(IB + IBV1 + IBV2)*dt;
      RNoVacc = RNoVacc + (muA*IA + muB*IB)*dt;
      RV1 = RV1 + (muA*IAV1 + muB*IBV1)*dt;

      preHospInter = preHosp + sIHRAlpha*IHR .*(IA*muA + (1.0-oneDoseAlpha)*IAV1*muA + (1.0-twoDoseAlpha)*IAV2*muA + sIHRDelta*IB*muB + sIHRDelta*(1.0-oneDoseDelta)*IBV1*muB + sIHRDelta*(1.0-twoDoseDelta)*IBV2*muB)*dt;
      preHospInter = preHospInter - preHosp*alpha*dt;

      hosp_compICUInter = hosp_compICU + preHosp .*pICU*probICU*alpha*dt - hosp_compICU*muHosp*dt;

      hosp_compDeathInter = hosp_compDeath + preHosp .*pDH*probDH*alpha*dt - hosp_compDeath*muDH*dt;

      ICU_compDeathInter = ICU_compDeath + hosp_compICU .*pDICU*probDICU*muHosp*dt - ICU_compDeath*muDICU*dt;

      incTotVacc1[t] = incTotVacc1[t] + dot_product((IAV1*muA + IBV1*muB)*dt,fD);
      incTotVacc2[t] = incTotVacc2[t] + dot_product((IAV2*muA + IBV2*muB)*dt,fD);
      incTotNoVacc[t] = incTotNoVacc[t] + dot_product((IA*muA + IB*muB)*dt,fD);

      incA[t] = incA[t] + sum(muA*(IA + IAV1 + IAV2)*dt);
      incB[t] = incB[t] + sum(muB*(IB + IBV1 + IBV2)*dt);
      incTot[t] = incTot[t] + sum((IA*muA + IAV1*muA + IAV2*muA + IB*muB + IBV1*muB + IBV2*muB)*dt);
      incTotD[t] = incTotD[t] + sum((IA*muA + IAV1*muA + IAV2*muA + IB*muB + IBV1*muB + IBV2*muB) .*fD*dt);
      incTotFraccVac[t] = incTotFraccVac[t] + dot_product((IA*muA + IAV1*muA + IAV2*muA + IB*muB + IBV1*muB + IBV2*muB)*dt, fD);


      HospTot[t] = HospTot[t] + sum(preHosp*alpha*dt);
      ICUTot[t] = ICUTot[t] + sum(hosp_compICU*muHosp*dt);
      DeathTot[t] = DeathTot[t] + sum(ICU_compDeath*muDICU*dt);

      incTotAge[:,t] = incTotAge[:,t] + (IA*muA + IAV1*muA + IAV2*muA + IB*muB + IBV1*muB + IBV2*muB)*dt .*fD;
      incTotAgeA[:,t] = incTotAgeA[:,t] + (IA*muA + IAV1*muA + IAV2*muA)*dt .*fD;
      incTotAgeB[:,t] = incTotAgeB[:,t] + (IB*muB + IBV1*muB + IBV2*muB)*dt .*fD;
      incTotAgeNoVacc[:,t] = incTotAgeNoVacc[:,t] + (IA*muA + IB*muB)*dt .*fD;
      incTotAgeV1[:,t] = incTotAgeV1[:,t] + (IAV1*muA + IBV1*muB)*dt .*fD;
      incTotAgeV2[:,t] = incTotAgeV2[:,t] + (IAV2*muA + IBV2*muB)*dt .*fD;


      HospAge[:,t] = HospAge[:,t] + preHosp*alpha*dt;

      //update variables
      SN = SNInter;
      SV1 = SV1Inter;
      SV2 = SV2Inter;

      SN_t[:,t] = SN_t[:,t] + SN/stepsDif;
      SV1_t[:,t] = SV1_t[:,t] + SV1/stepsDif;
      SV2_t[:,t] = SV2_t[:,t] + SV2/stepsDif;

      EA1 = EA1Inter;
      EA1V1 = EA1V1Inter;
      EA1V2 = EA1V2Inter;
      EA2 = EA2Inter;
      EA2V1 = EA2V1Inter;
      EA2V2 = EA2V2Inter;

      EB1 = EB1Inter;
      EB1V1 = EB1V1Inter;
      EB1V2 = EB1V2Inter;
      EB2 = EB2Inter;
      EB2V1 = EB2V1Inter;
      EB2V2 = EB2V2Inter;

      IA = IAInter;
      IAV1 = IAV1Inter;
      IAV2 = IAV2Inter;
      IB = IBInter;
      IBV1 = IBV1Inter;
      IBV2 = IBV2Inter;

      IA_t[:,t] = IA_t[:,t] + IA/stepsDif;
      IAV1_t[:,t] = IAV1_t[:,t] + IAV1/stepsDif;
      IAV2_t[:,t] = IAV2_t[:,t] + IAV2/stepsDif;

      IB_t[:,t] = IB_t[:,t] + IB/stepsDif;
      IBV1_t[:,t] = IBV1_t[:,t] + IBV1/stepsDif;
      IBV2_t[:,t] = IBV2_t[:,t] + IBV2/stepsDif;

      R = RInter;
      preHosp = preHospInter;
      hosp_compICU = hosp_compICUInter;
      hosp_compDeath = hosp_compDeathInter;
      ICU_compDeath = ICU_compDeathInter;
    }
    fractionDMod[t] = incB[t]/incTot[t];
    incTotVacc1Mod[t] = incTotVacc1[t]/incTotFraccVac[t];
    incTotVacc2Mod[t] = incTotVacc2[t]/incTotFraccVac[t];
    incTotNoVaccMod[t] = incTotNoVacc[t]/incTotFraccVac[t];

    Itot_t[:,t] = IA_t[:,t] + IAV1_t[:,t] + IAV2_t[:,t] + IB_t[:,t] + IBV1_t[:,t] + IBV2_t[:,t];
  }
}
model {
  delta_seq ~ binomial(n_seq, fractionDMod[(n_pre+1):(n_pre+n_variant)]);

  for(i in 1:M){
    cases[i,:] ~ neg_binomial_2(incTotAge[i, (n_pre+1):(n_pre+n_days)], phi);
    hosp_cases[i,:]  ~ neg_binomial_2(HospAge[i, (n_pre+1):(n_pre+n_days)], phi);

    for(j in 2:n_days){
      if(j < (t_boostYoung - n_pre) || j > (t_boostYoungClose-n_pre)){
        target += normal_lpdf( beta[i,j] - beta[i,j-1] | 0.0, var1);

      }
      else{
        target += normal_lpdf( beta[i,j] - beta[i,j-1] | 0.0, var2);
      }
    }
  }

  for(i in 2:n_days){
	beta[:,i] ~ cauchy(0.5, 2);
  }

  ICUAgg ~ neg_binomial_2(ICUTot[(n_pre+1):(n_pre+n_days)], phi);
  DeathAgg ~ neg_binomial_2(DeathTot[(n_pre+1):(n_pre+n_days)], phi);

  beta[:,1] ~ exponential(1.0/0.3);

  alpha_inv ~ normal(5.0, 5.0);
  muHosp_inv ~ normal(2.0, 2.0);

  probICU ~ normal(0.24, 0.12/1.96);
  probDICU ~ normal(0.67, 0.10/1.96);

  phi_inv ~ exponential(5.0);

  I0Tot ~ normal(5.2*sum(initialDistCases ./ fD), 5.2*sum(initialDistCases ./ fD)*0.2);
  target += I0Tot_b;
  I0B ~ normal(100.0, 80.0);
  target += I0B_b;

  muHosp_inv ~ normal(5.0, 2.0);
  muDICU_inv ~ normal(10.0, 4.0);

  fD ~ cauchy(0.5, 1.0);
  fB ~ normal(1.5, 0.5);
}

generated quantities {
  real predCases[n_days] = incTot[(n_pre+1):(n_pre+n_days)];
  real predCasesD[n_days] = incTotD[(n_pre+1):(n_pre+n_days)];
  real predCasesB[n_days] = incB[(n_pre+1):(n_pre+n_days)];
  real hospPred[n_days] = HospTot[(n_pre+1):(n_pre+n_days)];
  matrix[M, n_days] predCasesAge = incTotAge[:,(n_pre+1):(n_pre+n_days)];
  matrix[M, n_days] predCasesAgeA = incTotAgeA[:,(n_pre+1):(n_pre+n_days)];
  matrix[M, n_days] predCasesAgeB = incTotAgeB[:,(n_pre+1):(n_pre+n_days)];
  matrix[M, n_days] hospPredAge = HospAge[:,(n_pre+1):(n_pre+n_days)];

  vector[n_days] repTotal;
  vector[n_days] fDTime;

  vector[M] attackRate;
  matrix[M, n_days] repNumber;

  for(i in 1:M){
    attackRate[i] = sum(incTotAge[i,:])/fD[i]/ N[i];
   }


  for(t in 1:n_days){
    repNumber[:,t] = beta[:,t] .* (C*(SN_t[:,t+n_pre] + sA1*SV1_t[:,t+n_pre] + sA2*SV2_t[:,t+n_pre]))./N .* (IA_t[:,t+n_pre] + gA1*IAV1_t[:,t+n_pre] + gA2*IAV2_t[:,t+n_pre]) ./ Itot_t[:,t+n_pre]/muA;
    repNumber[:,t] = repNumber[:,t] + fB*beta[:,t] .* (C*(SN_t[:,t+n_pre] + sB1*SV1_t[:,t+n_pre] + sB2*SV2_t[:,t+n_pre]))./N .* (IB_t[:,t+n_pre] + gB1*IBV1_t[:,t+n_pre] + gB2*IBV2_t[:,t+n_pre]) ./ Itot_t[:,t+n_pre]/muB;

    repTotal[t] = sum(beta[:,t] .* (C*(SN_t[:,t+n_pre] + sA1*SV1_t[:,t+n_pre] + sA2*SV2_t[:,t+n_pre]))./N .* (IA_t[:,t+n_pre] + gA1*IAV1_t[:,t+n_pre] + gA2*IAV2_t[:,t+n_pre]) / sum(Itot_t[:,t+n_pre]))/muA;
    repTotal[t] = repTotal[t] + sum(fB*beta[:,t] .* (C*(SN_t[:,t+n_pre] + sB1*SV1_t[:,t+n_pre] + sB2*SV2_t[:,t+n_pre]))./N .* (IB_t[:,t+n_pre] + gB1*IBV1_t[:,t+n_pre] + gB2*IBV2_t[:,t+n_pre]) / sum(Itot_t[:,t+n_pre]))/muB;

    fDTime[t] = sum(fD .*Itot_t[:,t+n_pre])/sum(Itot_t[:,t+n_pre]);
  }



}
