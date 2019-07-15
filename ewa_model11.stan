// Individual estimate for everything, doesn't run

data{
  int N;
  int N_id;           // number of unique individuals
  int Session_ID[N];
  int sid[N];          // id within session
  int id[N];         // unique id across all sessions
  int trial[N];
  int group_id[N];
  int Choice[N];
  int Correct[N];
  real Payoff[N];
  int Experience[N];
  int farm[N];
  int new_farm[N];
  int n1[N];
  int n2[N];
  int n3[N];
  int n4[N];

  matrix[N,4] nmat;

  // age vars

  matrix[N,3] age_models;
  matrix[N,3] choice_models;
}

parameters{
  real logit_phi; //mean phi
  real log_L; //mean L

  vector[20] logit_sigma; // mean logit weight of social learning for each level of experience
  vector[20] mean_beta;        // strength of age bias for each level of experience
  vector[20] log_f; // strength of conformity for each level of experience
  vector[20] logit_kappa; // weight of age bias for each level of experience


  // varying effects clustered on individual
  // [1:20] logit_sigma
  // [21:40] mean_beta
  // [41 : 60] log_f
  // [61:80] logit_kappa
  // [81] logit phi
  // [82] log_L

    matrix[82,N_id] z_ind;
    vector<lower=0>[82] sigma_ind; // standard deviation of parameter values among individuals
    cholesky_factor_corr[82] L_Rho_ind;

}

transformed parameters{

    matrix[N_id,82] v_ind; // varying effects on al the parameters (centered)
    v_ind = ( diag_pre_multiply( sigma_ind , L_Rho_ind ) * z_ind )';
}

model{

  matrix[N_id,4] A; // attraction matrix

  logit_phi ~ normal(0,1);
  log_L ~ normal(0,1);
  logit_sigma ~ normal(0,1);
  log_f ~ normal(0,1);
  logit_kappa ~ normal(0,1);
  mean_beta ~ normal(0,1);

  // varying effects
  to_vector(z_ind) ~ normal(0,1);
  sigma_ind ~ exponential(1);
  L_Rho_ind ~ lkj_corr_cholesky(4);

  // initialize attraction scores

  for ( i in 1:N_id ) A[i,1:4] = rep_vector(0,4)';

  // loop over choices

  for ( i in 1:N ) {

    vector[4] pay;
    vector[4] pA;
    vector[4] pC;
    vector[4] pS;
    vector[4] p;
    real sigma;
    real beta;
    real f;
    real kappa;
    real phi;
    real lambda;

    if ( new_farm[i]==1 ) A[id[i],1:4] = rep_vector(0,4)';

    // first, what is log-prob of observed choice

    lambda = exp(log_L + v_ind[id[i],82]);

    pA = softmax( lambda * A[id[i],1:4]' );

    // second compute conformity thing

    if ( sum(nmat[i,:])==0 ) {
      p = pA;

    } else {

      // conformity

      f = exp(log_f[Experience[i]] + v_ind[id[i], (40+Experience[i])]);

      for ( j in 1:4 ) pC[j] = nmat[i,j]^f;

      pC = pC / sum(pC);

      //age bias

      for ( an_option in 1:4 ){

        pS[an_option] = 0;

        for ( a_model in 1:3 ) {

          if ( choice_models[i,a_model]==an_option )

          beta = mean_beta[Experience[i]] + v_ind[id[i], (20+Experience[i])];

          pS[an_option] = pS[an_option] + exp(beta * age_models[i,a_model]);

        }
      }

      pS = pS / sum(pS);

      // combine everything
      sigma = inv_logit(logit_sigma[Experience[i]] + v_ind[id[i],Experience[i]]);
      kappa = inv_logit(logit_kappa[Experience[i]] + v_ind[id[i],(60+Experience[i])]);

      p = (1-sigma)*pA + sigma*( (1-kappa)*pC + kappa*pS );

    }

    Choice[i] ~ categorical( p );

    // second, update attractions conditional on observed choice

    phi = inv_logit(logit_phi + v_ind[id[i],81]);

    pay[1:4] = rep_vector(0,4);
    pay[ Choice[i] ] = Payoff[i];

    A[ id[i] , 1:4 ] = ( (1-phi)*to_vector(A[ id[i] , 1:4 ]) + phi*pay)';

  }
  //i
}
