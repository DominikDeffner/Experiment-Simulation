// convex combo of conformity + age bias
// now social learning varies by Experience (time in region)

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
  vector[20] beta;        // strength of age bias for each level of experience
  vector<lower=0>[20] f; // strength of conformity for each level of experience
  vector<lower=0,upper=1>[20] kappa; // weight of age bias for each level of experience


  // varying effects clustered on individual
  // [1:20] logit_sigma
  // [21] logit_phi
  // [22] log_L
  // [23:42] beta

    matrix[22,N_id] z_ind;
    vector<lower=0>[22] sigma_ind; // standard deviation of parameter values among individuals
    cholesky_factor_corr[22] L_Rho_ind;

}

transformed parameters{

    matrix[N_id,22] v_ind; // varying effects on al the parameters (centered)
    v_ind = ( diag_pre_multiply( sigma_ind , L_Rho_ind ) * z_ind )';
}

model{

  matrix[N_id,4] A; // attraction matrix

  logit_phi ~ normal(0,1);
  log_L ~ normal(0,1);
  logit_sigma ~ normal(0,1);
  f ~ normal(1,0.5);
  kappa ~ beta(2,2);
  beta ~ normal(0,0.5);

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
    real phi;
    real L;

    if ( new_farm[i]==1 ) A[id[i],1:4] = rep_vector(0,4)';

    // first, what is log-prob of observed choice

    L = exp((log_L + v_ind[id[i],22]));

    pA = softmax( L * A[id[i],1:4]' );

    // second compute conformity thing

    if ( sum(nmat[i,:])==0 ) {
      p = pA;

    } else {

      // conformity

      for ( j in 1:4 ) pC[j] = nmat[i,j]^f[Experience[i]];

      pC = pC / sum(pC);

      //age bias

      for ( an_option in 1:4 ){

        pS[an_option] = 0;

        for ( a_model in 1:3 ) {

          if ( choice_models[i,a_model]==an_option )
          pS[an_option] = pS[an_option] + exp(beta[Experience[i]]  * age_models[i,a_model]);

        }
      }

      pS = pS / sum(pS);

      // combine everything
      sigma = inv_logit(logit_sigma[Experience[i]] + v_ind[id[i],Experience[i]]);

      p = (1-sigma)*pA + sigma*( (1-kappa[Experience[i]])*pC + kappa[Experience[i]]*pS );

    }

    Choice[i] ~ categorical( p );

    // second, update attractions conditional on observed choice

    phi = inv_logit(logit_phi + v_ind[id[i],21]);

    pay[1:4] = rep_vector(0,4);
    pay[ Choice[i] ] = Payoff[i];

    A[ id[i] , 1:4 ] = ( (1-phi)*to_vector(A[ id[i] , 1:4 ]) + phi*pay)';

  }
  //i
}
