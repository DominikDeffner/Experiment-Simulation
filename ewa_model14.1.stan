// Monotonic effect on sigma


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

real logit_sigma_first;
real logit_sigma_last;

real Gauss_beta_first;
real Gauss_beta_last;

real log_f_first;
real log_f_last;

real logit_kappa_first;
real logit_kappa_last;

real logit_phi;
real log_L;


// monotonic effect on Experience
simplex[19] delta_sigma;
simplex[19] delta_beta;
simplex[19] delta_f;
simplex[19] delta_kappa;


// Varying effects clustered on individual
matrix[10,N_id] z_ID;
vector<lower=0>[10] sigma_ID;       //SD of parameters among individuals
cholesky_factor_corr[10] Rho_ID;

}

transformed parameters{
matrix[N_id,10] v_ID; // varying effects on stuff
v_ID = ( diag_pre_multiply( sigma_ID , Rho_ID ) * z_ID )';
}


//1  logit_sigma_first;
//2 logit_sigma_last;

//3 beta_first;
//4 beta_last;

//5 log_f_first;
//6 log_f_last;

//7 logit_kappa_first;
//8 logit_kappa_last;

//9 logit_phi;
//10 log_L;


model{

matrix[N_id,4] A; // attraction matrix
vector[20] delta_sigma_container;
vector[20] delta_beta_container;
vector[20] delta_f_container;
vector[20] delta_kappa_container;

logit_phi ~  normal(0,1);
log_L ~  normal(0,1);

logit_sigma_first ~ normal(0,1);
logit_sigma_last ~ normal(0,1);

Gauss_beta_first ~ normal(0,1);
Gauss_beta_last ~ normal(0,1);

logit_kappa_first ~ normal(0,1);
logit_kappa_last ~ normal(0,1);

log_f_first ~ normal(0,1);
log_f_last ~ normal(0,1);


// monotonic Experience terms
delta_sigma ~ dirichlet( rep_vector(0.1,19) );
delta_sigma_container = append_row( 0 , delta_sigma);

delta_beta ~ dirichlet( rep_vector(0.1,19) );
delta_beta_container = append_row( 0 , delta_beta);

delta_f ~ dirichlet( rep_vector(0.1,19) );
delta_f_container = append_row( 0 , delta_f);

delta_kappa ~ dirichlet( rep_vector(0.1,19) );
delta_kappa_container = append_row( 0 , delta_kappa);

// varying effects
to_vector(z_ID) ~ normal(0,1);
sigma_ID ~ exponential(1);
Rho_ID ~ lkj_corr_cholesky(4);

// initialize attraction scores

for ( i in 1:N_id ) A[i,1:4] = rep_vector(0,4)';

// loop over choices

for ( i in 1:N ) {

  vector[4] pay;
  vector[4] pA;
  vector[4] pC;
  vector[4] pS;
  vector[4] p;

  real sigma_first;
  real sigma_last;
  real sigma;

  real beta_first;
  real beta_last;
  real beta;

  real kappa_first;
  real kappa_last;
  real kappa;

  real f_first;
  real f_last;
  real f;

  real L;
  real phi;


  if ( new_farm[i]==1 ) A[id[i],1:4] = rep_vector(0,4)';

  // first, what is log-prob of observed choice
  L =  exp(log_L + v_ID[id[i],10]);

  pA = softmax( L*A[id[i],1:4]' );

  // second compute conformity thing

  if ( sum(nmat[i,:])==0 ) {
    p = pA;

  } else {

    // conformity
    f_first = exp( log_f_first + v_ID[id[i],5] );
    f_last = exp( log_f_last + v_ID[id[i],6] );
    f = f_first - (f_first-f_last) * sum( delta_f_container[1:Experience[i]]);

    for ( j in 1:4 ) pC[j] = nmat[i,j]^f;

    pC = pC / sum(pC);

    //age bias

    beta_first = Gauss_beta_first + v_ID[id[i],3] ;
    beta_last =  Gauss_beta_last  + v_ID[id[i],4] ;
    beta = beta_first - (beta_first-beta_last) * sum( delta_beta_container[1:Experience[i]]);

    for ( an_option in 1:4 ){

      pS[an_option] = 0;

      for ( a_model in 1:3 ) {

        if ( choice_models[i,a_model]==an_option )

        pS[an_option] = pS[an_option] + exp(beta*age_models[i,a_model]);

      }
    }

    pS = pS / sum(pS);

    // combine everything
    sigma_first = inv_logit( logit_sigma_first + v_ID[id[i],1] );
    sigma_last = inv_logit( logit_sigma_last + v_ID[id[i],2] );
    sigma = sigma_first - (sigma_first-sigma_last) * sum( delta_sigma_container[1:Experience[i]]);

    kappa_first = inv_logit( logit_kappa_first + v_ID[id[i],7] );
    kappa_last  = inv_logit( logit_kappa_last  + v_ID[id[i],8] );
    kappa = kappa_first - (kappa_first-kappa_last) * sum( delta_kappa_container[1:Experience[i]]);

    p = (1-sigma)*pA + sigma*( (1-kappa)*pC + kappa*pS );

  }

  Choice[i] ~ categorical( p );

  // second, update attractions conditional on observed choice

  pay[1:4] = rep_vector(0,4);
  pay[ Choice[i] ] = Payoff[i];

  phi =  inv_logit(logit_phi + v_ID[id[i],9]);

  A[ id[i] , 1:4 ] = ( (1-phi)*to_vector(A[ id[i] , 1:4 ]) + phi*pay)';


}
//i
}
