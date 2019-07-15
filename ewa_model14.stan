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
real<lower=0,upper=1> phi;
real<lower=0> L;

real logit_sigma_first;
real logit_sigma_last;

real<lower=0> f; // strength of conformity
real<lower=0,upper=1> kappa; // weight of age bias
real beta; // strength of age bias

// monotonic effect on Experience
simplex[19] delta_sigma;

// Varying effects clustered on individual
matrix[2,N_id] z_sigma;
vector<lower=0>[2] sigma_ID;       //SD of parameters among individuals
cholesky_factor_corr[2] Rho_ID;

}

transformed parameters{
matrix[N_id,2] v_sigma; // varying effects on stuff
v_sigma = ( diag_pre_multiply( sigma_ID , Rho_ID ) * z_sigma )';
}

model{

matrix[N_id,4] A; // attraction matrix
vector[20] delta_sigma_container;

phi ~ beta(2,2);
L ~ exponential(1);
logit_sigma_first ~ normal(0,1);
logit_sigma_last ~ normal(0,1);

f ~ normal(1,0.5);
kappa ~ beta(2,2);
beta ~ normal(0,0.5);

// monotonic Experience terms
delta_sigma ~ dirichlet( rep_vector(2,19) );
delta_sigma_container = append_row( 0 , delta_sigma);

// varying effects
to_vector(z_sigma) ~ normal(0,1);
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

  if ( new_farm[i]==1 ) A[id[i],1:4] = rep_vector(0,4)';

  // first, what is log-prob of observed choice

  pA = softmax( L*A[id[i],1:4]' );

  // second compute conformity thing

  if ( sum(nmat[i,:])==0 ) {
    p = pA;

  } else {

    // conformity

    for ( j in 1:4 ) pC[j] = nmat[i,j]^f;

    pC = pC / sum(pC);

    //age bias

    for ( an_option in 1:4 ){

      pS[an_option] = 0;

      for ( a_model in 1:3 ) {

        if ( choice_models[i,a_model]==an_option )
        pS[an_option] = pS[an_option] + exp(beta*age_models[i,a_model]);

      }
    }

    pS = pS / sum(pS);

    // combine everything
    sigma_first = inv_logit( logit_sigma_first + v_sigma[id[i],1] );
    sigma_last = inv_logit( logit_sigma_last + v_sigma[id[i],2] );


    sigma = sigma_first - (sigma_first-sigma_last) * sum( delta_sigma_container[1:Experience[i]]);


    p = (1-sigma)*pA + sigma*( (1-kappa)*pC + kappa*pS );

  }

  Choice[i] ~ categorical( p );

  // second, update attractions conditional on observed choice

  pay[1:4] = rep_vector(0,4);
  pay[ Choice[i] ] = Payoff[i];
  A[ id[i] , 1:4 ] = ( (1-phi)*to_vector(A[ id[i] , 1:4 ]) + phi*pay)';


}
//i
}
