//Try to add varying effects on GP

functions{
    matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real sq_sigma) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + sq_sigma;
          for (j in (i + 1):N) {
            K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sq_alpha + sq_sigma;
        return K;
    }
}

data{
int N;
int N_id;
int Session_ID[N];
int id[N];
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

//Frequency of each option chosen
matrix[N,4] nmat;

//Experience distance matrix
matrix[20,20] expmat;

// age vars

matrix[N,3] age_models;
matrix[N,3] choice_models;
}

parameters{
real<lower=0,upper=1> phi;
real<lower=0> L;
real<lower=0> f;                    // strength of conformity
real<lower=0,upper=1> kappa;    // weight of age bias
real beta;                         // strength of age bias

real mean_sigma; // mean weight of social info for average level of experience

//Gaussian process stuff
real log_etasq;   // max covariance in Gaussian process
real log_rhosq;   // rate of decline
real log_sigmasq; // additional variance of main diagonal

matrix[N_id, 20] dev_sigma;   // Average deviations for levels of experience

// Varying effects clustered on individual
matrix[4,N_id] z_GP;

//[1,] mean_sigma : mean weight of social info
//[2,] etasq: max covariance in Gaussian process
//[3,] rhosq: rate of decline
//[4,] sigmasq: additional variance of main diagonal

vector<lower=0>[4] sigma_ID;       //SD of parameters among individuals
cholesky_factor_corr[4] Rho_ID;

}

transformed parameters{

  matrix[N_id,4] v_GP; // varying effects on stuff
  v_GP = ( diag_pre_multiply( sigma_ID , Rho_ID ) * z_GP )';

}

model{

matrix[N_id,4] A; // attraction matrix
matrix[20, 20] Kmat[N_id];  // Array of N_id 20 x 20 matrices to store person specific covariances

mean_sigma ~ normal(0,1);
phi ~ beta(2,2);
L ~ exponential(1);
f ~ normal(1,0.5);
kappa ~ beta(2,2);
beta ~ normal(0,0.5);
log_rhosq ~ normal(0,1.5);
log_etasq ~ normal(0,1.5);
log_sigmasq ~ normal(0,1.5);

//varying effects
to_vector(z_GP) ~ normal(0,1);
sigma_ID ~ exponential(1);
Rho_ID ~ lkj_corr_cholesky(4);

for ( i in 1:N_id) {
   Kmat[i] = cov_GPL2(expmat, exp(log_etasq + v_GP[i,2]) , exp(log_rhosq + v_GP[i,3]), exp(log_sigmasq + v_GP[i,4]));

    dev_sigma[i,] ~ multi_normal(rep_vector(0,20) , Kmat[i]);
}

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
// First part adds person-specific deviation to global mean sigma
sigma = inv_logit((mean_sigma + v_GP[id[i],1]) + dev_sigma[id[i], Experience[i]]);

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
