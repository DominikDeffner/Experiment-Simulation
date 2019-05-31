
// Attempt at simplest Gaussian process

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

matrix[N,4] nmat;
matrix[20,20] expmat;

// age vars

matrix[N,3] age_models;
matrix[N,3] choice_models;
}

parameters{
real<lower=0,upper=1> phi;
real<lower=0> L;
real mean_sigma;           //mean weight of social info
vector[20] dev_sigma;
real<lower=0> f;                    // strength of conformity
real<lower=0,upper=1> kappa;    // weight of age bias
real beta;                         // strength of age bias
<lower=0> etasq; // max covariance in Gaussian process
real<lower=0> rhosq; // rate of decline
real<lower=0> sigmasq; // additional variance of main diagonal

}

model{

matrix[N_id,4] A; // attraction matrix
matrix[20,20] Kmat; // Covariance matrix for parameters

mean_sigma ~ normal(0,1);
phi ~ beta(2,2);
L ~ exponential(1);
f ~ normal(1,0.5);
kappa ~ beta(2,2);
beta ~ normal(0,0.5);
rhosq ~ exponential( 0.5 );
etasq ~ exponential( 2 );
sigmasq ~ exponential(2);


Kmat = cov_GPL2(expmat, etasq, rhosq, sigmasq);
dev_sigma ~ multi_normal( rep_vector(0,20) , Kmat );


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

sigma = inv_logit(mean_sigma + dev_sigma[Experience[i]]);

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
