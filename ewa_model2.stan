// Simplest model with individual and social learning

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
}


parameters{
   real phi;
   real<lower=0> L;
   real<lower=0,upper=1> sigma;
   real<lower=0> f;
}

model{

matrix[N_id,4] A; // attraction matrix

phi ~ beta(2,2);
L ~ exponential(1);
sigma ~ beta(2,2);
f ~ normal(1,0.5);

// initialize attraction scores

for ( i in 1:N_id ) A[i,1:4] = rep_vector(0,4)';

// loop over choices

for ( i in 1:N ) {
vector[4] pay;
vector[4] pA;
vector[4] pS;
vector[4] p;

if ( new_farm[i]==1 ) A[id[i],1:4] = rep_vector(0,4)';

// first, what is log-prob of observed choice

pA = softmax( L*A[id[i],1:4]' );

// second compute conformity thing

if ( sum(nmat[i,:])==0 ) {
p = pA;
} else {
for ( j in 1:4 ) pS[j] = nmat[i,j]^f;
pS = pS / sum(pS);

// combine the two

p = (1-sigma)*pA + sigma*pS;
}

Choice[i] ~ categorical( p );

// second, update attractions conditional on observed choice

pay[1:4] = rep_vector(0,4);
pay[ Choice[i] ] = Payoff[i];
A[ id[i] , 1:4 ] = ( (1-phi)*to_vector(A[ id[i] , 1:4 ]) + phi*pay)';

}
//i
}
