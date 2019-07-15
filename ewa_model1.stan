// Individual learning only

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
}

parameters{
real phi;
real<lower=0> L;
}

model{

matrix[N_id,4] A; // attraction matrix

phi ~ beta(2,2);
L ~ exponential(1);

// initialize attraction scores
	for ( i in 1:N_id ) A[i,1:4] = rep_vector(0,4)';

// loop over choices

for ( i in 1:N ) {
vector[4] pay;
vector[4] p;

if ( new_farm[i]==1 ) A[id[i],1:4] = rep_vector(0,4)';

// first, what is log-prob of observed choice

p = softmax( L*A[id[i],1:4]' );
Choice[i] ~ categorical( p );

// second, update attractions conditional on observed choice

pay[1:4] = rep_vector(0,4);
pay[ Choice[i] ] = Payoff[i];
A[ id[i] , 1:4 ] = ( (1-phi)*to_vector(A[ id[i] , 1:4 ]) + phi*pay)';

}//i
}
