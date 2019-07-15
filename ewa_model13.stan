//All functional + varying effects

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

    real log_lambda;
    real logit_phi;

    real a_sigma; // mean weight of social info for average level of experience
    real b_sigma;
    real a_beta;  // mean strength of age bias for average level of experience
    real b_beta;
    real a_f;
    real b_f;
    real a_kappa;
    real b_kappa;


    // Varying effects clustered on individual
    matrix[10,N_id] z_GP;

    //[1,] a_sigma : mean weight of social info
    //[2,] b_sigma: max covariance in Gaussian process

    //[3, ] a_beta
    //[4, ] b_beta

    //[5, ] a_f
    //[6, ] b_f

    //[7, ] a_kappa
    //[8, ] b_kappa

    //[9,] log_lambda
    //[10,] logit_phi

    vector<lower=0>[10] sigma_ID;       //SD of parameters among individuals
    cholesky_factor_corr[10] Rho_ID;

}

transformed parameters{
    matrix[N_id,10] v_GP; // varying effects on stuff
    v_GP = ( diag_pre_multiply( sigma_ID , Rho_ID ) * z_GP )';
}

model{

    matrix[N_id,4] A; // attraction matrix

    logit_phi ~ normal(0,1);
    log_lambda ~ normal(0,1);
    a_sigma ~ normal(0,1);
    b_sigma ~ normal(0,1);
    a_beta ~ normal(0,1);
    b_beta ~ normal(0,1);
    a_f ~ normal(0,1);
    b_f ~ normal(0,1);
    a_kappa ~ normal(0,1);
    b_kappa ~ normal(0,1);

    //varying effects
    to_vector(z_GP) ~ normal(0,1);
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
    real sigma;
    real beta;
    real f;
    real kappa;
    real phi;
    real lambda;

    if ( new_farm[i]==1 ) A[id[i],1:4] = rep_vector(0,4)';

    // first, what is log-prob of observed choice
    lambda = exp(log_lambda + v_GP[id[i],9]);

    pA = softmax( lambda*A[id[i],1:4]' );

    // second compute conformity everything

    if ( sum(nmat[i,:])==0 ) {
        p = pA;
    } else {

    // conformity
        f = exp((a_f + v_GP[id[i],5]) + (b_f + v_GP[id[i],6])*Experience[i]);

        for ( j in 1:4 ) pC[j] = nmat[i,j]^f;

        pC = pC / sum(pC);

    //age bias

    // Global Mean + Person-specific deviation in average beta + person-specific deviation
    // for level of experience

        beta =(a_beta + v_GP[id[i],3]) + (b_beta + v_GP[id[i],4])*Experience[i];

        for ( an_option in 1:4 ){

            pS[an_option] = 0;

        for ( a_model in 1:3 ) {

            if ( choice_models[i,a_model]==an_option )
                pS[an_option] = pS[an_option] + exp(beta*age_models[i,a_model]);
            }
        }

        pS = pS / sum(pS);

        // combine everything
        // Global Mean + Person-specific deviation in average sigma + person-specific deviation
        // for level of experience

        sigma = inv_logit((a_sigma + v_GP[id[i],1]) + (b_sigma+ v_GP[id[i],2])*Experience[i]);
        kappa = inv_logit((a_kappa + v_GP[id[i],7]) + (b_kappa + v_GP[id[i],8])*Experience[i]);

        p = (1-sigma)*pA + sigma*( (1-kappa)*pC + kappa*pS );

    }

    Choice[i] ~ categorical( p );

    // second, update attractions conditional on observed choice
    phi = inv_logit(logit_phi + v_GP[id[i],10]);

    pay[1:4] = rep_vector(0,4);
    pay[ Choice[i] ] = Payoff[i];
    A[ id[i] , 1:4 ] = ( (1-phi)*to_vector(A[ id[i] , 1:4 ]) + phi*pay)';
    }
}
