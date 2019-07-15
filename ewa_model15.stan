// GP on sigma, monotonic on rest


//Function for Gaussian process kernel
functions{
    matrix cov_IDL2(matrix x, real sq_alpha, real sq_rho, real sq_sigma) {
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

    real logit_phi;
    real log_L;

    real mean_sigma; // mean weight of social info for average level of experience

    //Gaussian process stuff; on the log/logit scale to make them positive/between 0-1 in the end
    real log_etasq_sigma;   // max covariance in Gaussian process
    real log_rhosq_sigma;   // rate of decline
    real log_sigmasq_sigma; // additional variance of main diagonal


    real Gauss_beta_first;
    real Gauss_beta_last;

    real log_f_first;
    real log_f_last;

    real logit_kappa_first;
    real logit_kappa_last;

    // GP for sigma
    matrix[N_id, 20] dev_sigma;   // Average deviations for levels of experience

    // monotonic effect on Experience
    simplex[19] delta_beta;
    simplex[19] delta_f;
    simplex[19] delta_kappa;


    // Varying effects clustered on individual
    matrix[12,N_id] z_GP;

    // 1 mean_sigma : mean weight of social info

    // 2 log_etasq_sigma: max covariance in Gaussian process
    // 3 log_rhosq_sigma: rate of decline
    // 4 log_sigmasq_sigma: additional variance of main diagonal

   // 5 Gauss_beta_first
   // 6 Gauss_beta_last

   // 7 log_f_first
   // 8 log_f_last

   // 9 logit_kappa_first
   // 10 logit_kappa_last

    // 11 log_L
    // 12 logit_phi

    vector<lower=0>[12] sigma_ID;       //SD of parameters among individuals
    cholesky_factor_corr[12] Rho_ID;

}

transformed parameters{
    matrix[N_id,12] v_ID; // varying effects on stuff
    v_ID = ( diag_pre_multiply( sigma_ID , Rho_ID ) * z_GP )';
}

model{

    matrix[N_id,4] A; // attraction matrix


    matrix[20, 20] Kmat_sigma[N_id];  // Array of N_id 20 x 20 matrices to store person specific covariances

    vector[20] delta_beta_container;
    vector[20] delta_f_container;
    vector[20] delta_kappa_container;


    mean_sigma ~ normal(0,1);

    logit_phi ~ normal(0,1);
    log_L ~ normal(0,1);

    log_rhosq_sigma ~ normal(-1,1);
    log_etasq_sigma ~ normal(-1,1);
    log_sigmasq_sigma ~ normal(-1,1);

    Gauss_beta_first ~ normal(0,0.5);
    Gauss_beta_last ~ normal(0,0.5);

    logit_kappa_first ~ normal(0,1);
    logit_kappa_last ~ normal(0,1);

    log_f_first ~ normal(0,0.5);
    log_f_last ~ normal(0,0.5);


    // monotonic Experience terms
    delta_beta ~ dirichlet( rep_vector(0.1,19) );
    delta_beta_container = append_row( 0 , delta_beta);

    delta_f ~ dirichlet( rep_vector(0.1,19) );
    delta_f_container = append_row( 0 , delta_f);

    delta_kappa ~ dirichlet( rep_vector(0.1,19) );
    delta_kappa_container = append_row( 0 , delta_kappa);

    //varying effects
    to_vector(z_GP) ~ normal(0,1);
    sigma_ID ~ exponential(1);
    Rho_ID ~ lkj_corr_cholesky(4);

    // Gaussian process fun
    for ( i in 1:N_id) {
       Kmat_sigma[i] = cov_IDL2(expmat, exp(log_etasq_sigma + v_ID[i,2]) , exp(log_rhosq_sigma + v_ID[i,3]), exp(log_sigmasq_sigma + v_ID[i,4]));
       dev_sigma[i,] ~ multi_normal(rep_vector(0,20) , Kmat_sigma[i]);
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

    real beta_first;
    real beta_last;
    real beta;

    real kappa_first;
    real kappa_last;
    real kappa;

    real f_first;
    real f_last;
    real f;

    real lambda;
    real phi;

    if ( new_farm[i]==1 ) A[id[i],1:4] = rep_vector(0,4)';

    // first, what is log-prob of observed choice
    lambda = exp(log_L + v_ID[id[i],11]);

    pA = softmax( lambda*A[id[i],1:4]' );

    // second compute conformity everything

    if ( sum(nmat[i,:])==0 ) {
        p = pA;
    } else {

    // conformity
        f_first = exp( log_f_first + v_ID[id[i],7] );
        f_last = exp( log_f_last + v_ID[id[i],8] );
        f = f_first - (f_first-f_last) * sum( delta_f_container[1:Experience[i]]);

        for ( j in 1:4 ) pC[j] = nmat[i,j]^f;

        pC = pC / sum(pC);

    //age bias

    // Global Mean + Person-specific deviation in average beta + person-specific deviation
    // for level of experience

        beta_first = Gauss_beta_first + v_ID[id[i],5] ;
        beta_last =  Gauss_beta_last  + v_ID[id[i],6] ;
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
        // Global Mean + Person-specific deviation in average sigma + person-specific deviation
        // for level of experience

        sigma = inv_logit(mean_sigma + v_ID[id[i],1] + dev_sigma[id[i], Experience[i]]);

        kappa_first = inv_logit( logit_kappa_first + v_ID[id[i],9] );
        kappa_last  = inv_logit( logit_kappa_last  + v_ID[id[i],10] );

        kappa = kappa_first - (kappa_first-kappa_last) * sum( delta_kappa_container[1:Experience[i]]);
        p = (1-sigma)*pA + sigma*( (1-kappa)*pC + kappa*pS );

    }

    Choice[i] ~ categorical( p );

    // second, update attractions conditional on observed choice
    phi = inv_logit(logit_phi + v_ID[id[i],12]);

    pay[1:4] = rep_vector(0,4);
    pay[ Choice[i] ] = Payoff[i];
    A[ id[i] , 1:4 ] = ( (1-phi)*to_vector(A[ id[i] , 1:4 ]) + phi*pay)';
    }
}
