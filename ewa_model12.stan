//GP + varying effects on sigma
// logistic on dev_kappa
//log on f
//gaussian on beta


//Function for Gaussian process kernel
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

    real log_lambda;
    real logit_phi;

    real mean_sigma; // mean weight of social info for average level of experience
    real a_beta;  // mean strength of age bias for average level of experience
    real b_beta;
    real a_f;
    real b_f;
    real a_kappa;
    real b_kappa;


    //Gaussian process stuff; on the log/logit scale to make them positive/between 0-1 in the end
    real log_etasq_sigma;   // max covariance in Gaussian process
    real log_rhosq_sigma;   // rate of decline
    real log_sigmasq_sigma; // additional variance of main diagonal


    matrix[N_id, 20] dev_sigma;   // Average deviations for levels of experience


    // Varying effects clustered on individual
    matrix[12,N_id] z_GP;

    //[1,] mean_sigma : mean weight of social info

    //[2,] log_etasq_sigma: max covariance in Gaussian process
    //[3,] log_rhosq_sigma: rate of decline
    //[4,] log_sigmasq_sigma: additional variance of main diagonal

    //[5, ] a_beta
    //[6, ] b_beta

    //[7, ] a_f
    //[8, ] b_f

    //[9, ] a_kappa
    //[10, ] b_kappa

    //[11,] log_lambda
    //[12,] logit_phi

    vector<lower=0>[12] sigma_ID;       //SD of parameters among individuals
    cholesky_factor_corr[12] Rho_ID;

}

transformed parameters{
    matrix[N_id,12] v_GP; // varying effects on stuff
    v_GP = ( diag_pre_multiply( sigma_ID , Rho_ID ) * z_GP )';
}

model{

    matrix[N_id,4] A; // attraction matrix

    matrix[20, 20] Kmat_sigma[N_id];  // Array of N_id 20 x 20 matrices to store person specific covariances

    mean_sigma ~ normal(0,1);

    logit_phi ~ normal(0,1);
    log_lambda ~ normal(0,1);

    a_beta ~ normal(0,1);
    b_beta ~ normal(0,1);
    a_f ~ normal(0,1);
    b_f ~ normal(0,1);
    a_kappa ~ normal(0,1);
    b_kappa ~ normal(0,1);



    log_rhosq_sigma ~ normal(-1,1);
    log_etasq_sigma ~ normal(-1,1);
    log_sigmasq_sigma ~ normal(-1,1);


    //varying effects
    to_vector(z_GP) ~ normal(0,1);
    sigma_ID ~ exponential(1);
    Rho_ID ~ lkj_corr_cholesky(4);

    // Gaussian process fun
    for ( i in 1:N_id) {
       Kmat_sigma[i] = cov_GPL2(expmat, exp(log_etasq_sigma + v_GP[i,2]) , exp(log_rhosq_sigma + v_GP[i,3]), exp(log_sigmasq_sigma + v_GP[i,4]));
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
    real beta;
    real f;
    real kappa;
    real phi;
    real lambda;

    if ( new_farm[i]==1 ) A[id[i],1:4] = rep_vector(0,4)';

    // first, what is log-prob of observed choice
    lambda = exp(log_lambda + v_GP[id[i],11]);

    pA = softmax( lambda*A[id[i],1:4]' );

    // second compute conformity everything

    if ( sum(nmat[i,:])==0 ) {
        p = pA;
    } else {

    // conformity
        f = exp((a_f + v_GP[id[i],7]) + (b_f + v_GP[id[i],8])*Experience[i]);

        for ( j in 1:4 ) pC[j] = nmat[i,j]^f;

        pC = pC / sum(pC);

    //age bias

    // Global Mean + Person-specific deviation in average beta + person-specific deviation
    // for level of experience

        beta =(a_beta + v_GP[id[i],5]) + (b_beta + v_GP[id[i],6])*Experience[i];

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

        sigma = inv_logit(mean_sigma + v_GP[id[i],1] + dev_sigma[id[i], Experience[i]]);
        kappa = inv_logit((a_kappa + v_GP[id[i],9]) + (b_kappa + v_GP[id[i],10])*Experience[i]);

        p = (1-sigma)*pA + sigma*( (1-kappa)*pC + kappa*pS );

    }

    Choice[i] ~ categorical( p );

    // second, update attractions conditional on observed choice
    phi = inv_logit(logit_phi + v_GP[id[i],12]);

    pay[1:4] = rep_vector(0,4);
    pay[ Choice[i] ] = Payoff[i];
    A[ id[i] , 1:4 ] = ( (1-phi)*to_vector(A[ id[i] , 1:4 ]) + phi*pay)';
    }
}
