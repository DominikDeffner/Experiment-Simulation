
library(rethinking)

# Load data

setwd("C:/Users/dominik_deffner/Documents/GitHub/Experiment-Simulation")  

Result <- read_excel("")


#Prepare data and create new variables from simulation output

dat <- as.list(Result)
dat$farm <- as.integer( ceiling( dat$trial / 25 ) )
dat$new_farm <- sapply( 1:nrow(Result) , function(i) {
	if ( i>1 ) {
		return( ifelse(dat$farm[i]!=dat$farm[i-1],1,0) )
	} else {
		return(0)
	}
} )


# prep id's
dat$group_id <- dat$group.id
dat$group.id <- NULL
dat$sid <- dat$id
dat$id <- (dat$Session_ID - 1)*8 + dat$id
dat$N_id <- max(dat$id)
dat$N <- nrow(Result)

#prep matrices with choices and ages
nmat <- cbind( dat$n1 , dat$n2 , dat$n3 , dat$n4 )
choice_models <- cbind(dat$Choice1, dat$Choice2, dat$Choice3)
age_models    <- cbind(dat$Age1, dat$Age2, dat$Age3)

choice_models[is.na(choice_models)] <- 0
age_models[is.na(age_models)] <- 0

dat$nmat <- nmat
dat$choice_models <- choice_models
dat$age_models <- age_models


#Experience matric for Gaussian processes
dat$expmat <- matrix(nrow = 20, ncol = 20)
for (i in 1:20) {
  for (j in 1:20) {
    dat$expmat[i,j] <- abs(i-j)
  }
}



m1 <- stan( file="ewa_model1.stan" , data=dat , chains=1 ) #Individual learning only

m2 <- stan( file="ewa_model2.stan" , data=dat , chains=1, iter = 500 )  #Simplest model with individual and social learning

m3 <- stan( file="ewa_model3.stan" , data=dat , chains=1, iter = 500 )  #Social Learning with frequency and experience bias

m4 <- stan( file="ewa_model4.stan" , data=dat , chains=1, iter = 500 )  #Logistic regression on weight of social learning (sigma); parameters vary by individual

m5 <- stan( file="ewa_model5.stan" , data=dat , chains=1, iter = 4000 )

m6 <- stan( file="ewa_model6.stan" , data=dat , chains=1, iter = 1000 )

m7 <- stan( file="ewa_model7.stan" , data=dat , chains=1, iter = 500 ) # Gaussian process on sigma

m8 <- stan( file="ewa_model8.stan" , data=dat , chains=1, iter = 1000 ) # Gaussian process on sigma with varying effects

m9 <- stan( file="ewa_model9.stan" , data=dat , chains=1, iter = 1000 ) # GP + varying effects on sigma and beta

m10 <- stan( file="ewa_model10.stan" , data=dat , chains=4,cores = 4, iter = 1000)  #Full GP model

m11 <- stan( file="ewa_model11.stan" , data=dat , chains=1,cores = 1, iter = 500)  #Individual estimate for everything

m12 <- stan( file="ewa_model12.stan" , data=dat , chains=1,cores = 1, iter = 500)  #Independent intercepts all over

m13 <- stan( file="ewa_model13.stan" , data=dat , chains=1,cores = 1, iter = 500) # //GP + varying effects on sigma, functions on rest

m14 <- stan( file="ewa_model14.1.stan" , data=dat , chains=1,cores = 1, iter = 500, control = list(adapt_delta=0.8))  #Monotonic effects

m15 <- stan( file="ewa_model15.stan" , data=dat , chains=1,cores = 1, iter = 500, control = list(adapt_delta=0.8))  #GP on sigma, Monotonic effects on rest


