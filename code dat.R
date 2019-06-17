# code data

setwd("C:/Users/dominik_deffner/Documents/GitHub/Experiment-Simulation")  #MPI

#setwd("C:/Users/Dominik/Documents/GitHub/Experiment-Simulation")   #Laptop

dat <- as.list(Result)
dat$farm <- as.integer( ceiling( dat$trial / 25 ) )
dat$new_farm <- sapply( 1:nrow(Result) , function(i) {
	if ( i>1 ) {
		return( ifelse(dat$farm[i]!=dat$farm[i-1],1,0) )
	} else {
		return(0)
	}
} )

dat$group_id <- dat$group.id
dat$group.id <- NULL

# prep id's
dat$sid <- dat$id
dat$id <- (dat$Session_ID - 1)*8 + dat$id
dat$N_id <- max(dat$id)
dat$N <- nrow(Result)

nmat <- cbind( dat$n1 , dat$n2 , dat$n3 , dat$n4 )
choice_models <- cbind(dat$Choice1, dat$Choice2, dat$Choice3)
age_models    <- cbind(dat$Age1, dat$Age2, dat$Age3)

choice_models[is.na(choice_models)] <- 0
age_models[is.na(age_models)] <- 0

dat$nmat <- nmat
dat$choice_models <- choice_models
dat$age_models <- age_models

dat$expmat <- matrix(nrow = 20, ncol = 20)
for (i in 1:20) {
  for (j in 1:20) {
    dat$expmat[i,j] <- abs(i-j)
  }
}

library(rethinking)
library(scales)


#m0 <- stan( file="ewa_model1.stan" , data=dat , chains=1 )

#m1 <- stan( file="ewa_model2.stan" , data=dat , chains=1, iter = 500 )  #now with social cues

#m2 <- stan( file="ewa_model3.stan" , data=dat , chains=1, iter = 500 )  #convex combination of age and frequency

#m3 <- stan( file="ewa_model4.stan" , data=dat , chains=1, iter = 500 )

#m5 <- stan( file="ewa_model5.stan" , data=dat , chains=1, iter = 4000 )

#m6 <- stan( file="ewa_model6.stan" , data=dat , chains=1, iter = 1000 )

#m7 <- stan( file="ewa_model7.stan" , data=dat , chains=1, iter = 500 )

#m8 <- stan( file="ewa_model8.stan" , data=dat , chains=1, iter = 1000 )


m10 <- stan( file="ewa_model10.stan" , data=dat , chains=4,cores = 4, iter = 10000)



m <- m10samp
m <- extract.samples(m10)

#GP curves for sigma
m_Mat_sigma <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    m_Mat_sigma[j,i] <- inv_logit(mean(m$mean_sigma) + mean(m$z_GP[,1,j]) + mean(m$dev_sigma[,j,i]))
    
  }
}
Mean_M_mat_sigma <- apply(m_Mat_sigma,2,mean)

#GP curves for f

m_Mat_f <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    m_Mat_f[j,i] <- exp(mean(m$mean_f) + mean(m$z_GP[,3,j]) + mean(m$dev_f[,j,i]))
    
  }
}
Mean_m_Mat_f <- apply(m_Mat_f,2,mean)


#GP curves for beta

m_Mat_beta <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    m_Mat_beta[j,i] <-mean(m$mean_beta) + mean(m$z_GP[,2,j]) + mean(m$dev_beta[,j,i])
    
  }
}
Mean_m_Mat_beta <- apply(m_Mat_beta,2,mean)


#GP curves for kappa

m_Mat_kappa <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    m_Mat_kappa[j,i] <- inv_logit(mean(m$mean_kappa) + mean(m$z_GP[,4,j]) + mean(m$dev_beta[,j,i]))
    
  }
}
Mean_m_Mat_kappa <- apply(m_Mat_kappa,2,mean)

#Plot it

par(mfrow= c(2,2), 
    mar=c(2,4,2,1), 
    oma=c(3,0,0,0))

plot(Mean_M_mat_sigma, ylim=c(0,1), ylab=expression(sigma),type="n", xlab="Experience", lwd=3, col="black")
for (j in 1:dat$N_id) {
  lines(m_Mat_sigma[j,], type="l", ylim=c(0,1), ylab=expression(sigma), xlab="Experience", lwd=0.5,col="lightgrey")
}
lines(Mean_M_mat_sigma, type="l", ylim=c(0,1), ylab=expression(sigma), xlab="Experience", lwd=2, col="black")



plot(Mean_M_mat, ylim=c(0,2), ylab="f",type="n", xlab="Experience", lwd=3, col="black")
for (j in 1:dat$N_id) {
  lines(m_Mat_f[j,], type="l", ylim=c(0,2), ylab=expression(sigma), xlab="Experience", lwd=0.5,col="lightgrey")
}
lines(Mean_m_Mat_f, type="l", ylim=c(0,2), ylab=expression(sigma), xlab="Experience", lwd=2, col="black")
abline(h=1, lty=2)


plot(Mean_m_Mat_beta, ylim=c(-3,3), ylab=expression(beta),type="n", xlab="Experience", lwd=3, col="black")
for (j in 1:dat$N_id) {
  lines(m_Mat_beta[j,], type="l", ylim=c(-3,3), ylab=expression(sigma), xlab="Experience", lwd=0.5,col="lightgrey")
}
lines(Mean_m_Mat_beta, type="l", ylim=c(-3,3), ylab=expression(sigma), xlab="Experience", lwd=2, col="black")
abline(h=0, lty=2)


plot(Mean_m_Mat_kappa, ylim=c(0,1), ylab=expression(kappa),type="n", xlab="Experience", lwd=3, col="black")
for (j in 1:dat$N_id) {
  lines(m_Mat_kappa[j,], type="l", ylim=c(0,1), ylab=expression(sigma), xlab="Experience", lwd=0.5,col="lightgrey")
}
lines(Mean_m_Mat_kappa, type="l", ylim=c(0,1), ylab=expression(sigma), xlab="Experience", lwd=2, col="black")

mtext(side = 1, line = 1.2 , "Experience", outer = TRUE)
