# code data

#setwd("C:/Users/dominik_deffner/Documents/GitHub/Experiment-Simulation")  #MPI

setwd("C:/Users/Dominik/Documents/GitHub/Experiment-Simulation")   #Laptop

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


m0 <- stan( file="ewa_model1.stan" , data=dat , chains=1 )

m1 <- stan( file="ewa_model2.stan" , data=dat , chains=1, iter = 500 )

m2 <- stan( file="ewa_model3.stan" , data=dat , chains=1, iter = 500 )

m3 <- stan( file="ewa_model4.stan" , data=dat , chains=1, iter = 500 )

m5 <- stan( file="ewa_model5.stan" , data=dat , chains=1, iter = 500 )

m6 <- stan( file="ewa_model6.stan" , data=dat , chains=1, iter = 1000 )

m7 <- stan( file="ewa_model7.stan" , data=dat , chains=1, iter = 500 )


sum <- precis(m6, depth = 3)
plot(inv_logit(sum$mean[3:22]), type="l", ylim=c(0,1), ylab=expression(sigma), xlab="Experience", lwd=2)
polygon(c(1:20,20:1), c(inv_logit(sum$`94.5%`[3:22]), rev(inv_logit(sum$`5.5%`[3:22]))), col=alpha("blue",alpha = 0.2), border = NA, ylim=c(0,1))
curve(0.7 *exp(-0.1*(x-1)), 1,20, ylim=c(0,1), add = TRUE, lty=2, ylab = "", xlab = "", lwd=2)