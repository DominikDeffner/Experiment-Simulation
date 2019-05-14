# code data

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
dat$uid <- (dat$Session_ID - 1)*8 + dat$id
dat$N_id <- max(dat$uid)
dat$N <- nrow(Result)

nmat <- cbind( dat$n1 , dat$n2 , dat$n3 , dat$n4 )
choice_models <- cbind(dat$Choice1, dat$Choice2, dat$Choice3)
age_models    <- cbind(dat$Age1, dat$Age2, dat$Age3)

choice_models[is.na(choice_models)] <- 0
age_models[is.na(age_models)] <- 0

dat$nmat <- nmat
dat$choice_models <- choice_models
dat$age_models <- age_models

library(rethinking)

m0 <- stan( file="ewa_model1.stan" , data=dat , chains=1 )

m1 <- stan( file="ewa_model2.stan" , data=dat , chains=1 )

m2 <- stan( file="ewa_model3.stan" , data=dat , chains=1 )
