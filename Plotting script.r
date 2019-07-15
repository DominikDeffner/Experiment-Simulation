


#Full GP model (EWA 10)
{
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

#GP curves for lambda

m_Mat_l <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    m_Mat_l[j,i] <- exp(mean(m$mean_L) + mean(m$z_GP[,17,j])) 
    
  }
}
Mean_m_Mat_l <- apply(m_Mat_l,2,mean)


#GP curves for phi

m_Mat_phi <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    m_Mat_phi[j,i] <- inv_logit(mean(m$mean_phi) + mean(m$z_GP[,18,j])) 
    
  }
}
Mean_m_Mat_phi <- apply(m_Mat_phi,2,mean)

#Plot it

par(mfrow= c(3,2), 
    mar=c(2,4,2,1), 
    oma=c(3,0,0,0))

graphics.off()
png("GP_sigma.png", res = 1200, height = 10, width = 12, units = "cm")

plot(Mean_M_mat_sigma, ylim=c(0,1), ylab=expression(sigma),type="n", xlab="Experience", lwd=3, col="black")
for (j in 1:dat$N_id) {
  lines(m_Mat_sigma[j,], type="l", ylim=c(0,1), ylab=expression(sigma), xlab="Experience", lwd=0.5,col="lightgrey")
}
lines(Mean_M_mat_sigma, type="l", ylim=c(0,1), ylab=expression(sigma), xlab="Experience", lwd=2, col="black")
curve(0.9 * exp(-1 * (x-1)) ,1,20,add = TRUE, lty=2, lwd=2)
curve(0.4 * exp(-0 * (x-1)) ,1,20,add = TRUE, lty=2,lwd=2)
dev.off()


plot(Mean_m_Mat_f, ylim=c(0,2), ylab="f",type="n", xlab="Experience", lwd=3, col="black")
for (j in 1:dat$N_id) {
  lines(m_Mat_f[j,], type="l", ylim=c(0,2), ylab=expression(sigma), xlab="Experience", lwd=0.5,col="lightgrey")
}
lines(Mean_m_Mat_f, type="l", ylim=c(0,2), ylab=expression(sigma), xlab="Experience", lwd=2, col="black")
abline(h=1, lty=2)
curve(3 * exp(-1 * x) + 1 ,1,20,add = TRUE, lty=2,lwd=2)


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

plot(Mean_m_Mat_l, ylim=c(0,1), ylab=expression(lambda),type="n", xlab="Experience", lwd=3, col="black")
for (j in 1:dat$N_id) {
  lines(m_Mat_l[j,], type="l", ylim=c(0,1), ylab=expression(sigma), xlab="Experience", lwd=0.5,col="lightgrey")
}
lines(Mean_m_Mat_l, type="l", ylim=c(0,1), ylab=expression(sigma), xlab="Experience", lwd=2, col="black")

plot(Mean_m_Mat_phi, ylim=c(0,1), ylab=expression(phi),type="n", xlab="Experience", lwd=3, col="black")
for (j in 1:dat$N_id) {
  lines(m_Mat_phi[j,], type="l", ylim=c(0,1), ylab=expression(sigma), xlab="Experience", lwd=0.5,col="lightgrey")
}
lines(Mean_m_Mat_phi, type="l", ylim=c(0,1), ylab=expression(sigma), xlab="Experience", lwd=2, col="black")


mtext(side = 1, line = 1.2 , "Experience", outer = TRUE)


#GP parameters
#rho
par(mfrow= c(2,2), 
    mar=c(2,4,2,1), 
    oma=c(3,0,0,0))
dens(exp(m$log_rhosq_sigma), xlim=c(0,0.1))
dens(exp(rnorm(25000, -1,1)), add=TRUE, lty=2)
dens(exp(m$log_rhosq_f),xlim=c(0,5))
dens(exp(rnorm(25000, -1,1)),add=TRUE,lty=2,xlim=c(0,5))

dens(exp(m$log_rhosq_beta),xlim=c(0,5))
dens(exp(rnorm(25000, -1,1)),add=TRUE,lty=2,xlim=c(0,5))

dens(exp(m$log_rhosq_kappa),xlim=c(0,5))
dens(exp(rnorm(25000, -1,1)),add=TRUE,lty=2,xlim=c(0,5))

#etasq

par(mfrow= c(2,2), 
    mar=c(2,4,2,1), 
    oma=c(3,0,0,0))
dens(exp(m$log_etasq_sigma),xlim=c(0,5))
dens(exp(rnorm(25000, -1,1)), add=TRUE, lty=2)
dens(exp(m$log_etasq_f),xlim=c(0,5))
dens(exp(rnorm(25000, -1,1)),add=TRUE,lty=2,xlim=c(0,5))

dens(exp(m$log_etasq_beta),xlim=c(0,5))
dens(exp(rnorm(25000, -1,1)),add=TRUE,lty=2,xlim=c(0,5))

dens(exp(m$log_etasq_kappa),xlim=c(0,5))
dens(exp(rnorm(25000, -1,1)),add=TRUE,lty=2,xlim=c(0,5))


#sigmasq

par(mfrow= c(2,2), 
    mar=c(2,4,2,1), 
    oma=c(3,0,0,0))
dens(exp(m$log_sigmasq_sigma),xlim=c(0,5))
dens(exp(rnorm(25000, -1,1)), add=TRUE, lty=2)
dens(exp(m$log_sigmasq_f),xlim=c(0,5))
dens(exp(rnorm(25000, -1,1)),add=TRUE,lty=2,xlim=c(0,5))

dens(exp(m$log_sigmasq_beta),xlim=c(0,5))
dens(exp(rnorm(25000, -1,1)),add=TRUE,lty=2,xlim=c(0,5))

dens(exp(m$log_sigmasq_kappa),xlim=c(0,5))
dens(exp(rnorm(25000, -1,1)),add=TRUE,lty=2,xlim=c(0,5))
}

#Individual estimates for sigma, binary for rest (EWA 11.1)
{
#Curve for sigma
m_Mat_sigma <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    m_Mat_sigma[j,i] <- inv_logit(mean(m$logit_sigma[,i]) + mean(m$z_ind[,i,j]))
    
  }
}
Mean_M_mat_sigma <- apply(m_Mat_sigma,2,mean)


# curves for f

m_Mat_f <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    if (i==1){
      m_Mat_f[j,i] <- exp(mean(m$log_f[,1]) + mean(m$z_ind[,23,j]))
    }else
      m_Mat_f[j,i] <- exp(mean(m$log_f[,2]) + mean(m$z_ind[,24,j]))
    
  }
}
Mean_m_Mat_f <- apply(m_Mat_f,2,mean)




m_Mat_beta <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    if (i==1){
    m_Mat_beta[j,i] <-mean(m$mean_beta[,1]) + mean(m$z_ind[,21,j])
    }else
    m_Mat_beta[j,i] <-mean(m$mean_beta[,2]) + mean(m$z_ind[,22,j])
    
  }
}
Mean_m_Mat_beta <- apply(m_Mat_beta,2,mean)




m_Mat_kappa <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    if (i==1){
      m_Mat_kappa[j,i] <- inv_logit(mean(m$logit_kappa[,1]) + mean(m$z_ind[,25,j]))
    }else
      m_Mat_kappa[j,i] <- inv_logit(mean(m$logit_kappa[,2]) + mean(m$z_ind[,26,j]))
    
  }
}
Mean_m_Mat_kappa <- apply(m_Mat_kappa,2,mean)


#curves for lambda

m_Mat_l <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    m_Mat_l[j,i] <- exp(mean(m$log_L) + mean(m$z_ind[,27,j])) 
    
  }
}
Mean_m_Mat_l <- apply(m_Mat_l,2,mean)


#curves for phi

m_Mat_phi <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    m_Mat_phi[j,i] <- inv_logit(mean(m$logit_phi) + mean(m$z_ind[,27,j])) 
    
  }
}
Mean_m_Mat_phi <- apply(m_Mat_phi,2,mean)

}

#GP for sigma and binary for rest
{
#GP curves for sigma

m_Mat_sigma <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    m_Mat_sigma[j,i] <- inv_logit(mean(m$mean_sigma) + mean(m$v_GP[,j,1]) + mean(m$dev_sigma[,j,i]))
    
  }
}
Mean_M_mat_sigma <- apply(m_Mat_sigma,2,mean)

# curves for f

m_Mat_f <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    if (i==1){
      m_Mat_f[j,i] <- exp(mean(m$log_f[,1]) + mean(m$v_GP[,j,7]))
    }else
      m_Mat_f[j,i] <- exp(mean(m$log_f[,2]) + mean(m$v_GP[,j,8]))
    
  }
}
Mean_m_Mat_f <- apply(m_Mat_f,2,mean)




m_Mat_beta <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    if (i==1){
      m_Mat_beta[j,i] <-mean(m$mean_beta[,1]) + mean(m$v_GP[,j,5])
    }else
      m_Mat_beta[j,i] <-mean(m$mean_beta[,2]) + mean(m$v_GP[,j,6])
    
  }
}
Mean_m_Mat_beta <- apply(m_Mat_beta,2,mean)




m_Mat_kappa <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    if (i==1){
      m_Mat_kappa[j,i] <- inv_logit(mean(m$logit_kappa[,1]) + mean(m$v_GP[,j,9]))
    }else
      m_Mat_kappa[j,i] <- inv_logit(mean(m$logit_kappa[,2]) + mean(m$v_GP[,j,10]))
    
  }
}
Mean_m_Mat_kappa <- apply(m_Mat_kappa,2,mean)


#curves for lambda

m_Mat_l <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    m_Mat_l[j,i] <- exp(mean(m$log_lambda) + mean(m$v_GP[,j,11])) 
    
  }
}
Mean_m_Mat_l <- apply(m_Mat_l,2,mean)


#curves for phi

m_Mat_phi <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    m_Mat_phi[j,i] <- inv_logit(mean(m$logit_phi) + mean(m$v_GP[,j,12])) 
    
  }
}
Mean_m_Mat_phi <- apply(m_Mat_phi,2,mean)
}


#Monotonic effects

{
#Sigma
Delta_sigmas <- c()

for (i in 1:19) {
  Delta_sigmas[i] <- mean(m$delta_sigma[,i])
}
Delta_sigmas <- c(0, Delta_sigmas)

m_Mat_sigma <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    m_Mat_sigma[j,i] <- inv_logit(mean(m$logit_sigma_first) + mean(m$v_ID[,j,1])) -
                                    (inv_logit(mean(m$logit_sigma_first) + mean(m$v_ID[,j,1]))- inv_logit(mean(m$logit_sigma_last) + mean(m$v_ID[,j,2])))*
                                      sum(Delta_sigmas[1:i])
  }
}
Mean_M_mat_sigma <- apply(m_Mat_sigma,2,mean)




#f
Delta_f <- c()

for (i in 1:19) {
  Delta_f[i] <- mean(m$delta_f[,i])
}
Delta_f <- c(0, Delta_f)

m_Mat_f <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    m_Mat_f[j,i] <- exp(mean(m$log_f_first) + mean(m$v_ID[,j,5])) -
                                    (exp(mean(m$log_f_first) + mean(m$v_ID[,j,5])) - exp(mean(m$log_f_last) + mean(m$v_ID[,j,6]))) *
                                    sum(Delta_f[1:i])
  }
}
Mean_m_Mat_f <- apply(m_Mat_f,2,mean)



#beta
Delta_beta <- c()

for (i in 1:19) {
  Delta_beta[i] <- mean(m$delta_beta[,i])
}
Delta_beta <- c(0, Delta_beta)

m_Mat_beta <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    m_Mat_beta[j,i] <- (mean(m$Gauss_beta_first) + mean(m$v_ID[,j,3])) -
      ((mean(m$Gauss_beta_first) + mean(m$v_ID[,j,3])) - (mean(m$Gauss_beta_last) + mean(m$v_ID[,j,4]))) *
      sum(Delta_beta[1:i])
  }
}
Mean_m_Mat_beta <- apply(m_Mat_beta,2,mean)





#kappa
Delta_kappa <- c()

for (i in 1:19) {
  Delta_kappa[i] <- mean(m$delta_kappa[,i])
}
Delta_kappa <- c(0, Delta_kappa)

m_Mat_kappa <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    m_Mat_kappa[j,i] <- exp(mean(m$logit_kappa_first) + mean(m$v_ID[,j,7])) -
      (exp(mean(m$logit_kappa_first) + mean(m$v_ID[,j,7])) - exp(mean(m$logit_kappa_last) + mean(m$v_ID[,j,8]))) *
      sum(Delta_kappa[1:i])
  }
}
Mean_m_Mat_kappa <- apply(m_Mat_kappa,2,mean)


#curves for lambda

m_Mat_l <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    m_Mat_l[j,i] <- exp(mean(m$log_L) + mean(m$v_ID[,j,10])) 
    
  }
}
Mean_m_Mat_l <- apply(m_Mat_l,2,mean)


#curves for phi

m_Mat_phi <- matrix(0, nrow = dat$N_id, ncol = 20)
for (j in 1:dat$N_id) {
  for (i in 1:20) {
    m_Mat_phi[j,i] <- inv_logit(mean(m$logit_phi) + mean(m$v_ID[,j,9])) 
    
  }
}
Mean_m_Mat_phi <- apply(m_Mat_phi,2,mean)



}


#GP + monotonic effects

{
  
  #GP curves for sigma
  m_Mat_sigma <- matrix(0, nrow = dat$N_id, ncol = 20)
  for (j in 1:dat$N_id) {
    for (i in 1:20) {
      m_Mat_sigma[j,i] <- inv_logit(mean(m$mean_sigma) + mean(m$v_ID[,j,1]) + mean(m$dev_sigma[,j,i]))
      
    }
  }
  Mean_M_mat_sigma <- apply(m_Mat_sigma,2,mean)
  
  
  
  #f
  Delta_f <- c()
  
  for (i in 1:19) {
    Delta_f[i] <- mean(m$delta_f[,i])
  }
  Delta_f <- c(0, Delta_f)
  
  m_Mat_f <- matrix(0, nrow = dat$N_id, ncol = 20)
  for (j in 1:dat$N_id) {
    for (i in 1:20) {
      m_Mat_f[j,i] <- exp(mean(m$log_f_first) + mean(m$v_ID[,j,7])) -
        (exp(mean(m$log_f_first) + mean(m$v_ID[,j,7])) - exp(mean(m$log_f_last) + mean(m$v_ID[,j,8]))) *
        sum(Delta_f[1:i])
    }
  }
  Mean_m_Mat_f <- apply(m_Mat_f,2,mean)
  
  
  
  #beta
  Delta_beta <- c()
  
  for (i in 1:19) {
    Delta_beta[i] <- mean(m$delta_beta[,i])
  }
  Delta_beta <- c(0, Delta_beta)
  
  m_Mat_beta <- matrix(0, nrow = dat$N_id, ncol = 20)
  for (j in 1:dat$N_id) {
    for (i in 1:20) {
      m_Mat_beta[j,i] <- (mean(m$Gauss_beta_first) + mean(m$v_ID[,j,5])) -
        ((mean(m$Gauss_beta_first) + mean(m$v_ID[,j,5])) - (mean(m$Gauss_beta_last) + mean(m$v_ID[,j,6]))) *
        sum(Delta_beta[1:i])
    }
  }
  Mean_m_Mat_beta <- apply(m_Mat_beta,2,mean)
  
  
  
  
  
  #kappa
  Delta_kappa <- c()
  
  for (i in 1:19) {
    Delta_kappa[i] <- mean(m$delta_kappa[,i])
  }
  Delta_kappa <- c(0, Delta_kappa)
  
  m_Mat_kappa <- matrix(0, nrow = dat$N_id, ncol = 20)
  for (j in 1:dat$N_id) {
    for (i in 1:20) {
      m_Mat_kappa[j,i] <- exp(mean(m$logit_kappa_first) + mean(m$v_ID[,j,9])) -
        (exp(mean(m$logit_kappa_first) + mean(m$v_ID[,j,9])) - exp(mean(m$logit_kappa_last) + mean(m$v_ID[,j,10]))) *
        sum(Delta_kappa[1:i])
    }
  }
  Mean_m_Mat_kappa <- apply(m_Mat_kappa,2,mean)
  
  
  #curves for lambda
  
  m_Mat_l <- matrix(0, nrow = dat$N_id, ncol = 20)
  for (j in 1:dat$N_id) {
    for (i in 1:20) {
      m_Mat_l[j,i] <- exp(mean(m$log_L) + mean(m$v_ID[,j,11])) 
      
    }
  }
  Mean_m_Mat_l <- apply(m_Mat_l,2,mean)
  
  
  #curves for phi
  
  m_Mat_phi <- matrix(0, nrow = dat$N_id, ncol = 20)
  for (j in 1:dat$N_id) {
    for (i in 1:20) {
      m_Mat_phi[j,i] <- inv_logit(mean(m$logit_phi) + mean(m$v_ID[,j,12])) 
      
    }
  }
  Mean_m_Mat_phi <- apply(m_Mat_phi,2,mean)
  
  
  
}



#Plot it all 


graphics.off()

png("ExperiencePlotsGPMon.png", res = 1200, height = 15, width = 12, units = "cm")


par(mfrow= c(3,2), 
    mar=c(2,4,2,1), 
    oma=c(3,0,0,0))

plot(Mean_M_mat_sigma, ylim=c(0,1), ylab=expression(sigma),type="n", xlab="Experience", lwd=3, col="black")
for (j in 1:dat$N_id) {
  lines(m_Mat_sigma[j,], type="l", ylim=c(0,1), ylab=expression(sigma), xlab="Experience", lwd=0.5,col="lightgrey")
}
lines(Mean_M_mat_sigma, type="l", ylim=c(0,1), ylab=expression(sigma), xlab="Experience", lwd=2, col="black")



plot(Mean_m_Mat_f, ylab="f",type="n", xlab="Experience", lwd=3, col="black", ylim=c(0,3))
for (j in 1:dat$N_id) {
  lines(m_Mat_f[j,], type="l",  ylab=expression(sigma), xlab="Experience", lwd=0.5,col="lightgrey")
}
lines(Mean_m_Mat_f, type="l", ylab=expression(sigma), xlab="Experience", lwd=2, col="black")
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


plot(Mean_m_Mat_l, ylim=c(0,1), ylab=expression(lambda),type="n", xlab="Experience", lwd=3, col="black")
for (j in 1:dat$N_id) {
  lines(m_Mat_l[j,], type="l", ylim=c(0,1), ylab=expression(sigma), xlab="Experience", lwd=0.5,col="lightgrey")
}
lines(Mean_m_Mat_l, type="l", ylim=c(0,1), ylab=expression(sigma), xlab="Experience", lwd=2, col="black")

plot(Mean_m_Mat_phi, ylim=c(0,1), ylab=expression(phi),type="n", xlab="Experience", lwd=3, col="black")
for (j in 1:dat$N_id) {
  lines(m_Mat_phi[j,], type="l", ylim=c(0,1), ylab=expression(sigma), xlab="Experience", lwd=0.5,col="lightgrey")
}
lines(Mean_m_Mat_phi, type="l", ylim=c(0,1), ylab=expression(sigma), xlab="Experience", lwd=2, col="black")


mtext(side = 1, line = 1.2 , "Experience", outer = TRUE)

dev.off()
