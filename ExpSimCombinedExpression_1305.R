

# Simple social learning experiment with migration
# use simple reinforcement recursion:
# A_i,t = (1-phi)*A_i,t-1 + phi*P_i,t
# where
# A_i,t is latent attraction to option i at time t
# phi is individual specific updating parameter
# P_i is payoff from option i at time t
# choice probability given by softmax
# Pr(1|A_t) = exp( L*A_1,t ) / Sum( exp(L*A_t) )
# where L is a scaling parameter to add noise to choice.
# As L -> 0, choice becomes noisy.
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))


library(rethinking) # for utilities
library(Rlab)

 Sim_fct <- function(Tmax, Group_size, N_group, N_sessions, L, Payoff_Better, Payoff_Worse, Hard_SD = 2,
                     Easy_SD = 0.5, phi,sigma, f, b){
 
 #Create big output matrix
 Result_Overall <- c()
 
 N_part= Group_size*N_group
 
 for (session in 1:N_sessions) {
   
   
   #Create simulated participants divided into two groups who have attraction scores A_i
   Homunculi <- data.frame(id=1:8, L=NA, phi=NA, sigma=NA, f=NA, b=NA, Group=c(rep(1,4), rep(2,4)), A1=NA, A2=NA, A3= NA, A4=NA)
   
   
   #Assign individual-specific parameter values to include some variation
   
   Homunculi$L <-     rtruncnorm(8, a=0, b=+Inf, mean=L, sd=0.02) 
   Homunculi$phi <-   rtruncnorm(8, a=0, b=1, mean=phi, sd=0.05)
   Homunculi$sigma <- sigma#rtruncnorm(8, a=0, b=1, mean=sigma, sd=0.05)
   Homunculi$f <-     rtruncnorm(8, a=0, b=+Inf, mean=f, sd=0.01)
   Homunculi$b <-     rtruncnorm(8, a=-Inf, b=+Inf, mean=b, sd=0.01)
   
   
   #Create output matrix to record choices and payoffs participants received and Attraction scores at each point in time
   
   Result <- data.frame(Session_ID=session,id=c(rep(1,100),rep(2,100),rep(3,100),rep(4,100),
                             rep(5,100),rep(6,100),rep(7,100),rep(8,100)),
                        trial= rep(1:Tmax, N_part), group.id=NA, Choice=NA, Correct=NA, Payoff=NA, Experience=NA, TimeSinceChange=NA, Exp_Stable=NA, SD_Payoff=NA,
                        n1=NA, n2=NA, n3=NA, n4=NA, 
                        Choice1=NA,Choice2=NA,Choice3=NA,
                        Age1=NA, Age2=NA, Age3=NA
                        )
   
   #Define payoff schedule, change every 25 rounds
   Hard_Phases <- sample(1:4, 2)
   
   Times_Env_Change <- c(25,50,75)
   
   if (1 %in% Hard_Phases){
     Result$SD_Payoff[which(Result$trial %in% 1:25)] <- Hard_SD
   } else {
     Result$SD_Payoff[which(Result$trial %in% 1:25)] <- Easy_SD
   }
   
   if (2 %in% Hard_Phases){
     Result$SD_Payoff[which(Result$trial %in% 26:50)] <- Hard_SD
   } else {
     Result$SD_Payoff[which(Result$trial %in% 26:50)] <- Easy_SD
   }
   
   if (3 %in% Hard_Phases){
     Result$SD_Payoff[which(Result$trial %in% 51:75)] <- Hard_SD
   } else {
     Result$SD_Payoff[which(Result$trial %in% 51:75)] <- Easy_SD
   }
   
   if (4 %in% Hard_Phases){
     Result$SD_Payoff[which(Result$trial %in% 76:100)] <- Hard_SD
   } else {
     Result$SD_Payoff[which(Result$trial %in% 76:100)] <- Easy_SD
   }
   
   
   #4 different options
   Crops <- sample(1:4)
   
   #Which crop is optimal in both groups in all 4 phases
   
   #Define Optimal choice for all 100 rounds
   Better_choices <- matrix(nrow = Tmax, ncol = N_group)
   
   Better_choices[,1] <- c(rep(Crops[1],25), rep(Crops[3],25),rep(Crops[2],25), rep(Crops[4],25))
   Better_choices[,2] <- c(rep(Crops[2],25), rep(Crops[4],25),rep(Crops[1],25), rep(Crops[3],25))
   
   
   #Define migration shedule, migration every five rounds. Define when individuals from 1st group migrate, always switch with same "partner"
   Migration <- rep(c(0,0,0,0,1,0,0,0,0,2,0,0,0,0,3,0,0,0,0,4),5)
   
   #Initialize attractions, we assume there are no preexisting differences
   Homunculi$A1 <- 0
   Homunculi$A2 <- 0
   Homunculi$A3 <- 0
   Homunculi$A4 <- 0
   
   
   # Start simulation loop
   for (i in 1:Tmax) {
     
     #Loop over all individuals
     
     for (x in Homunculi$id){
       
       #Assign new vars for individual-specific parameter values, just to make code more readable
       
       sigma_Ind <- Homunculi$sigma[which(Homunculi$id==x)]
       phi_Ind <-   Homunculi$phi[which(Homunculi$id==x)]
       f_Ind <-     Homunculi$f[which(Homunculi$id==x)]
       b_Ind <-     Homunculi$b[which(Homunculi$id==x)]
       L_Ind <-     Homunculi$L[which(Homunculi$id==x)]
       
       
       #Record group_id in output matrix
       Result$group.id[which(Result$id==x & Result$trial==i)] <- Homunculi$Group[which(Homunculi$id==x)]
       
       #Record experience in present region
       #Individual has not migrated yet
       if (length(unique(na.omit(Result$group.id[which(Result$id==x)])))==1){
         Result$Experience[which(Result$id==x & Result$trial==i)] <- i
         
         #Individual has migrated already
       } else {
         Result$Experience[which(Result$id==x & Result$trial==i)] <- i- max(na.omit(Result$trial[which(Result$id==x & 
                                                                                                         Result$group.id != Result$group.id[which(Result$id==x & Result$trial==i)])]))
       }
       
       #Record experience in present environment so actual experience with one payoff structure
       #Environment has not changed yet
       if (i <= min(Times_Env_Change)){
         Result$TimeSinceChange[which(Result$id==x & Result$trial==i)] <- i
         
         #Environment has changed
       } else {
         Result$TimeSinceChange[which(Result$id==x & Result$trial==i)] <- i- max(Times_Env_Change[which(Times_Env_Change<i)])
       }
       
       if (Result$TimeSinceChange[which(Result$id==x & Result$trial==i)] >= Result$Experience[which(Result$id==x & Result$trial==i)]) {
         Result$Exp_Stable[which(Result$id==x & Result$trial==i)] <- Result$Experience[which(Result$id==x & Result$trial==i)]
       } else {
         Result$Exp_Stable[which(Result$id==x & Result$trial==i)] <- Result$TimeSinceChange[which(Result$id==x & Result$trial==i)]
       }
       
       
       
       #Find choices of group members in previous round
       
       #Get other current group members
       Group_members <- Homunculi$id[which(Homunculi$Group == Homunculi$Group[which(Homunculi$id == x)] & Homunculi$id != x)]
       
       #Indicates whether there is a newly arrived member, which cannot be copied
       
       New_Member <- Result$id[which(Result$id %in% Group_members & Result$trial==i & Result$Experience==1)]
       
       #Individuals that chose option 1
       n1 <- length(which(Result$Choice[which(Result$trial==i-1  & Result$id %in% Group_members & Result$id %not in% New_Member)] == 1))
       
       #Individuals that chose option 2
       n2 <- length(which(Result$Choice[which(Result$trial==i-1  & Result$id %in% Group_members & Result$id %not in% New_Member)] == 2))
       
       #Individuals that chose option 3
       n3 <- length(which(Result$Choice[which(Result$trial==i-1  & Result$id %in% Group_members & Result$id %not in% New_Member)] == 3))
       
       #Individuals that chose option 4
       n4 <- length(which(Result$Choice[which(Result$trial==i-1  & Result$id %in% Group_members & Result$id %not in% New_Member)] == 4))
       
       
       Result$n1[which(Result$id==x & Result$trial==i)] <- n1
       Result$n2[which(Result$id==x & Result$trial==i)] <- n2
       Result$n3[which(Result$id==x & Result$trial==i)] <- n3
       Result$n4[which(Result$id==x & Result$trial==i)] <- n4
       
       Result$Choice1[which(Result$id==x & Result$trial==i)] <- Result$Choice[which(Result$trial==i-1  & Result$id %in% Group_members & Result$id %not in% New_Member)][1]
       Result$Choice2[which(Result$id==x & Result$trial==i)] <- Result$Choice[which(Result$trial==i-1  & Result$id %in% Group_members & Result$id %not in% New_Member)][2]
       Result$Choice3[which(Result$id==x & Result$trial==i)] <- Result$Choice[which(Result$trial==i-1  & Result$id %in% Group_members & Result$id %not in% New_Member)][3]
       
       Result$Age1[which(Result$id==x & Result$trial==i)] <- Result$Experience[which(Result$trial==i-1  & Result$id %in% Group_members & Result$id %not in% New_Member)][1]
       Result$Age2[which(Result$id==x & Result$trial==i)] <- Result$Experience[which(Result$trial==i-1  & Result$id %in% Group_members & Result$id %not in% New_Member)][2]
       Result$Age3[which(Result$id==x & Result$trial==i)] <- Result$Experience[which(Result$trial==i-1  & Result$id %in% Group_members & Result$id %not in% New_Member)][3]
       
       
      {       
       #Experience of individuals that chose option 1
       if (n1 == 0){
         age1 <- 0
       } else {
         age1 <- sum(Result$Experience[which(Result$trial==i-1 & Result$id %in% Group_members & Result$Choice==1)])/n1
       }
       
       #Experience of individuals that chose option 2
       if (n2 == 0){
         age2 <- 0
       } else {
         age2 <- sum(Result$Experience[which(Result$trial==i-1 & Result$id %in% Group_members & Result$Choice==2)])/n2
       }    
       
       #Experience of individuals that chose option 3
       if (n3 == 0){
         age3 <- 0
       } else {
         age3 <- sum(Result$Experience[which(Result$trial==i-1 & Result$id %in% Group_members & Result$Choice==3)])/n3
       }  
       
       #Experience of individuals that chose option 4
       if (n4 == 0){
         age4 <- 0
       } else {
         age4 <- sum(Result$Experience[which(Result$trial==i-1 & Result$id %in% Group_members & Result$Choice==4)])/n4
       }  
       }
       
       #Generate choice based on previous attraction score and choice of other group members
       
       
       #Sum of attractions, frequencies and experiences
       Sum_Attr <- exp(L_Ind*Homunculi$A1[which(Homunculi$id==x)]) + 
         exp(L_Ind*Homunculi$A2[which(Homunculi$id==x)]) +
         exp(L_Ind*Homunculi$A3[which(Homunculi$id==x)]) +
         exp(L_Ind*Homunculi$A4[which(Homunculi$id==x)])
       
       
       
       #1st round
       #There are no choices to copy, so basic reinforcement equation
       
       if(i==1){
         Prob_1 <- exp(L_Ind*Homunculi$A1[which(Homunculi$id==x)]) / Sum_Attr
         
         Prob_2 <- exp(L_Ind*Homunculi$A2[which(Homunculi$id==x)]) / Sum_Attr
         
         Prob_3 <- exp(L_Ind*Homunculi$A3[which(Homunculi$id==x)]) / Sum_Attr
         
         Prob_4 <- exp(L_Ind*Homunculi$A4[which(Homunculi$id==x)]) / Sum_Attr
         
         #From 2nd round, combination of private and social information
         
       } else {
         
         Prob_1 <- (1-sigma_Ind)*(exp(L_Ind*Homunculi$A1[which(Homunculi$id==x)]))/ Sum_Attr +
           sigma_Ind* (((n1^f_Ind)*exp(b_Ind*age1))/((n1^f_Ind)*exp(b_Ind*age1)+(n2^f_Ind)*exp(b_Ind*age2)+(n3^f_Ind)*exp(b_Ind*age3)+(n4^f_Ind)*exp(b_Ind*age4)))
         
         Prob_2 <- (1-sigma_Ind)*(exp(L_Ind*Homunculi$A2[which(Homunculi$id==x)]))/ Sum_Attr +
           sigma_Ind* (((n2^f_Ind)*exp(b_Ind*age2))/((n1^f_Ind)*exp(b_Ind*age1)+(n2^f_Ind)*exp(b_Ind*age2)+(n3^f_Ind)*exp(b_Ind*age3)+(n4^f_Ind)*exp(b_Ind*age4)))
         
         Prob_3 <- (1-sigma_Ind)*(exp(L_Ind*Homunculi$A3[which(Homunculi$id==x)]))/ Sum_Attr +
           sigma_Ind* (((n3^f_Ind)*exp(b_Ind*age3))/((n1^f_Ind)*exp(b_Ind*age1)+(n2^f_Ind)*exp(b_Ind*age2)+(n3^f_Ind)*exp(b_Ind*age3)+(n4^f_Ind)*exp(b_Ind*age4)))
         
         Prob_4 <- (1-sigma_Ind)*(exp(L_Ind*Homunculi$A4[which(Homunculi$id==x)]))/ Sum_Attr +
           sigma_Ind* (((n4^f_Ind)*exp(b_Ind*age4))/((n1^f_Ind)*exp(b_Ind*age1)+(n2^f_Ind)*exp(b_Ind*age2)+(n3^f_Ind)*exp(b_Ind*age3)+(n4^f_Ind)*exp(b_Ind*age4)))
       }
       
       
       #Make choice proportional to attraction scores and social information
       Result$Choice[which(Result$id==x & Result$trial==i)] <- sample(c(1:4), size = 1, prob = c(Prob_1, Prob_2, Prob_3, Prob_4))
       
       
       #Generate a payoff based on choice
       
       #Individual chose better option
       if (Result$Choice[which(Result$id==x & Result$trial==i)] == Better_choices[i,Homunculi$Group[which(Homunculi$id == x)]]){
         Result$Payoff[which(Result$id==x & Result$trial==i)] <- round(rtruncnorm(1, a=0, b=+Inf, mean=Payoff_Better, sd=Result$SD_Payoff[which(Result$id==x & Result$trial==i)])) 
         Result$Correct[which(Result$id==x & Result$trial==i)] <- 1
         
         #Individual chose worse option
       } else {
         Result$Payoff[which(Result$id==x & Result$trial==i)] <- round(rtruncnorm(1, a=0, b=+Inf, mean=Payoff_Worse, sd=Result$SD_Payoff[which(Result$id==x & Result$trial==i)]))
         Result$Correct[which(Result$id==x & Result$trial==i)] <- 0
       }
       
       
       #Update attraction scores based on payoff
       
       #Individual chose 1
       if (Result$Choice[which(Result$id==x & Result$trial==i)]==1){
         Homunculi$A1[which(Homunculi$id==x)] <- (1-phi_Ind)*Homunculi$A1[which(Homunculi$id==x)] + phi_Ind*Result$Payoff[which(Result$id==x & Result$trial==i)]
         Homunculi$A2[which(Homunculi$id==x)] <- (1-phi_Ind)*Homunculi$A2[which(Homunculi$id==x)] + phi_Ind*0
         Homunculi$A3[which(Homunculi$id==x)] <- (1-phi_Ind)*Homunculi$A3[which(Homunculi$id==x)] + phi_Ind*0
         Homunculi$A4[which(Homunculi$id==x)] <- (1-phi_Ind)*Homunculi$A4[which(Homunculi$id==x)] + phi_Ind*0
       }
       
       #Individual chose 2
       if (Result$Choice[which(Result$id==x & Result$trial==i)]==2){
         Homunculi$A1[which(Homunculi$id==x)] <- (1-phi_Ind)*Homunculi$A1[which(Homunculi$id==x)] + phi_Ind*0
         Homunculi$A2[which(Homunculi$id==x)] <- (1-phi_Ind)*Homunculi$A2[which(Homunculi$id==x)] + phi_Ind*Result$Payoff[which(Result$id==x & Result$trial==i)]
         Homunculi$A3[which(Homunculi$id==x)] <- (1-phi_Ind)*Homunculi$A3[which(Homunculi$id==x)] + phi_Ind*0
         Homunculi$A4[which(Homunculi$id==x)] <- (1-phi_Ind)*Homunculi$A4[which(Homunculi$id==x)] + phi_Ind*0
       }
       
       #Individual chose 3
       if (Result$Choice[which(Result$id==x & Result$trial==i)]==3){
         Homunculi$A1[which(Homunculi$id==x)] <- (1-phi_Ind)*Homunculi$A1[which(Homunculi$id==x)] + phi_Ind*0
         Homunculi$A2[which(Homunculi$id==x)] <- (1-phi_Ind)*Homunculi$A2[which(Homunculi$id==x)] + phi_Ind*0
         Homunculi$A3[which(Homunculi$id==x)] <- (1-phi_Ind)*Homunculi$A3[which(Homunculi$id==x)] + phi_Ind*Result$Payoff[which(Result$id==x & Result$trial==i)]
         Homunculi$A4[which(Homunculi$id==x)] <- (1-phi_Ind)*Homunculi$A4[which(Homunculi$id==x)] + phi_Ind*0
       }
       
       #Individual chose 4
       if (Result$Choice[which(Result$id==x & Result$trial==i)]==4){
         Homunculi$A1[which(Homunculi$id==x)] <- (1-phi_Ind)*Homunculi$A1[which(Homunculi$id==x)] + phi_Ind*0
         Homunculi$A2[which(Homunculi$id==x)] <- (1-phi_Ind)*Homunculi$A2[which(Homunculi$id==x)] + phi_Ind*0
         Homunculi$A3[which(Homunculi$id==x)] <- (1-phi_Ind)*Homunculi$A3[which(Homunculi$id==x)] + phi_Ind*0
         Homunculi$A4[which(Homunculi$id==x)] <- (1-phi_Ind)*Homunculi$A4[which(Homunculi$id==x)] + phi_Ind*Result$Payoff[which(Result$id==x & Result$trial==i)]
       }
       
     }
     
     #Migration
     #Always the same migration partner
     
     if (Migration[i] == 1){
       Old1 <- Homunculi$Group[which(Homunculi$id==1)]
       Old5 <- Homunculi$Group[which(Homunculi$id==5)]
       Homunculi$Group[which(Homunculi$id==1)] <- Old5
       Homunculi$Group[which(Homunculi$id==5)] <- Old1
     }
     if (Migration[i] == 2){
       Old2 <- Homunculi$Group[which(Homunculi$id==2)]
       Old6 <- Homunculi$Group[which(Homunculi$id==6)]
       Homunculi$Group[which(Homunculi$id==2)] <- Old6
       Homunculi$Group[which(Homunculi$id==6)] <- Old2
     }
     if (Migration[i] == 3){
       Old3 <- Homunculi$Group[which(Homunculi$id==3)]
       Old7 <- Homunculi$Group[which(Homunculi$id==7)]
       Homunculi$Group[which(Homunculi$id==3)] <- Old7
       Homunculi$Group[which(Homunculi$id==7)] <- Old3
     }
     if (Migration[i] == 4){
       Old4 <- Homunculi$Group[which(Homunculi$id==4)]
       Old8 <- Homunculi$Group[which(Homunculi$id==8)]
       Homunculi$Group[which(Homunculi$id==4)] <- Old8
       Homunculi$Group[which(Homunculi$id==8)] <- Old4
     }
     
   }
   Result_Overall <- rbind(Result_Overall, Result)
   
   # return(Result_Overall)
 }
 return(Result_Overall)
}

 
 Result <- Sim_fct(Tmax=100,               #Number of trials
                   Group_size=4,
                   N_group=2,
                   N_sessions=10,
                   L=0.2,                   #Noise of coices, impact determined by size of payoff
                   Payoff_Better=15,
                   Payoff_Worse=10,
                   Hard_SD = 2,
                   Easy_SD = 0.5,                   
                   phi=0.4,                 #Updating parameter; weight of recent choices
                   sigma=0,               #Reliance on SL
                   f=2,                     #strength of conformity bias
                   b=0.3)                     #strength of age bias

 

 