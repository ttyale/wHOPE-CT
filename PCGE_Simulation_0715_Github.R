##------------------------------------------------------------------------------
##
## Title: Simulation Code for the Paper "Impact of complex, partially nested clustering in a three-arm individually randomized group treatment trial: A case study with the wHOPE trial"  
## Description: The code contains (1) a data generating process for the Primary Education Group Treatment(PC-GE) arm in the wHOPE trial and (2) Bayesian estimation for intraclass correlation coefficients(ICC).  
## Author: Tony Tong 
## Date: 2021-07-15
## Version 0.0.0
## Email: tgy.yale@gmail.com

##load required packages
require(MASS)
require(mvtnorm)
#####

## Simulation function with both the data generating process and the Bayesian estimation through a Gibbs sampler
## rho1 and rho2 are session and therapist level ICCs
## nsim is the number of simulations
## v and u are mapping matrices for patients to sessions and sessions to therapists 

sim<-function(rho1=0.01, rho2=0.10, v="v1",u="u1",nsim=1000){
 
  #store the output
  S<-5000  # chain length
  BETA<-SIGMA2<-TAU1<-TAU2<-ICC1<-ICC2<-matrix(NA,nsim,S) 
  
  n=275   #number of patients
  site<-5    #number of sites
  session<-200 #number of sessions
  therapist<-11   #number of therapists
  
  for (ss in 1:nsim){
   
    #generate pseudo data
    V<-matrix(0,ncol=session,nrow=n)
    U<-matrix(0,ncol=therapist,nrow=session)
    
    #assign weights
    if(v=="v1"){
      while (min(colSums(V))==0) {
        V<-matrix(0,ncol=session,nrow=n)
        for (j in 1:site){
          for (i in (55*(j-1)+1):(55*j)){
               V[i,sample(((j-1)*40+1):(j*40),5)]<-0.2
          }
        }
      }
    } 
      
    if (v=="v2"){
      while (min(colSums(V))==0) {
        V<-matrix(0,ncol=session,nrow=n)
        for (j in 1:site){
          lat_idx<-c(sample(((j-1)*40+16):(j*40),25), sample(((j-1)*40+16):(j*40),25))
          for (i in (55*(j-1)+1):(55*j)){
             V[i, sample(((j-1)*40+1):(j*40-25),4)]<-0.2
             V[i, lat_idx[i-55*(j-1)]]<-0.2
          }
        }
      }
    }  
      
    if (u=="u1"){
      U<-matrix(0,ncol=therapist,nrow=session)
      #pcgep<-c(1,1,1,1,1)
      #pcgeb<-c(1,1,1,2,1)
      #site 1-3
      for (j in 1:3){
        pidx<-sample(((j-1)*40+1):(j*40),36)
        sidx<-setdiff(((j-1)*40+1):(j*40),pidx)
        U[pidx,(j*2-1)]<-1
        U[sidx,(j*2)]<-1
      }
      
      #site 4 with 3 therapists
      pidx<-sample(((4-1)*40+1):(4*40),36)
      sidx<-setdiff(((4-1)*40+1):(4*40),pidx)
      sidx1<-sample(sidx,2)
      sidx2<-setdiff(sidx,sidx1)
      U[pidx,7]<-1
      U[sidx1,8]<-1
      U[sidx2,9]<-1
      
      #site 5
      pidx<-sample(((5-1)*40+1):(5*40),36)
      sidx<-setdiff(((5-1)*40+1):(5*40),pidx)
      U[pidx,10]<-1
      U[sidx,11]<-1
    }
    
    if (u=="u2"){
      U<-matrix(0,ncol=therapist,nrow=session)
      #pcgep<-c(1,1,1,1,1)
      #pcgeb<-c(1,1,1,2,1)
      #site 1-3
      for (j in 1:3){
        pidx<-sample(((j-1)*40+1):(j*40),24)
        sidx<-setdiff(((j-1)*40+1):(j*40),pidx)
        U[pidx,(j*2-1)]<-1
        U[sidx,(j*2)]<-1
      }
      
      #site 4 with 3 therapists
      pidx<-sample(((4-1)*40+1):(4*40),24)
      sidx<-setdiff(((4-1)*40+1):(4*40),pidx)
      sidx1<-sample(sidx,8)
      sidx2<-setdiff(sidx,sidx1)
      U[pidx,7]<-1
      U[sidx1,8]<-1
      U[sidx2,9]<-1
      
      #site 5
      pidx<-sample(((5-1)*40+1):(5*40),24)
      sidx<-setdiff(((5-1)*40+1):(5*40),pidx)
      U[pidx,10]<-1
      U[sidx,11]<-1
    }  
      
    
    #intercept and no site-effect for now
    X<-rep(1,n)
    beta<--1.5 #true beta
    
    #session effect is rho1= c(0.01,0.10); therapist effect is rho2=c(0.01,0.10), residual is sigma2=1--tau1-tau2
    tau1=rho1*5.9
    tau2=rho2*5.9
    sigma2=5.9-tau1-tau2
    
    #Generate D the session effect and C the therapist effect
    D<-rnorm(session, 0, sqrt(tau1))
    C<-rnorm(therapist,0,sqrt(tau2))
    
    #Genrate Y
    Y<-X*beta+V%*%((U%*%C+D))+rnorm(n,0,sqrt(sigma2))
    
    #####Fit the model with Gibbs sampler
  
    #non-informative priors
    beta0=rnorm(1,0,1000)
    sigma20=1000
    a0=0.001
    b0=0.001
    c0=0.001
    d0=0.001
    e0=0.001
    f0=0.001
    g0=0.001
    h0=0.001
    
    #Save output: tau1 is pie2; tau2 is phi2
    betap<-sigma2p<-tau1p<-tau2p<-rep(NA,S)  
    
    #initial values
    n=length(Y)
    t<-dim(U)[2]
    s<-dim(U)[1]
    p<-1
    
    beta<-0
    sigma2=1
    tau1=1
    tau2=1
    
    D<-rnorm(s,0,1)
    C<-rnorm(t,0,1)
    
    
    for (i in 1:S){
      #update beta
      VBeta<-solve(crossprod(X)/sigma2+1/sigma20)
      MBeta<-VBeta%*%(t(X)%*%(Y-V%*%D)/sigma2+beta0/sigma20) 
      beta<-c(rnorm(1,MBeta,VBeta))
      betap[i]<-beta
      
      #update D
      VD<-solve(crossprod(V)/sigma2+diag(s)/tau1)
      MD<-VD%*%(t(V)%*%(Y-X*beta)/sigma2+U%*%C/tau1)
      D<-c(rmvnorm(1,MD,VD))
      
      #update C
      VC<-solve(crossprod(U)/tau1+diag(t)/tau2)
      MC<-VC%*%(t(U)%*%D)/tau1
      C<-c(rmvnorm(1,MC,VC))
      
      #update tau1 (pie2)
      tau1<-1/rgamma(1,shape=e0+s/2,rate=f0+1/2*crossprod(D-U%*%C))
      tau1p[i]<-tau1
      
      #update tau2 (phi2)
      tau2<-1/rgamma(1,shape=g0+t/2,rate=h0+1/2*crossprod(C))
      tau2p[i]<-tau2
      
      #update sigma2
      sigma2<-1/rgamma(1,shape=a0+n/2,rate=b0+0.5*crossprod(Y-X*beta-V%*%D))
      sigma2p[i]<-sigma2
    
      #if (i%%100==0) print(i)
    }
    
    #store chain
    BETA[ss,]<-betap
    SIGMA2[ss,]<-sigma2p
    TAU1[ss,]<-tau1p
    TAU2[ss,]<-tau2p
    ICC1[ss,]<-tau1p/(tau1p+tau2p+sigma2p)
    ICC2[ss,]<-tau2p/(tau1p+tau2p+sigma2p)
    
    if (ss%%100==0) print(ss)
    
  }

  return(list(BETA=BETA,SIGMA2=SIGMA2, TAU1=TAU1, TAU2=TAU2, ICC1=ICC1, ICC2=ICC2))
  
}

##-----------------------------------------------------------------------------
##run the simulation
##set directory to save simulation output
#setwd("")

#V1U2
set.seed(20201202)
start_time <- Sys.time()
saveRDS(sim(rho1=0.01, rho2=0.01, v="v1", u="u2",nsim=1000),"PCGE_V1U2_ICC001_001.RDS")
end_time <- Sys.time()
end_time - start_time

set.seed(20201202)
start_time <- Sys.time()
saveRDS(sim(rho1=0.01, rho2=0.10, v="v1", u="u2",nsim=1000),"PCGE_V1U2_ICC001_010.RDS")
end_time <- Sys.time()
end_time - start_time

set.seed(20201202)
start_time <- Sys.time()
saveRDS(sim(rho1=0.10, rho2=0.01, v="v1", u="u2",nsim=1000),"PCGE_V1U2_ICC010_001.RDS")
end_time <- Sys.time()
end_time - start_time

set.seed(20201202)
start_time <- Sys.time()
saveRDS(sim(rho1=0.10, rho2=0.10, v="v1", u="u2",nsim=1000),"PCGE_V1U2_ICC010_010.RDS")
end_time <- Sys.time()
end_time - start_time


#V2U2
set.seed(20201202)
start_time <- Sys.time()
saveRDS(sim(rho1=0.01, rho2=0.01, v="v2", u="u2",nsim=1000),"PCGE_V2U2_ICC001_001.RDS")
end_time <- Sys.time()
end_time - start_time

set.seed(20201202)
start_time <- Sys.time()
saveRDS(sim(rho1=0.01, rho2=0.10, v="v2", u="u2",nsim=1000),"PCGE_V2U2_ICC001_010.RDS")
end_time <- Sys.time()
end_time - start_time

set.seed(20201202)
start_time <- Sys.time()
saveRDS(sim(rho1=0.10, rho2=0.01, v="v2", u="u2",nsim=1000),"PCGE_V2U2_ICC010_001.RDS")
end_time <- Sys.time()
end_time - start_time

set.seed(20201202)
start_time <- Sys.time()
saveRDS(sim(rho1=0.10, rho2=0.10, v="v2", u="u2",nsim=1000),"PCGE_V2U2_ICC010_010.RDS")
end_time <- Sys.time()
end_time - start_time


#V1U1
set.seed(20201202)
start_time <- Sys.time()
saveRDS(sim(rho1=0.01, rho2=0.01, v="v1", u="u1",nsim=1000),"PCGE_V1U1_ICC001_001.RDS")
end_time <- Sys.time()
end_time - start_time

set.seed(20201202)
start_time <- Sys.time()
saveRDS(sim(rho1=0.01, rho2=0.10, v="v1", u="u1",nsim=1000),"PCGE_V1U1_ICC001_010.RDS")
end_time <- Sys.time()
end_time - start_time

set.seed(20201202)
start_time <- Sys.time()
saveRDS(sim(rho1=0.10, rho2=0.01, v="v1", u="u1",nsim=1000),"PCGE_V1U1_ICC010_001.RDS")
end_time <- Sys.time()
end_time - start_time

set.seed(20201202)
start_time <- Sys.time()
saveRDS(sim(rho1=0.10, rho2=0.10, v="v1", u="u1",nsim=1000),"PCGE_V1U1_ICC010_010.RDS")
end_time <- Sys.time()
end_time - start_time


#V2U1
set.seed(20201202)
start_time <- Sys.time()
saveRDS(sim(rho1=0.01, rho2=0.01, v="v2", u="u1",nsim=1000),"PCGE_V2U1_ICC001_001.RDS")
end_time <- Sys.time()
end_time - start_time

set.seed(20201202)
start_time <- Sys.time()
saveRDS(sim(rho1=0.01, rho2=0.10, v="v2", u="u1",nsim=1000),"PCGE_V2U1_ICC001_010.RDS")
end_time <- Sys.time()
end_time - start_time

set.seed(20201202)
start_time <- Sys.time()
saveRDS(sim(rho1=0.10, rho2=0.01, v="v2", u="u1",nsim=1000),"PCGE_V2U1_ICC010_001.RDS")
end_time <- Sys.time()
end_time - start_time

set.seed(20201202)
start_time <- Sys.time()
saveRDS(sim(rho1=0.10, rho2=0.10, v="v2", u="u1",nsim=1000),"PCGE_V2U1_ICC010_010.RDS")
end_time <- Sys.time()
end_time - start_time


##------------------------------------------------------------------------------
##check posterior
#load the data and name it x1
x1-readRDS("")

burnin=1500

#posterior median and mean
median(x1$BETA[1,-c(1:burnin)])
mean(x1$BETA[1,-c(1:burnin)])

median(x1$ICC1[1,-c(1:burnin)])
mean(x1$ICC1[1,-c(1:burnin)])

median(x1$ICC2[1,-c(1:burnin)])
mean(x1$ICC2[1,-c(1:burnin)])









