##------------------------------------------------------------------------------
##
## Title: Simulation Code for the Paper "Impact of complex, partially nested clustering in a three-arm individually randomized group treatment trial: A case study with the wHOPE trial"  
## Description: The code contains (1) a data generating process for the Whole Health Team(WHT) arm in the wHOPE trial and (2) Bayesian estimation for intraclass correlation coefficients(ICC).  
## Author: Tony Tong 
## Date: 2021-07-15
## Version 1.1.0
## Email: tgy.yale@gmail.com

##load required packages
require(MASS)
require(mvtnorm)


##Simulation function with both the data generating process and the Bayesian estimation through a Gibbs sampler.
##condition option is the mapping matrices of different conditions, avaliable from "w1" to "w6". 
##nsim is the number of simulations

sim<-function(rho=0.01,condition="w6",nsim=1000){

  S<-5000 #chain length
  BETA<-SIGMA2<-TAU2<-ICC<-matrix(NA,nsim,S) #create empty output matrices to store key parameters 
  
  
  n=275   #number of patients
  site<-5    #number of sites
  therapist1<-10   #number of WH coaches
  therapist2<-8    #number of PCP
  therapist3<-9    #number of CIH
  therapist4<-4    #number of therapists
  p1=0.75          #probability of patient encountering partial clustering
  
  for (ss in 1:nsim){
   
    ## Generate mapping matrices under six conditions from w6-w1. Descriptions can be found in the paper draft.
    #w6: all therapists with max encounters
    if (condition=="w6"){ 
      
      #WHC weights: For situation with two therapists, primary has weight of pw=5/6 for WHC and 3/4 for the other roles.  
      W1<-matrix(0,ncol=therapist1,nrow=n)
      #assign weights
      for (i in 1:55){
        if (rbinom(1,1,p1)==1) W1[i,1]=12 else W1[i,1:2]=c(5/6*12,(1-5/6)*12)
      }
      for (i in 56:110){
        W1[i,3]=12
      }
      for (i in 111:165){
        if (rbinom(1,1,p1)==1) W1[i,sample(4:6,1)]=12 else W1[i,sample(4:6,2)]=c(5/6*12,(1-5/6)*12)
      }
      for (i in 166:220){
        if (rbinom(1,1,p1)==1) W1[i,sample(7:8,1)]=12 else W1[i,sample(7:8,2)]=c(5/6*12,(1-5/6)*12)
      }
      for (i in 221:275){
        if (rbinom(1,1,p1)==1) W1[i,9]=12 else W1[i,9:10]=c(5/6*12,(1-5/6)*12)
      }
     
      #PCP weights: 3/4 chance of encounter no backup PCP if available, if encountering, only one is allowed.
      W2<-matrix(0,ncol=therapist2,nrow=n)
      #assign weights
      for (i in 1:55){
        if (rbinom(1,1,p1)==1) W2[i,1]=4 else W2[i,1:2]=c(3,1)
      }
      for (i in 56:110){
        if (rbinom(1,1,p1)==1) W2[i,3]=4 else W2[i,3:4]=c(3,1)
      }
      for (i in 111:165){
        W2[i,5]=4
      }
      for (i in 166:220){
        if (rbinom(1,1,p1)==1) W2[i,sample(6:7,1)]=4 else W2[i,sample(6:7,2)]=c(3,1)
      }
      for (i in 221:275){
        W2[i,8]=4
      }
      
      #CIH weights: 3/4 chance of encounter no backup PCP if available, if encountering, only one is allowed.
      W3<-matrix(0,ncol=therapist3,nrow=n)
      #assign weights
      for (i in 1:55){
        if (rbinom(1,1,p1)==1) W3[i,1]=4 else W3[i,1:2]=c(3,1)
      }
      for (i in 56:110){
        if (rbinom(1,1,p1)==1) W3[i,sample(3:4,1)]=4 else W3[i,c(sample(3:4,1),sample(5:6,1))]=c(3,1)
      }
      for (i in 111:165){
        W3[i,7]=4
      }
      for (i in 166:220){
        W3[i,8]=4
      }
      for (i in 221:275){
        W3[i,9]=4
      }
      
      #MHP weights: only zero or one each site
      W4<-matrix(0,ncol=therapist4,nrow=n)
      #assign weights
      for (i in 1:55){
        W4[i,1]=4
      }
      for (i in 111:165){
        W4[i,2]=4
      }
      for (i in 166:220){
        W4[i,3]=4
      }
      for (i in 221:275){
        W4[i,4]=4
      }
      
      W<-as.matrix(cbind(W1,W2,W3,W4))
      W<-t(apply(W,1, function(x){x/sum(x)})) #scale to row sum of 1
  
      therapist<-therapist1+therapist2+therapist3+therapist4 #total number of therapists
    }  
    
    
    #### min backup condition
    #w5: all therapists with min encounters
    if (condition=="w5"){
      #WHC weights: Each patient encounters 2 therapists with  p1= 3/4. For two therapists, primary has weight of pw=8/9.  
      W1<-matrix(0,ncol=therapist1,nrow=n)
      #assign weights
      for (i in 1:55){
        if (rbinom(1,1,p1)==1) W1[i,1]=9 else W1[i,1:2]=c(8/9*9,(1-8/9)*9)
      }
      for (i in 56:110){
        W1[i,3]=9
      }
      for (i in 111:165){
        if (rbinom(1,1,p1)==1) W1[i,sample(4:6,1)]=9 else W1[i,sample(4:6,2)]=c(8/9*9,(1-8/9)*9)
      }
      for (i in 166:220){
        if (rbinom(1,1,p1)==1) W1[i,sample(7:8,1)]=9 else W1[i,sample(7:8,2)]=c(8/9*9,(1-8/9)*9)
      }
      for (i in 221:275){
        if (rbinom(1,1,p1)==1) W1[i,9]=9 else W1[i,9:10]=c(8/9*9,(1-8/9)*9)
      }
      
      # still 4 encounters with PCP
      W2<-matrix(0,ncol=therapist2,nrow=n)
      #assign weights
      for (i in 1:55){
        if (rbinom(1,1,p1)==1) W2[i,1]=4 else W2[i,1:2]=c(3,1)
      }
      for (i in 56:110){
        if (rbinom(1,1,p1)==1) W2[i,3]=4 else W2[i,3:4]=c(3,1)
      }
      for (i in 111:165){
        W2[i,5]=4
      }
      for (i in 166:220){
        if (rbinom(1,1,p1)==1) W2[i,sample(6:7,1)]=4 else W2[i,sample(6:7,2)]=c(3,1)
      }
      for (i in 221:275){
        W2[i,8]=4
      }
      
      #CIH weights: one encounter for CIH over the whole time; 3/4 chance of encounter no backup PCP if available, if encountering, only one is allowed.
      W3<-matrix(0,ncol=therapist3,nrow=n)
      #assign weights
      for (i in 1:55){
        if (rbinom(1,1,p1)==1) W3[i,1]=1 else W3[i,1:2]=c(0,1)
      }
      for (i in 56:110){
        if (rbinom(1,1,p1)==1) W3[i,sample(3:4,1)]=1 else W3[i,c(sample(3:4,1),sample(5:6,1))]=c(0,1)
      }
      for (i in 111:165){
        W3[i,7]=1
      }
      for (i in 166:220){
        W3[i,8]=1
      }
      for (i in 221:275){
        W3[i,9]=1
      }
      
      #MHP weights: only zero or one each site
      W4<-matrix(0,ncol=therapist4,nrow=n)
      #assign weights
      for (i in 1:55){
        W4[i,1]=1
      }
      for (i in 111:165){
        W4[i,2]=1
      }
      for (i in 166:220){
        W4[i,3]=1
      }
      for (i in 221:275){
        W4[i,4]=1
      }
      
      W<-as.matrix(cbind(W1,W2,W3,W4))
      W<-t(apply(W,1, function(x){x/sum(x)})) #scale to row sum of 1
      
      therapist<-therapist1+therapist2+therapist3+therapist4
    }
    
    #w4: only primary therapists with max encounters
    ####### max primary condition w4
    #the chance of encountering co-primary is 5/6 out of 12 for WHT arm
    if (condition=="w4"){
      #WHC weights: Each patient encounters only primary
      W1<-matrix(0,ncol=8,nrow=n)
      #assign weights
      for (i in 1:55){
        W1[i,1]=12 
      }
      for (i in 56:110){
        W1[i,2]=12
      }
      for (i in 111:165){
        if (rbinom(1,1,p1)==1) W1[i,sample(3:5,1)]=12 else W1[i,sample(3:5,2)]=c(5/6*12,(1-5/6)*12)
      }
      for (i in 166:220){
        if (rbinom(1,1,p1)==1) W1[i,sample(6:7,1)]=12 else W1[i,sample(6:7,2)]=c(5/6*12,(1-5/6)*12)
      }
      for (i in 221:275){
        W1[i,2]=12
      }
      
      #PCP weights: 3/4 chance of encounter no backup PCP if available, if encountering, only one is allowed.
      W2<-matrix(0,ncol=6,nrow=n)
      #assign weights
      for (i in 1:55){
         W2[i,1]=4
      }
      for (i in 56:110){
         W2[i,2]=4
      }
      for (i in 111:165){
        W2[i,3]=4
      }
      for (i in 166:220){
        if (rbinom(1,1,p1)==1) W2[i,sample(4:5,1)]=4 else W2[i,sample(4:5,2)]=c(3,1)
      }
      for (i in 221:275){
        W2[i,6]=4
      }
      
      #CIH weights: 3/4 chance of encounter no backup PCP if available, if encountering, only one is allowed.
      W3<-matrix(0,ncol=6,nrow=n)
      #assign weights
      for (i in 1:55){
        W3[i,1]=4
      }
      for (i in 56:110){
        if (rbinom(1,1,p1)==1) W3[i,sample(2:3,1)]=4 else W3[i,sample(2:3,2)]=c(3,1)
      }
      for (i in 111:165){
        W3[i,4]=4
      }
      for (i in 166:220){
        W3[i,5]=4
      }
      for (i in 221:275){
        W3[i,6]=4
      }
      
      #MHP weights: only zero or one each site
      W4<-matrix(0,ncol=therapist4,nrow=n)
      #assign weights
      for (i in 1:55){
        W4[i,1]=4
      }
      for (i in 111:165){
        W4[i,2]=4
      }
      for (i in 166:220){
        W4[i,3]=4
      }
      for (i in 221:275){
        W4[i,4]=4
      }
      
      W<-as.matrix(cbind(W1,W2,W3,W4))
      W<-t(apply(W,1, function(x){x/sum(x)})) #scale to row sum of 1
      
      therapist<-24
    }
    
    #w3: only primary therapists with min encounters
    if (condition=="w3"){
      #WHC weights: Each patient encounters 2 therapists with  p1= 3/4. For two therapists, primary has weight of pw=8/9.  
      W1<-matrix(0,ncol=8,nrow=n)
      #assign weights
      for (i in 1:55){
        W1[i,1]=9
      }
      for (i in 56:110){
        W1[i,2]=9
      }
      for (i in 111:165){
        if (rbinom(1,1,p1)==1) W1[i,sample(3:5,1)]=9 else W1[i,sample(3:5,2)]=c(8/9*9,(1-8/9)*9)
      }
      for (i in 166:220){
        if (rbinom(1,1,p1)==1) W1[i,sample(6:7,1)]=9 else W1[i,sample(6:7,2)]=c(8/9*9,(1-8/9)*9)
      }
      for (i in 221:275){
        W1[i,8]=9
      }
      
      # still 4 encounters with PCP
      W2<-matrix(0,ncol=6,nrow=n)
      #assign weights
      for (i in 1:55){
        W2[i,1]=4
      }
      for (i in 56:110){
        W2[i,2]=4
      }
      for (i in 111:165){
        W2[i,3]=4
      }
      for (i in 166:220){
        if (rbinom(1,1,p1)==1) W2[i,sample(4:5,1)]=4 else W2[i,sample(4:5,2)]=c(3,1)
      }
      for (i in 221:275){
        W2[i,6]=4
      }
      
      #CIH weights: one encounter for CIH over the whole time
      W3<-matrix(0,ncol=6,nrow=n)
      #assign weights
      for (i in 1:55){
        W3[i,1]=1
      }
      for (i in 56:110){
        W3[i,sample(2:3,2)]=c(0,1)
      }
      for (i in 111:165){
        W3[i,4]=1
      }
      for (i in 166:220){
        W3[i,5]=1
      }
      for (i in 221:275){
        W3[i,6]=1
      }
      
      #MHP weights: only zero or one each site
      W4<-matrix(0,ncol=4,nrow=n)
      #assign weights
      for (i in 1:55){
        W4[i,1]=1
      }
      for (i in 111:165){
        W4[i,2]=1
      }
      for (i in 166:220){
        W4[i,3]=1
      }
      for (i in 221:275){
        W4[i,4]=1
      }
      
      W<-as.matrix(cbind(W1,W2,W3,W4))
      W<-t(apply(W,1, function(x){x/sum(x)})) #scale to row sum of 1
      
      therapist<-24
    }
    
    #w2:All Whole Health Coaches
    if (condition=="w2"){
      therapist=10
      W<-matrix(0,ncol=therapist,nrow=n)
      pw=0.875
      #assign weights
      for (i in 1:55){
        if (rbinom(1,1,p1)==1) W[i,1]=1 else W[i,1:2]=c(pw,1-pw)
      }
      for (i in 56:110){
        W[i,3]=1
      }
      for (i in 111:165){
        if (rbinom(1,1,p1)==1) W[i,sample(4:6,1)]=1 else W[i,sample(4:6,2)]=c(pw,1-pw)
      }
      for (i in 166:220){
        if (rbinom(1,1,p1)==1) W[i,sample(7:8,1)]=1 else W[i,sample(7:8,2)]=c(pw,1-pw)
      }
      for (i in 221:275){
        if (rbinom(1,1,p1)==1) W[i,9]=1 else W[i,9:10]=c(pw,1-pw)
      }
    }
    
    #w1: Only primary Whole Health Coaches
    if (condition=="w1"){
      pw=0.875
      therapist=8
      W<-matrix(0,ncol=therapist,nrow=n)
      #assign weights
      for (i in 1:55){
        W[i,1]=1
      }
      for (i in 56:110){
        W[i,2]=1
      }
      for (i in 111:165){
        if (rbinom(1,1,p1)==1) W[i,sample(3:5,1)]=1 else W[i,sample(3:5,2)]=c(pw,1-pw)
      }
      for (i in 166:220){
        if (rbinom(1,1,p1)==1) W[i,sample(6:7,1)]=1 else W[i,sample(6:7,2)]=c(pw,1-pw)
      }
      for (i in 221:275){
        W[i,8]=1
      }
    }
    
    # design matrix with only intercept and no site-effect for now
    X<-rep(1,n)
  
    #therapist effect's variance is rho=c(0.01,0.02,0.05, 0.10,0.15) multiply by outcome variance, which is 5.9.
    tau2=rho*5.9
    sigma2=5.9-tau2
    
    #Generate B, the therapist effect
    B<-rnorm(therapist,0,sqrt(tau2))
    
    #Generate outcome Y
    beta<--1.5 #true beta
    Y<-X*beta+W%*%B+rnorm(n,0,sqrt(sigma2))

        
    ##Fit the model with Gibbs sampler
    #non-informative priors for all parameters
    beta0=rnorm(1,0,1000)
    sigma20=1000
    a0=0.001
    b0=0.001
    c0=0.001
    d0=0.001
    
    #Save output
    betap<-sigma2p<-tau2p<-rep(NA,S)  
    
    #initial values
    n=length(Y)
    t<-dim(W)[2]
    p<-1
    
    beta<-0
    sigma2=1
    tau2=1
    
    B<-rnorm(t,0,1)
    
    for (i in 1:S){
      #update beta
      VBeta<-solve(crossprod(X)/sigma2+1/sigma20)
      MBeta<-VBeta%*%(t(X)%*%(Y-W%*%B)/sigma2+beta0/sigma20) 
      beta<-c(rnorm(1,MBeta,VBeta))
      betap[i]<-beta
      
      #update B
      VB<-solve(crossprod(W)/sigma2+diag(t)/tau2)
      MB<-VB%*%(t(W)%*%(Y-X*beta))/sigma2
      B<-c(rmvnorm(1,MB,VB))
      
      #update tau2
      tau2<-1/rgamma(1,shape=c0+t/2,rate=d0+1/2*crossprod(B))
      tau2p[i]<-tau2
      
      #update sigma2
      sigma2<-1/rgamma(1,shape=a0+n/2,rate=b0+0.5*crossprod(Y-X*beta-W%*%B))
      sigma2p[i]<-sigma2
      # print(i)
    }
    
    #store
    BETA[ss,]<-betap
    SIGMA2[ss,]<-sigma2p
    TAU2[ss,]<-tau2p
    ICC[ss,]<-tau2p/(tau2p+sigma2p)
    
    if (ss%%200==0) print(ss) #print after every 200 iterations
    
  }

  return(list(BETA=BETA,SIGMA2=SIGMA2,TAU2=TAU2, ICC=ICC))
  
}


#setwd("") #set working directory

#run simulation for each mapping matrix and ICC.
start_time <- Sys.time()
saveRDS(sim(rho=0.01,condition="w6",nsim=1000),"W6_WHT_ICC001.RDS")
saveRDS(sim(rho=0.02,condition="w6",nsim=1000),"W6_WHT_ICC002.RDS")
saveRDS(sim(rho=0.05,condition="w6",nsim=1000),"W6_WHT_ICC005.RDS")
saveRDS(sim(rho=0.10,condition="w6",nsim=1000),"W6_WHT_ICC010.RDS")
saveRDS(sim(rho=0.15,condition="w6",nsim=1000),"W6_WHT_ICC015.RDS")

end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
saveRDS(sim(rho=0.01,condition="w5",nsim=1000),"W5_WHT_ICC001.RDS")
saveRDS(sim(rho=0.02,condition="w5",nsim=1000),"W5_WHT_ICC002.RDS")
saveRDS(sim(rho=0.05,condition="w5",nsim=1000),"W5_WHT_ICC005.RDS")
saveRDS(sim(rho=0.10,condition="w5",nsim=1000),"W5_WHT_ICC010.RDS")
saveRDS(sim(rho=0.15,condition="w5",nsim=1000),"W5_WHT_ICC015.RDS")

end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
saveRDS(sim(rho=0.01,condition="w4",nsim=1000),"W4_WHT_ICC001.RDS")
saveRDS(sim(rho=0.02,condition="w4",nsim=1000),"W4_WHT_ICC002.RDS")
saveRDS(sim(rho=0.05,condition="w4",nsim=1000),"W4_WHT_ICC005.RDS")
saveRDS(sim(rho=0.10,condition="w4",nsim=1000),"W4_WHT_ICC010.RDS")
saveRDS(sim(rho=0.15,condition="w4",nsim=1000),"W4_WHT_ICC015.RDS")

end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
saveRDS(sim(rho=0.01,condition="w3",nsim=1000),"W3_WHT_ICC001.RDS")
saveRDS(sim(rho=0.02,condition="w3",nsim=1000),"W3_WHT_ICC002.RDS")
saveRDS(sim(rho=0.05,condition="w3",nsim=1000),"W3_WHT_ICC005.RDS")
saveRDS(sim(rho=0.10,condition="w3",nsim=1000),"W3_WHT_ICC010.RDS")
saveRDS(sim(rho=0.15,condition="w3",nsim=1000),"W3_WHT_ICC015.RDS")

end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
saveRDS(sim(rho=0.01,condition="w2",nsim=1000),"W2_WHT_ICC001.RDS")
saveRDS(sim(rho=0.02,condition="w2",nsim=1000),"W2_WHT_ICC002.RDS")
saveRDS(sim(rho=0.05,condition="w2",nsim=1000),"W2_WHT_ICC005.RDS")
saveRDS(sim(rho=0.10,condition="w2",nsim=1000),"W2_WHT_ICC010.RDS")
saveRDS(sim(rho=0.15,condition="w2",nsim=1000),"W2_WHT_ICC015.RDS")

end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
saveRDS(sim(rho=0.01,condition="w1",nsim=1000),"W1_WHT_ICC001.RDS")
saveRDS(sim(rho=0.02,condition="w1",nsim=1000),"W1_WHT_ICC002.RDS")
saveRDS(sim(rho=0.05,condition="w1",nsim=1000),"W1_WHT_ICC005.RDS")
saveRDS(sim(rho=0.10,condition="w1",nsim=1000),"W1_WHT_ICC010.RDS")
saveRDS(sim(rho=0.15,condition="w1",nsim=1000),"W1_WHT_ICC015.RDS")

end_time <- Sys.time()
end_time - start_time


##------------------------------------------------------------------------------
##Result summary
burnin=2000

#Scenario 1
x1<-readRDS("W1_WHT_ICC001.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W1_WHT_ICC002.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W1_WHT_ICC005.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W1_WHT_ICC010.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W1_WHT_ICC015.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))


#Scenario 2
x1<-readRDS("W2_WHT_ICC001.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W2_WHT_ICC002.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W2_WHT_ICC005.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W2_WHT_ICC010.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W2_WHT_ICC015.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))


#Scenario 3
x1<-readRDS("W3_WHT_ICC001.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W3_WHT_ICC002.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W3_WHT_ICC005.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W3_WHT_ICC010.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W3_WHT_ICC015.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))


#Scenario 4
x1<-readRDS("W4_WHT_ICC001.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W4_WHT_ICC002.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W4_WHT_ICC005.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W4_WHT_ICC010.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W4_WHT_ICC015.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))


#Scenario 5
x1<-readRDS("W5_WHT_ICC001.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W5_WHT_ICC002.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W5_WHT_ICC005.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W5_WHT_ICC010.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W5_WHT_ICC015.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

#Scenario 6
x1<-readRDS("W6_WHT_ICC001.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W6_WHT_ICC002.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W6_WHT_ICC005.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W6_WHT_ICC010.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

x1<-readRDS("W6_WHT_ICC015.RDS")
mean(rowMeans(x1$ICC[,-(1:burnin)]))
mean(apply(x1$ICC[,-(1:burnin)],1,median))

