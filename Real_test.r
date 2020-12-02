rm(list=ls())
library(MASS)
library(mvtnorm)
library(matrixcalc)
library(rlecuyer)
library(edgeR)
library(parallel)
# source("https://bioconductor.org/biocLite.R")
# biocLite("edgeR")
library("edgeR")
# library(CDROM)
# ,lib.loc="/storage/home/bzk18/R/x86_64-redhat-linux-gnu-library/3.4"
#Expected gene expression of Parent
ExpP <- function(ThetaP, ThetaA, Alpha, TimePC){
  return ( ThetaP + (ThetaA - ThetaP) * exp(- Alpha * TimePC) )
}

#Expected gene expression of Child
ExpC <- function(ThetaC, ThetaA, Alpha, TimePC){
  return (ThetaC + (ThetaA - ThetaC) * exp(- Alpha * TimePC) )
}

#Expected gene expression of Ancestor
ExpA <- function(ThetaA){
  return (ThetaA)
}

#Funciton to output log-likelihood for neochild model given parameters
logll.neochild <- function(x,other.data){
  if(all(is.finite(x)) & all(x>=0)){ 
    ObsX <- c(other.data$Obs$X_carcass, other.data$Obs$X_fHead, other.data$Obs$X_ovary, other.data$Obs$X_mHead, other.data$Obs$X_testis, other.data$Obs$X_accg)
    T.PC <- other.data$T.PC
    T.PCA <- other.data$T.PCA
    ExpX <- c()
    ExpX["carcass_P"] <- ExpP(ThetaP=x[1],ThetaA=x[1],Alpha=x[3],TimePC=T.PC)
    ExpX["carcass_C"] <- ExpC(ThetaC=x[2],ThetaA=x[1],Alpha=x[3],TimePC=T.PC)
    ExpX["carcass_A"] <- ExpA(ThetaA=x[1])
    
    ExpX["fHead_P"] <- ExpP(ThetaP=x[1],ThetaA=x[1],Alpha=x[4],TimePC=T.PC)
    ExpX["fHead_C"] <- ExpC(ThetaC=x[2],ThetaA=x[1],Alpha=x[4],TimePC=T.PC)
    ExpX["fHead_A"] <- ExpA(ThetaA=x[1])
    
    ExpX["ovary_P"] <- ExpP(ThetaP=x[1],ThetaA=x[1],Alpha=x[5],TimePC=T.PC)
    ExpX["ovary_C"] <- ExpC(ThetaC=x[2],ThetaA=x[1],Alpha=x[5],TimePC=T.PC)
    ExpX["ovary_A"] <- ExpA(ThetaA=x[1])
    
    ExpX["mHead_P"] <- ExpP(ThetaP=x[1],ThetaA=x[1],Alpha=x[6],TimePC=T.PC)
    ExpX["mHead_C"] <- ExpC(ThetaC=x[2],ThetaA=x[1],Alpha=x[6],TimePC=T.PC)
    ExpX["mHead_A"] <- ExpA(ThetaA=x[1])
    
    ExpX["testis_P"] <- ExpP(ThetaP=x[1],ThetaA=x[1],Alpha=x[7],TimePC=T.PC)
    ExpX["testis_C"] <- ExpC(ThetaC=x[2],ThetaA=x[1],Alpha=x[7],TimePC=T.PC)
    ExpX["testis_A"] <- ExpA(ThetaA=x[1])
    
    ExpX["accg_P"] <- ExpP(ThetaP=x[1],ThetaA=x[1],Alpha=x[8],TimePC=T.PC)
    ExpX["accg_C"] <- ExpC(ThetaC=x[2],ThetaA=x[1],Alpha=x[8],TimePC=T.PC)
    ExpX["accg_A"] <- ExpA(ThetaA=x[1])
    
    var <- sapply(3:8, function(i) x[9]/(2*x[i]))
    cov_t_t1 <- sapply(3:8, function(i) exp(- 2 * x[i] * T.PC) * var[i-2])
    cov_t_t2 <- sapply(3:8, function(i) exp(- 2 * x[i] * T.PC) * var[i-2])
    Var.mat <- diag(as.vector(sapply(var, function(x) rep(x,3))),nrow=length(ExpX),ncol=length(ExpX))
    j <- 1
    for (entry in seq(1,length(ExpX),by=3)){
      Var.mat[entry,entry+1] <- cov_t_t1[j]
      Var.mat[entry+1,entry] <- cov_t_t1[j]
      Var.mat[entry,entry+2] <- cov_t_t2[j]
      Var.mat[entry+2,entry] <- cov_t_t2[j]
      Var.mat[entry+1,entry+2] <- cov_t_t2[j]
      Var.mat[entry+2,entry+1] <- cov_t_t2[j]
      j <- j+ 1
    }
    
    if(any(!is.numeric(Var.mat))){
      return(NA)}
    else if(any(!is.finite(Var.mat))){
      return (NA)
    }
    else if( !is.positive.semi.definite(Var.mat)){return (NA)}
    else{
      #print(Var.mat)
      retval <- dmvnorm(x=ObsX,mean=ExpX, sigma=Var.mat,log=1)    
      if (!is.finite(retval))
      {return (NA)
      }else { return (retval) }
      
    }
  }
  else return(NA)
}



#Function to return negative log of MLE for neoparent model
logll.neoparent <- function(x,other.data){
  if(all(is.finite(x)) & all(x>=0)){ 
    ObsX <- c(other.data$Obs$X_carcass, other.data$Obs$X_fHead, other.data$Obs$X_ovary, other.data$Obs$X_mHead, other.data$Obs$X_testis, other.data$Obs$X_accg)
    T.PC <- other.data$T.PC
    T.PCA <- other.data$T.PCA
    ExpX <- c()
    ExpX["carcass_P"] <- ExpP(ThetaP=x[1],ThetaA=x[2],Alpha=x[3],TimePC=T.PC)
    ExpX["carcass_C"] <- ExpC(ThetaC=x[2],ThetaA=x[2],Alpha=x[3],TimePC=T.PC)
    ExpX["carcass_A"] <- ExpA(ThetaA=x[2])
    
    ExpX["fHead_P"] <- ExpP(ThetaP=x[1],ThetaA=x[2],Alpha=x[4],TimePC=T.PC)
    ExpX["fHead_C"] <- ExpC(ThetaC=x[2],ThetaA=x[2],Alpha=x[4],TimePC=T.PC)
    ExpX["fHead_A"] <- ExpA(ThetaA=x[2])
    
    ExpX["ovary_P"] <- ExpP(ThetaP=x[1],ThetaA=x[2],Alpha=x[5],TimePC=T.PC)
    ExpX["ovary_C"] <- ExpC(ThetaC=x[2],ThetaA=x[2],Alpha=x[5],TimePC=T.PC)
    ExpX["ovary_A"] <- ExpA(ThetaA=x[2])
    
    ExpX["mHead_P"] <- ExpP(ThetaP=x[1],ThetaA=x[2],Alpha=x[6],TimePC=T.PC)
    ExpX["mHead_C"] <- ExpC(ThetaC=x[2],ThetaA=x[2],Alpha=x[6],TimePC=T.PC)
    ExpX["mHead_A"] <- ExpA(ThetaA=x[2])
    
    ExpX["testis_P"] <- ExpP(ThetaP=x[1],ThetaA=x[2],Alpha=x[7],TimePC=T.PC)
    ExpX["testis_C"] <- ExpC(ThetaC=x[2],ThetaA=x[2],Alpha=x[7],TimePC=T.PC)
    ExpX["testis_A"] <- ExpA(ThetaA=x[2])
    
    ExpX["accg_P"] <- ExpP(ThetaP=x[1],ThetaA=x[2],Alpha=x[8],TimePC=T.PC)
    ExpX["accg_C"] <- ExpC(ThetaC=x[2],ThetaA=x[2],Alpha=x[8],TimePC=T.PC)
    ExpX["accg_A"] <- ExpA(ThetaA=x[2])
    
    var <- sapply(3:8, function(i) x[9]/(2*x[i]))
    cov_t_t1 <- sapply(3:8, function(i) exp(- 2 * x[i] * T.PC) * var[i-2])
    cov_t_t2 <- sapply(3:8, function(i) exp(- 2 * x[i] * T.PC) * var[i-2])
    Var.mat <- diag(as.vector(sapply(var, function(x) rep(x,3))),nrow=length(ExpX),ncol=length(ExpX))
    j <- 1
    for (entry in seq(1,length(ExpX),by=3)){
      Var.mat[entry,entry+1] <- cov_t_t1[j]
      Var.mat[entry+1,entry] <- cov_t_t1[j]
      Var.mat[entry,entry+2] <- cov_t_t2[j]
      Var.mat[entry+2,entry] <- cov_t_t2[j]
      Var.mat[entry+1,entry+2] <- cov_t_t2[j]
      Var.mat[entry+2,entry+1] <- cov_t_t2[j]
      j <- j+ 1
    }
    if(any(!is.numeric(Var.mat))){
      return(NA)}
    else if(any(!is.finite(Var.mat))){
      return (NA)
    }
    else if( !is.positive.semi.definite(Var.mat)){return (NA)}
    else{
      #print(Var.mat)
      retval <- dmvnorm(x=ObsX,mean=ExpX, sigma=Var.mat,log=1)    
      if (!is.finite(retval))
      {return (NA)
      }else { return (retval) }
      
    }
  }
  else return(NA)
}


#Function to output negative log-likelihood for spec model given parameters
logll.spec <- function(x,other.data){
  if(all(is.finite(x)) & all(x>=0)){ 
    ObsX <- c(other.data$Obs$X_carcass, other.data$Obs$X_fHead, other.data$Obs$X_ovary, other.data$Obs$X_mHead, other.data$Obs$X_testis, other.data$Obs$X_accg)
    T.PC <- other.data$T.PC
    T.PCA <- other.data$T.PCA
    ExpX <- c()
    ExpX["carcass_P"] <- ExpP(ThetaP=x[1],ThetaA=x[3],Alpha=x[4],TimePC=T.PC)
    ExpX["carcass_C"] <- ExpC(ThetaC=x[2],ThetaA=x[3],Alpha=x[4],TimePC=T.PC)
    ExpX["carcass_A"] <- ExpA(ThetaA=x[3])
    
    ExpX["fHead_P"] <- ExpP(ThetaP=x[1],ThetaA=x[3],Alpha=x[5],TimePC=T.PC)
    ExpX["fHead_C"] <- ExpC(ThetaC=x[2],ThetaA=x[3],Alpha=x[5],TimePC=T.PC)
    ExpX["fHead_A"] <- ExpA(ThetaA=x[3])
    
    ExpX["ovary_P"] <- ExpP(ThetaP=x[1],ThetaA=x[3],Alpha=x[6],TimePC=T.PC)
    ExpX["ovary_C"] <- ExpC(ThetaC=x[2],ThetaA=x[3],Alpha=x[6],TimePC=T.PC)
    ExpX["ovary_A"] <- ExpA(ThetaA=x[3])
    
    ExpX["mHead_P"] <- ExpP(ThetaP=x[1],ThetaA=x[3],Alpha=x[7],TimePC=T.PC)
    ExpX["mHead_C"] <- ExpC(ThetaC=x[2],ThetaA=x[3],Alpha=x[7],TimePC=T.PC)
    ExpX["mHead_A"] <- ExpA(ThetaA=x[3])
    
    ExpX["testis_P"] <- ExpP(ThetaP=x[1],ThetaA=x[3],Alpha=x[8],TimePC=T.PC)
    ExpX["testis_C"] <- ExpC(ThetaC=x[2],ThetaA=x[3],Alpha=x[8],TimePC=T.PC)
    ExpX["testis_A"] <- ExpA(ThetaA=x[3])
    
    ExpX["accg_P"] <- ExpP(ThetaP=x[1],ThetaA=x[3],Alpha=x[9],TimePC=T.PC)
    ExpX["accg_C"] <- ExpC(ThetaC=x[2],ThetaA=x[3],Alpha=x[9],TimePC=T.PC)
    ExpX["accg_A"] <- ExpA(ThetaA=x[3])
    
    var <- sapply(4:9, function(i) x[10]/(2*x[i]))
    cov_t_t1 <- sapply(4:9, function(i) exp(- 2 * x[i] * T.PC) * var[i-3])
    cov_t_t2 <- sapply(4:9, function(i) exp(- 2 * x[i] * T.PC) * var[i-3])
    Var.mat <- diag(as.vector(sapply(var, function(x) rep(x,3))),nrow=length(ExpX),ncol=length(ExpX))
    j <- 1
    for (entry in seq(1,length(ExpX),by=3)){
      Var.mat[entry,entry+1] <- cov_t_t1[j]
      Var.mat[entry+1,entry] <- cov_t_t1[j]
      Var.mat[entry,entry+2] <- cov_t_t2[j]
      Var.mat[entry+2,entry] <- cov_t_t2[j]
      Var.mat[entry+1,entry+2] <- cov_t_t2[j]
      Var.mat[entry+2,entry+1] <- cov_t_t2[j]
      j <- j+ 1
    }
    if(any(!is.numeric(Var.mat))){
      return(NA)}
    else if(any(!is.finite(Var.mat))){
      return (NA)
    }
    else if( !is.positive.semi.definite(Var.mat)){return (NA)}
    else{
      #print(Var.mat)
      retval <- dmvnorm(x=ObsX,mean=ExpX, sigma=Var.mat,log=1)    
      if (!is.finite(retval))
      {return (NA)
      }else { return (retval) }
      
    }
  }
  else return(NA)
}

logll.sub <- function(x,other.data){
  if(all(is.finite(x)) & all(x>=0)){ 
    ObsX <- c(other.data$Obs$X_carcass, other.data$Obs$X_fHead, other.data$Obs$X_ovary, other.data$Obs$X_mHead, other.data$Obs$X_testis, other.data$Obs$X_accg)
    T.PC <- other.data$T.PC
    T.PCA <- other.data$T.PCA
    x <- c(x[1:2],x[1]+x[2],x[3:length(x)])
    ExpX <- c()
    ExpX["carcass_P"] <- ExpP(ThetaP=x[1],ThetaA=x[3],Alpha=x[4],TimePC=T.PC)
    ExpX["carcass_C"] <- ExpC(ThetaC=x[2],ThetaA=x[3],Alpha=x[4],TimePC=T.PC)
    ExpX["carcass_A"] <- ExpA(ThetaA=x[3])
    
    ExpX["fHead_P"] <- ExpP(ThetaP=x[1],ThetaA=x[3],Alpha=x[5],TimePC=T.PC)
    ExpX["fHead_C"] <- ExpC(ThetaC=x[2],ThetaA=x[3],Alpha=x[5],TimePC=T.PC)
    ExpX["fHead_A"] <- ExpA(ThetaA=x[3])
    
    ExpX["ovary_P"] <- ExpP(ThetaP=x[1],ThetaA=x[3],Alpha=x[6],TimePC=T.PC)
    ExpX["ovary_C"] <- ExpC(ThetaC=x[2],ThetaA=x[3],Alpha=x[6],TimePC=T.PC)
    ExpX["ovary_A"] <- ExpA(ThetaA=x[3])
    
    ExpX["mHead_P"] <- ExpP(ThetaP=x[1],ThetaA=x[3],Alpha=x[7],TimePC=T.PC)
    ExpX["mHead_C"] <- ExpC(ThetaC=x[2],ThetaA=x[3],Alpha=x[7],TimePC=T.PC)
    ExpX["mHead_A"] <- ExpA(ThetaA=x[3])
    
    ExpX["testis_P"] <- ExpP(ThetaP=x[1],ThetaA=x[3],Alpha=x[8],TimePC=T.PC)
    ExpX["testis_C"] <- ExpC(ThetaC=x[2],ThetaA=x[3],Alpha=x[8],TimePC=T.PC)
    ExpX["testis_A"] <- ExpA(ThetaA=x[3])
    
    ExpX["accg_P"] <- ExpP(ThetaP=x[1],ThetaA=x[3],Alpha=x[9],TimePC=T.PC)
    ExpX["accg_C"] <- ExpC(ThetaC=x[2],ThetaA=x[3],Alpha=x[9],TimePC=T.PC)
    ExpX["accg_A"] <- ExpA(ThetaA=x[3])
    
    var <- sapply(4:9, function(i) x[10]/(2*x[i]))
    cov_t_t1 <- sapply(4:9, function(i) exp(- 2 * x[i] * T.PC) * var[i-3])
    cov_t_t2 <- sapply(4:9, function(i) exp(- 2 * x[i] * T.PC) * var[i-3])
    Var.mat <- diag(as.vector(sapply(var, function(x) rep(x,3))),nrow=length(ExpX),ncol=length(ExpX))
    j <- 1
    for (entry in seq(1,length(ExpX),by=3)){
      Var.mat[entry,entry+1] <- cov_t_t1[j]
      Var.mat[entry+1,entry] <- cov_t_t1[j]
      Var.mat[entry,entry+2] <- cov_t_t2[j]
      Var.mat[entry+2,entry] <- cov_t_t2[j]
      Var.mat[entry+1,entry+2] <- cov_t_t2[j]
      Var.mat[entry+2,entry+1] <- cov_t_t2[j]
      j <- j+ 1
    }
    if(any(!is.numeric(Var.mat))){
      return(NA)}
    else if(any(!is.finite(Var.mat))){
      return (NA)
    }
    else if( !is.positive.semi.definite(Var.mat)){return (NA)}
    else{
      retval <- dmvnorm(x=ObsX,mean=ExpX, sigma=Var.mat,log=1)    
      if (!is.finite(retval))
      {return (NA)
      }else { return (retval) }
    }
  }
  else return(NA)
}

#Function to output negative log-likelihood for cons model given parameters
logll.cons <- function(x, other.data){
  if(all(is.finite(x)) & all(x>=0)){ 
    ObsX <- c(other.data$Obs$X_carcass, other.data$Obs$X_fHead, other.data$Obs$X_ovary, other.data$Obs$X_mHead, other.data$Obs$X_testis, other.data$Obs$X_accg)
    T.PC <- other.data$T.PC
    T.PCA <- other.data$T.PCA
    ExpX <- c()
    ExpX["carcass_P"] <- ExpP(ThetaP=x[1],ThetaA=x[1],Alpha=x[2],TimePC=T.PC)
    ExpX["carcass_C"] <- ExpC(ThetaC=x[1],ThetaA=x[1],Alpha=x[2],TimePC=T.PC)
    ExpX["carcass_A"] <- ExpA(ThetaA=x[1])
    
    ExpX["fHead_P"] <- ExpP(ThetaP=x[1],ThetaA=x[1],Alpha=x[3],TimePC=T.PC)
    ExpX["fHead_C"] <- ExpC(ThetaC=x[1],ThetaA=x[1],Alpha=x[3],TimePC=T.PC)
    ExpX["fHead_A"] <- ExpA(ThetaA=x[1])
    
    ExpX["ovary_P"] <- ExpP(ThetaP=x[1],ThetaA=x[1],Alpha=x[4],TimePC=T.PC)
    ExpX["ovary_C"] <- ExpC(ThetaC=x[1],ThetaA=x[1],Alpha=x[4],TimePC=T.PC)
    ExpX["ovary_A"] <- ExpA(ThetaA=x[1])
    
    ExpX["mHead_P"] <- ExpP(ThetaP=x[1],ThetaA=x[1],Alpha=x[5],TimePC=T.PC)
    ExpX["mHead_C"] <- ExpC(ThetaC=x[1],ThetaA=x[1],Alpha=x[5],TimePC=T.PC)
    ExpX["mHead_A"] <- ExpA(ThetaA=x[1])
    
    ExpX["testis_P"] <- ExpP(ThetaP=x[1],ThetaA=x[1],Alpha=x[6],TimePC=T.PC)
    ExpX["testis_C"] <- ExpC(ThetaC=x[1],ThetaA=x[1],Alpha=x[6],TimePC=T.PC)
    ExpX["testis_A"] <- ExpA(ThetaA=x[1])
    
    ExpX["accg_P"] <- ExpP(ThetaP=x[1],ThetaA=x[1],Alpha=x[7],TimePC=T.PC)
    ExpX["accg_C"] <- ExpC(ThetaC=x[1],ThetaA=x[1],Alpha=x[7],TimePC=T.PC)
    ExpX["accg_A"] <- ExpA(ThetaA=x[1])
    var <- sapply(2:7, function(i) x[8]/(2*x[i]))
    cov_t_t1 <- sapply(2:7, function(i) exp(- 2 * x[i] * T.PC) * var[i-1])
    cov_t_t2 <- sapply(2:7, function(i) exp(- 2 * x[i] * T.PC) * var[i-1])
    Var.mat <- diag(as.vector(sapply(var, function(x) rep(x,3))),nrow=length(ExpX),ncol=length(ExpX))
    j <- 1
    for (entry in seq(1,length(ExpX),by=3)){
      Var.mat[entry,entry+1] <- cov_t_t1[j]
      Var.mat[entry+1,entry] <- cov_t_t1[j]
      Var.mat[entry,entry+2] <- cov_t_t2[j]
      Var.mat[entry+2,entry] <- cov_t_t2[j]
      Var.mat[entry+1,entry+2] <- cov_t_t2[j]
      Var.mat[entry+2,entry+1] <- cov_t_t2[j]
      j <- j+ 1
    }
    if(any(!is.numeric(Var.mat))){
      return(NA)}
    else if(any(!is.finite(Var.mat))){
      return (NA)
    }
    else if( !is.positive.semi.definite(Var.mat)){return (NA)}
    else{
      #print(Var.mat)
      retval <- dmvnorm(x=ObsX,mean=ExpX, sigma=Var.mat,log=1)  
      if (!is.finite(retval))
      {return (NA)
      }else { return (retval) }
      
    }
  }
  else return(NA)
}


#clear screen and time my code
cat("\014")  
t1 <- Sys.time()
eps = 0
sig=0.05
sig2 <- 0.05
# Generating Simulated Data
PC_ind = c(rep(c(1,3),times=6) + rep(3*(0:5),each=2))
Thetamin = -3
Thetamax = 3
Alphamax = 3
Sigma2max = 3
load("Mean_Cor_Real_Data.Rdata")

data = read.csv("OU_table_drosophila_abs", sep = "")
indices = 1:nrow(data)
#Getting data of all genes in mel and pse for normalizing between A and P,C.
mel = read.table("mel_ExpressionTable.norm.newIDs",sep="\t",header=F,row.names=NULL,skip=1)
pse = read.table("pse_ExpressionTable.norm.newIDs",sep="\t",header=T,row.names=NULL)
mel$V11=NULL
colnames(mel) = colnames(pse)
med.norm = c()
med.norm[1] = median(pse$carcass) - median(mel$carcass)
med.norm[2] = median(pse$fHead) - median(mel$fHead)
med.norm[3] = median(pse$ovary) - median(mel$ovary)
med.norm[4] = median(pse$mHead) - median(mel$mHead)
med.norm[5] = median(pse$testis) - median(mel$testis)
med.norm[6] = median(pse$accg) - median(mel$accg)

#Getting the data
TimesPC <- data[indices,9]
TimesPCA <- data[indices,10]
TimesPC[TimesPC==0] <- 10**-5
TimesPCA[TimesPCA==0] <- 10**-5

ObsXP <- data[indices,17]
ObsXC <- data[indices,11]
ObsXA <- data[indices,23]
ObsX_carcass <- cbind("ObsXP"=ObsXP,"ObsXC"= ObsXC,"ObsXA"= ObsXA)

ObsXP <- data[indices,18]
ObsXC <- data[indices,12]
ObsXA <- data[indices,24]
ObsX_fHead <- cbind("ObsXP"=ObsXP,"ObsXC"= ObsXC,"ObsXA"= ObsXA)

ObsXP <- data[indices,19]
ObsXC <- data[indices,13]
ObsXA <- data[indices,25]
ObsX_ovary <- cbind("ObsXP"=ObsXP,"ObsXC"= ObsXC,"ObsXA"= ObsXA)

ObsXP <- data[indices,20]
ObsXC <- data[indices,14]
ObsXA <- data[indices,26]
ObsX_mHead <- cbind("ObsXP"=ObsXP,"ObsXC"= ObsXC,"ObsXA"= ObsXA)

ObsXP <- data[indices,21]
ObsXC <- data[indices,15]
ObsXA <- data[indices,27]
ObsX_testis <- cbind("ObsXP"=ObsXP,"ObsXC"= ObsXC,"ObsXA"= ObsXA)

ObsXP <- data[indices,22]
ObsXC <- data[indices,16]
ObsXA <- data[indices,28]
ObsX_accg <- cbind("ObsXP"=ObsXP,"ObsXC"= ObsXC,"ObsXA"= ObsXA)

X_all <- list("carcass"=ObsX_carcass,"fHead"=ObsX_fHead,"ovary"=ObsX_ovary,"mHead"=ObsX_mHead,"testis"=ObsX_testis,"accg"=ObsX_accg)

# #Normalizing by median within each of the tissues:
for(i in 1:length(X_all)){
  for(j in 1:2){X_all[[i]][,j] = X_all[[i]][,j] + med.norm[i]}
}


#Library size normalization
Xmat = X_all[[1]]
for(i in 2:length(X_all)){
  Xmat = cbind(Xmat,X_all[[i]])
}

cnorm <- calcNormFactors(Xmat,method="upperquartile",p=0.75)
for(i in 1:length(cnorm)){
  Xmat[,i] = Xmat[,i]*cnorm[i]
}
X_all <- list("carcass"=Xmat[,1:3],"fHead"=Xmat[,4:6],"ovary"=Xmat[,7:9],"mHead"=Xmat[,10:12],"testis"=Xmat[,13:15],"accg"=Xmat[,16:18])

# #Transforming the data
# for(i in 1:length(X_all)){
#   for(j in 1:3){
#     X_all[[i]][,j] = log(X_all[[i]][,j]+1)
#   }
# }
data.class = as.character(data[,1])
data.class = as.factor(data.class)
summary(data.class)


#Data.frame to store the results of Parameters and gene names
Thetas <- c("ThetaP","ThetaC","ThetaA")
tissues <- c("carcass","fHead","ovary","mHead","testis","accg")
Alphas <- paste("Alpha",tissues,sep="_")
All_Params_neoc <- c(Thetas[-3],Alphas)
cols.names <- c("Parent","Child","Ancestor",All_Params_neoc,"Sigma2", "log.likelihood","Pvalues_cons")
Neochild.Model <- matrix(data=0,nrow=length(indices),ncol=length(cols.names))
colnames(Neochild.Model) <- cols.names
Neochild.Model <- data.frame(Neochild.Model)

All_Params_neop <- c(Thetas[-3],Alphas)
cols.names <- c(All_Params_neop,"Sigma2", "log.likelihood","Pvalues_cons")
Neoparent.Model <- matrix(data=0,nrow=length(indices),ncol=length(cols.names))
colnames(Neoparent.Model) <- cols.names
Neoparent.Model <- data.frame(Neoparent.Model)

All_Params_cons <- c("Theta",Alphas)
cols_cons <- c(All_Params_cons,"Sigma2","log.likelihood")
Cons.Model <- matrix(data=0,nrow=length(indices),ncol=length(cols_cons))
colnames(Cons.Model) <- cols_cons
Cons.Model <- data.frame(Cons.Model)

# Thetas <- as.vector(t(sapply(Params,function(x) paste(x,tissues,sep="_"))))
All_Params_spec <- c(Thetas,Alphas)
cols.names <- c(All_Params_spec, "Sigma2", "log.likelihood", "Pvalues_neoc","Pvalues_neop","Pvalues_cons","Pvalues_sub")
Spec.Model <- matrix(data=0,nrow=length(indices),ncol=length(cols.names))
colnames(Spec.Model) <- cols.names
Spec.Model <- data.frame(Spec.Model)

All_Params_sub <- c(Thetas[-3],Alphas)
cols.names <- c(All_Params_sub, "Sigma2", "log.likelihood")
Sub.Model <- matrix(data=0,nrow=length(indices),ncol=length(cols.names))
colnames(Sub.Model) <- cols.names
Sub.Model <- data.frame(Sub.Model)

eps <- .Machine$double.eps
cons.lower <- c(0,rep(eps,7))
neochild.lower <- c(0,0,rep(eps,7))
neoparent.lower <- c(0,0,rep(eps,7))
spec.lower <- c(0,0,0,rep(eps,7))
# cons.lower <- c(0,rep(eps,2))
# neochild.lower <- c(0,0,rep(eps,2))
# neoparent.lower <- c(0,0,rep(eps,2))
# sub.lower <- c(0,0,rep(eps,2))
# spec.lower <- c(0,0,0,rep(eps,2))
up.a <- 10**12
up.theta <- log10(range(X_all)[2]+1)
low.theta <- log10(range(X_all)[1]+1)

# cons.upper <- c(up.a,rep(up.a,2))
# neochild.upper <- c(up.a,up.a,rep(up.a,2))
# neoparent.upper <- c(up.a,up.a,rep(up.a,2))
# neoparent.upper <- c(up.a,up.a,rep(up.a,2))
# spec.upper <- c(up.a,up.a,up.a,rep(up.a,2))
cons.upper <- c(up.theta,rep(up.a,7))
neochild.upper <- c(up.theta,up.theta,rep(up.a,7))
neoparent.upper <- c(up.theta,up.theta,rep(up.a,7))
spec.upper <- c(up.theta,up.theta,up.theta,rep(up.a,7))
n <- 100
Pred.Class = rep("cons",length(indices))
sig <- 0.05/3
sig2 <- 0.05
#Setting check variables
#Looping over each neochild duplicate pair

get.data <- function(i){
  Obs <- list("X_carcass"=X_all$carcass[i,],"X_fHead"=X_all$fHead[i,],"X_ovary"=X_all$ovary[i,],"X_mHead"=X_all$mHead[i,],"X_testis"=X_all$testis[i,],"X_accg"=X_all$accg[i,])
  other.data <-  list("Obs"=Obs,"T.PC"=TimesPC[i],"T.PCA"=TimesPCA[i])
  return(other.data)
}

no_cores = detectCores(TRUE,TRUE)-1
cl=makeCluster(no_cores,type="FORK")
varlist2 =  c("i","neochild.choice","neoparent.choice","cons.choice","spec.choice","sub.choice")
for(i in 1:length(indices)){
  set.seed(NULL)
  #choices is a list to compare and get the best answer from optim
  neochild.choice <- rep(0,length(neochild.lower)+1)
  cons.choice <- rep(0,length(cons.lower)+1)
  neoparent.choice <- rep(0,length(neoparent.lower)+1)
  sub.choice <- rep(0,length(spec.lower))
  spec.choice <- rep(0,length(spec.lower)+1)
  #Running optim n times to try to get a global minimium
  clusterExport(cl, varlist = varlist2,envir=.GlobalEnv)
  clusterSetRNGStream(cl,iseed=NULL)
  neochild.runs <- clusterApply(cl,1:n,function(j){
    while(1){
      neochild.init <- c(10**runif(2,low.theta,up.theta),10**runif(7,-Alphamax,Alphamax))
      if (!is.na(logll.neochild(neochild.init,get.data(i)))){break}
    }
    optim(neochild.init, fn=logll.neochild,other.data=get.data(i), method="Nelder-Mead",control=list("fnscale"=-1))
  })
  
  clusterSetRNGStream(cl,iseed=NULL)
  neoparent.runs <- clusterApply(cl,1:n, function(j){
    while(1){
      neoparent.init <- c(10**runif(2,low.theta,up.theta),10**runif(7,-Alphamax,Alphamax))
      if (!is.na(logll.neoparent(neoparent.init,get.data(i)))){break}
    }
    optim(par=neoparent.init, fn=logll.neoparent,other.data=get.data(i), method="Nelder-Mead",control=list("fnscale"=-1))
  })
  
  clusterSetRNGStream(cl,iseed=NULL)
  cons.runs <- clusterApply(cl,1:n,function(j){
    while(1){
      cons.init <- c(10**runif(1,low.theta,up.theta),10**runif(7,-Alphamax,Alphamax))
      if (!is.na(logll.cons(cons.init,get.data(i)))){break}
    }
    optim(par=cons.init, fn=logll.cons,other.data=get.data(i), method="Nelder-Mead",control=list("fnscale"=-1))
  })
  
  clusterSetRNGStream(cl,iseed=NULL)
  spec.runs <- clusterApply(cl,1:n, function(j){
    while(1){
      spec.init <- c(10**runif(3,low.theta,up.theta),10**runif(7,-Alphamax,Alphamax))
      if (!is.na(logll.spec(spec.init,get.data(i)))){break}
    }
    optim(spec.init, fn=logll.spec,other.data=get.data(i), method="Nelder-Mead",control=list("fnscale"=-1))
  })
  
  
  
  gc()
  neochild.runs=matrix(sapply(1:length(neochild.runs), function(x) unlist(neochild.runs[[x]][1:2])),ncol=length(neochild.lower)+1,byrow=T)
  neoc.ind = which.max(neochild.runs[,length(neochild.lower)+1])
  neochild.choice[length(neochild.choice)] = neochild.runs[neoc.ind,length(neochild.lower)+1]
  neochild.choice[1:length(neochild.lower)] = neochild.runs[neoc.ind,1:length(neochild.lower)]
  
  neoparent.runs=matrix(sapply(1:length(neoparent.runs), function(x) unlist(neoparent.runs[[x]][1:2])),ncol=length(neoparent.lower)+1,byrow=T)
  neop.ind = which.max(neoparent.runs[,length(neoparent.lower)+1])
  neoparent.choice[length(neoparent.choice)] = neoparent.runs[neop.ind,length(neoparent.lower)+1]
  neoparent.choice[1:length(neoparent.lower)] = neoparent.runs[neop.ind,1:length(neoparent.lower)]
  
  cons.runs=matrix(sapply(1:length(cons.runs), function(x) unlist(cons.runs[[x]][1:2])),ncol=length(cons.lower)+1, byrow=T)
  cons.ind = which.max(cons.runs[,length(cons.lower)+1])
  cons.choice[length(cons.choice)] = cons.runs[cons.ind,length(cons.lower)+1]
  cons.choice[1:length(cons.lower)] = cons.runs[cons.ind,1:length(cons.lower)]
  
  spec.runs=matrix(sapply(1:length(spec.runs), function(x) unlist(spec.runs[[x]][1:2])),ncol=length(spec.lower)+1, byrow=T)
  spec.ind = which.max(spec.runs[,length(spec.lower)+1])
  spec.choice[length(spec.choice)] = spec.runs[spec.ind,length(spec.lower)+1]
  spec.choice[1:length(spec.lower)] = spec.runs[spec.ind,1:length(spec.lower)]
  
  Neochild.Model[i,c(All_Params_neoc,"Sigma2","log.likelihood")] <- neochild.choice
  Neoparent.Model[i,c(All_Params_neop,"Sigma2","log.likelihood")] <- neoparent.choice
  Cons.Model[i,c(All_Params_cons,"Sigma2","log.likelihood")] <- cons.choice
  Spec.Model[i,c(All_Params_spec,"Sigma2","log.likelihood")] <- spec.choice
  
  
  LR.neoc_v_cons <- 2*(Neochild.Model[i,"log.likelihood"]  - Cons.Model[i,"log.likelihood"])
  LR.neop_v_cons <- 2*(Neoparent.Model[i,"log.likelihood"]  - Cons.Model[i,"log.likelihood"])
  LR.spec_v_neoc <- 2*(Spec.Model[i,"log.likelihood"]  - Neochild.Model[i,"log.likelihood"])
  LR.spec_v_neop <- 2*(Spec.Model[i,"log.likelihood"] - Neoparent.Model[i,"log.likelihood"])
  LR.spec_v_cons <- 2*(Spec.Model[i,"log.likelihood"] - Cons.Model[i,"log.likelihood"])
  
  #Resetting Pvalues of 0 to eps
  Pval.neoc_v_cons <- pchisq(LR.neoc_v_cons,df=1,lower.tail = FALSE)
  if(Pval.neoc_v_cons==0) Pval.neoc_v_cons=eps
  Pval.neop_v_cons <- pchisq(LR.neop_v_cons,df=1,lower.tail = FALSE)
  if(Pval.neop_v_cons==0) Pval.neop_v_cons=eps
  Pval.spec_v_neoc <- pchisq(LR.spec_v_neoc,df=1,lower.tail=FALSE)
  if(Pval.spec_v_neoc==0) Pval.spec_v_neoc=eps
  Pval.spec_v_neop <- pchisq(LR.spec_v_neop,df=1,lower.tail=FALSE)
  if(Pval.spec_v_neop==0) Pval.spec_v_neop=eps
  Pval.spec_v_cons <- pchisq(LR.spec_v_cons,df=2,lower.tail=FALSE)
  if(Pval.spec_v_cons==0) Pval.spec_v_cons=eps
  
  Neochild.Model[i,"Pvalues_cons"] <- Pval.neoc_v_cons
  Neoparent.Model[i,"Pvalues_cons"] <- Pval.neop_v_cons
  Spec.Model[i,"Pvalues_neoc"] <- Pval.spec_v_neoc
  Spec.Model[i,"Pvalues_neop"] <- Pval.spec_v_neop
  Spec.Model[i,"Pvalues_cons"] <- Pval.spec_v_cons
  
  #Predictions
  if(Neochild.Model$Pvalues_cons[i]<sig & Neochild.Model$Pvalues_cons[i]>0 & Neoparent.Model$Pvalues_cons[i]>=sig & Spec.Model$Pvalues_neoc[i]>=sig){
    Pred.Class[i]= "neochild"}
  if(Neochild.Model$Pvalues_cons[i]>=sig & Neoparent.Model$Pvalues_cons[i]<sig & Neoparent.Model$Pvalues_cons[i]>0 & Spec.Model$Pvalues_neop[i]>=sig){
    Pred.Class[i]= "neoparent"}
  if(Spec.Model$Pvalues_neoc[i]<sig & Spec.Model$Pvalues_neop[i]<sig & Spec.Model$Pvalues_cons[i]<sig & Spec.Model$Pvalues_neoc[i]>0 & Spec.Model$Pvalues_neop[i]>0 & Spec.Model$Pvalues_cons[i]>0){
    Pred.Class[i]= "spec"
    
    clusterSetRNGStream(cl,iseed=NULL)
    sub.runs <- clusterApply(cl,1:n, function(j){
      while(1){
        sub.init <- c(10**runif(2,low.theta,up.theta),10**runif(7,-Alphamax,Alphamax))
        if (!is.na(logll.sub(sub.init,get.data(i)))){break}
      }
      optim(sub.init, fn=logll.sub,other.data=get.data(i), method="Nelder-Mead",control=list("fnscale"=-1))
    })
    sub.runs=matrix(sapply(1:length(sub.runs), function(x) unlist(sub.runs[[x]][1:2])),ncol=length(sub.choice), byrow=T)
    sub.ind = which.max(sub.runs[,length(sub.choice)])
    sub.choice = sub.runs[sub.ind,]
    Sub.Model[i,c(All_Params_sub,"Sigma2","log.likelihood")] <- sub.choice
    LR.spec_v_sub  <- 2*(Spec.Model[i,"log.likelihood"] - Sub.Model[i,"log.likelihood"])
    Pval.spec_v_sub <- pchisq(LR.spec_v_sub,df=1,lower.tail=FALSE)
    if(Pval.spec_v_sub==0) Pval.spec_v_sub=eps
    Spec.Model[i,"Pvalues_sub"] <- Pval.spec_v_sub
    if(Pval.spec_v_sub>sig2){Pred.Class[i]= "sub"}
  }
  print(i)
}
stopCluster(cl)
gc()

Pred.Class = factor(Pred.Class,levels=c("cons","neochild","neoparent","spec","sub"))
print(table(Pred.Class,data.class))
save(Pred.Class,data.class,Neochild.Model,Neoparent.Model,Cons.Model, Spec.Model,Sub.Model, file="Est_Models_A_v12_nolog.Rdata")
t2 <- Sys.time()
print(t2 - t1)
cat("Exited Successfully!","\n")
# Version12: Alphamax and Sigma2max are set to 3
