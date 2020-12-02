rm(list=ls())
library(MASS)
library(mvtnorm)
library(matrixcalc)
library(rlecuyer)
library(parallel)
# library(CDROM)
# ,lib.loc="/storage/home/bzk18/R/x86_64-redhat-linux-gnu-library/3.1"
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

#Funciton to output negative log-likelihood for neochild model given parameters
logll.neochild <- function(x,other.data,cor.mat=mean_cor){
  if(all(is.finite(x)) & all(x>=0)){ 
    ObsX <- c(other.data$Obs$X_carcass, other.data$Obs$X_fHead, other.data$Obs$X_ovary, other.data$Obs$X_mHead, other.data$Obs$X_testis, other.data$Obs$X_accg)
    # ObsX = other.data$Obs
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
    Var.mat <- diag(var,nrow=length(ExpX),ncol=length(ExpX))
    entry<- 1
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
    j <- 1
    for(entry in seq(1,length(ExpX)-3,3)){
      k<-1
      for(a in seq(entry+3,length(ExpX),3)){
        Var.mat[entry,a] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry,a+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry,a+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+1,a] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+1,a+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+1,a+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+2,a] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+2,a+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+2,a+2] <- cor.mat[j,j+k] * var[j]
        
        Var.mat[a,entry] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+1,entry] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+2,entry] <- cor.mat[j,j+k] * var[j]
        Var.mat[a,entry+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+1,entry+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+2,entry+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[a,entry+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+1,entry+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+2,entry+2] <- cor.mat[j,j+k] * var[j]
        
        k <- k+1
      }
      j <- j+1
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

#Function to return negative log of MLE for neoparent model
logll.neoparent <- function(x,other.data,cor.mat=mean_cor){
  if(all(is.finite(x)) & all(x>=0)){
    ObsX <- c(other.data$Obs$X_carcass, other.data$Obs$X_fHead, other.data$Obs$X_ovary, other.data$Obs$X_mHead, other.data$Obs$X_testis, other.data$Obs$X_accg)
    # ObsX = other.data$Obs
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
    Var.mat <- diag(var,nrow=length(ExpX),ncol=length(ExpX))
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
    j <- 1
    for(entry in seq(1,length(ExpX)-3,3)){
      k<-1
      for(a in seq(entry+3,length(ExpX),3)){
        Var.mat[entry,a] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry,a+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry,a+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+1,a] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+1,a+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+1,a+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+2,a] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+2,a+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+2,a+2] <- cor.mat[j,j+k] * var[j]
        
        Var.mat[a,entry] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+1,entry] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+2,entry] <- cor.mat[j,j+k] * var[j]
        Var.mat[a,entry+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+1,entry+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+2,entry+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[a,entry+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+1,entry+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+2,entry+2] <- cor.mat[j,j+k] * var[j]
        
        k <- k+1
      }
      j <- j+1
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
  }else return(NA)
}

#Function to return log of MLE for Subfunctionalization model
logll.sub <- function(x,other.data,cor.mat=mean_cor){
  if(all(is.finite(x)) & all(x>=0)){ 
    ObsX <- c(other.data$Obs$X_carcass, other.data$Obs$X_fHead, other.data$Obs$X_ovary, other.data$Obs$X_mHead, other.data$Obs$X_testis, other.data$Obs$X_accg)
    # ObsX = other.data$Obs
    T.PC <- other.data$T.PC
    T.PCA <- other.data$T.PCA
    x <- c(x[1:2],x[1]+x[2],x[3:4])
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
    Var.mat <- diag(var,nrow=length(ExpX),ncol=length(ExpX))
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
    j <- 1
    for(entry in seq(1,length(ExpX)-3,3)){
      k<-1
      for(a in seq(entry+3,length(ExpX),3)){
        Var.mat[entry,a] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry,a+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry,a+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+1,a] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+1,a+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+1,a+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+2,a] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+2,a+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+2,a+2] <- cor.mat[j,j+k] * var[j]
        
        Var.mat[a,entry] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+1,entry] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+2,entry] <- cor.mat[j,j+k] * var[j]
        Var.mat[a,entry+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+1,entry+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+2,entry+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[a,entry+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+1,entry+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+2,entry+2] <- cor.mat[j,j+k] * var[j]
        
        k <- k+1
      }
      j <- j+1
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


#Function to output log-likelihood for specialization model given parameters
logll.spec <- function(x,other.data,cor.mat=mean_cor){
  if(all(is.finite(x)) & all(x>=0)){ 
    ObsX <- c(other.data$Obs$X_carcass, other.data$Obs$X_fHead, other.data$Obs$X_ovary, other.data$Obs$X_mHead, other.data$Obs$X_testis, other.data$Obs$X_accg)
    # ObsX = other.data$Obs
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
    Var.mat <- diag(var,nrow=length(ExpX),ncol=length(ExpX))
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
    j <- 1
    for(entry in seq(1,length(ExpX)-3,3)){
      k<-1
      for(a in seq(entry+3,length(ExpX),3)){
        Var.mat[entry,a] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry,a+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry,a+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+1,a] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+1,a+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+1,a+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+2,a] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+2,a+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+2,a+2] <- cor.mat[j,j+k] * var[j]
        
        Var.mat[a,entry] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+1,entry] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+2,entry] <- cor.mat[j,j+k] * var[j]
        Var.mat[a,entry+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+1,entry+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+2,entry+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[a,entry+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+1,entry+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+2,entry+2] <- cor.mat[j,j+k] * var[j]
        
        k <- k+1
      }
      j <- j+1
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
logll.cons <- function(x,other.data,cor.mat=mean_cor){
  if(all(is.finite(x)) & all(x>=0)){ 
    ObsX <- c(other.data$Obs$X_carcass, other.data$Obs$X_fHead, other.data$Obs$X_ovary, other.data$Obs$X_mHead, other.data$Obs$X_testis, other.data$Obs$X_accg)
    # ObsX = other.data$Obs
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
    j <- 1
    for(entry in seq(1,length(ExpX)-3,3)){
      k<-1
      for(a in seq(entry+3,length(ExpX),3)){
        Var.mat[entry,a] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry,a+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry,a+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+1,a] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+1,a+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+1,a+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+2,a] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+2,a+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+2,a+2] <- cor.mat[j,j+k] * var[j]
        
        Var.mat[a,entry] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+1,entry] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+2,entry] <- cor.mat[j,j+k] * var[j]
        Var.mat[a,entry+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+1,entry+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+2,entry+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[a,entry+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+1,entry+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+2,entry+2] <- cor.mat[j,j+k] * var[j]
        
        k <- k+1
      }
      j <- j+1
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

Obs.neochild <- function(x,cor.mat=mean_cor){
  T.PC <- x$T.PC
  T.PCA <- x$T.PCA
  x = c(x$P,x$C,x$Alpha,x$Sigma2)
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
  Var.mat <- diag(var,nrow=length(ExpX),ncol=length(ExpX))
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
  j <- 1
  for(entry in seq(1,length(ExpX)-3,3)){
    k<-1
    for(a in seq(entry+3,length(ExpX),3)){
      Var.mat[entry,a] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry,a+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry,a+2] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+1,a] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+1,a+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+1,a+2] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+2,a] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+2,a+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+2,a+2] <- cor.mat[j,j+k] * var[j]
      
      Var.mat[a,entry] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+1,entry] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+2,entry] <- cor.mat[j,j+k] * var[j]
      Var.mat[a,entry+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+1,entry+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+2,entry+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[a,entry+2] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+1,entry+2] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+2,entry+2] <- cor.mat[j,j+k] * var[j]
      
      k <- k+1
    }
    j <- j+1
  }
  if(any(!is.numeric(Var.mat))){
    return(NA)}
  else if(any(!is.finite(Var.mat))){
    return (NA)
  }
  else if( !is.positive.semi.definite(Var.mat)){return (NA)}
  else{
    retval <- dmvnorm(x=ExpX,mean=ExpX, sigma=Var.mat,log=1)    
    if (!is.finite(retval))
    {return (NA)
    }else { return(c(ExpX,retval)) }
  }
  
}

#Function to return negative log of MLE for neoparent model
Obs.neoparent <- function(x,cor.mat=mean_cor){
  T.PC <- x$T.PC
  T.PCA <- x$T.PCA
  x = c(x$P,x$C,x$Alpha,x$Sigma2)
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
  Var.mat <- diag(var,nrow=length(ExpX),ncol=length(ExpX))
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
  j <- 1
  for(entry in seq(1,length(ExpX)-3,3)){
    k<-1
    for(a in seq(entry+3,length(ExpX),3)){
      Var.mat[entry,a] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry,a+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry,a+2] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+1,a] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+1,a+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+1,a+2] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+2,a] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+2,a+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+2,a+2] <- cor.mat[j,j+k] * var[j]
      
      Var.mat[a,entry] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+1,entry] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+2,entry] <- cor.mat[j,j+k] * var[j]
      Var.mat[a,entry+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+1,entry+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+2,entry+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[a,entry+2] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+1,entry+2] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+2,entry+2] <- cor.mat[j,j+k] * var[j]
      
      k <- k+1
    }
    j <- j+1
  }
  if(any(!is.numeric(Var.mat))){
    return(NA)}
  else if(any(!is.finite(Var.mat))){
    return (NA)
  }
  else if( !is.positive.semi.definite(Var.mat)){return (NA)}
  else{
    retval <- dmvnorm(x=ExpX,mean=ExpX, sigma=Var.mat,log=1)    
    if (!is.finite(retval))
    {return (NA)
    }else { return(c(ExpX,retval)) }
    
  }
}

#Function to output negative log-likelihood for spec model given parameters
Obs.spec <- function(x,cor.mat=mean_cor){
  T.PC <- x$T.PC
  T.PCA <- x$T.PCA
  x = c(x$P,x$C,x$A,x$Alpha,x$Sigma2)
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
  Var.mat <- diag(var,nrow=length(ExpX),ncol=length(ExpX))
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
  j <- 1
  for(entry in seq(1,length(ExpX)-3,3)){
    k<-1
    for(a in seq(entry+3,length(ExpX),3)){
      Var.mat[entry,a] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry,a+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry,a+2] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+1,a] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+1,a+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+1,a+2] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+2,a] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+2,a+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+2,a+2] <- cor.mat[j,j+k] * var[j]
      
      Var.mat[a,entry] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+1,entry] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+2,entry] <- cor.mat[j,j+k] * var[j]
      Var.mat[a,entry+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+1,entry+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+2,entry+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[a,entry+2] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+1,entry+2] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+2,entry+2] <- cor.mat[j,j+k] * var[j]
      
      k <- k+1
    }
    j <- j+1
  }
  if(any(!is.numeric(Var.mat))){
    return(NA)}
  else if(any(!is.finite(Var.mat))){
    return (NA)
  }
  else if( !is.positive.semi.definite(Var.mat)){return (NA)}
  else{
    retval <- dmvnorm(x=ExpX,mean=ExpX, sigma=Var.mat,log=1)    
    if (!is.finite(retval))
    {return (NA)
    }else { return(c(ExpX,retval)) }
  }
  
}

#Function to optimize Theta3 for MLE for simulating data for Sub class
get.spec.theta <- function(x,other.data,cor.mat=mean_cor){
  # & abs(sum(other.data$params[1:2])-x)<=1
  if(is.finite(x) & (x>=0) ){ 
    ObsX = other.data$Obs
    T.PC <- other.data$T.PC
    T.PCA <- other.data$T.PCA
    x = c(other.data$params[1:2],x,other.data$params[3:4])
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
    Var.mat <- diag(var,nrow=length(ExpX),ncol=length(ExpX))
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
    j <- 1
    for(entry in seq(1,length(ExpX)-3,3)){
      k<-1
      for(a in seq(entry+3,length(ExpX),3)){
        Var.mat[entry,a] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry,a+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry,a+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+1,a] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+1,a+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+1,a+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+2,a] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+2,a+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[entry+2,a+2] <- cor.mat[j,j+k] * var[j]
        
        Var.mat[a,entry] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+1,entry] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+2,entry] <- cor.mat[j,j+k] * var[j]
        Var.mat[a,entry+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+1,entry+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+2,entry+1] <- cor.mat[j,j+k] * var[j]
        Var.mat[a,entry+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+1,entry+2] <- cor.mat[j,j+k] * var[j]
        Var.mat[a+2,entry+2] <- cor.mat[j,j+k] * var[j]
        
        k <- k+1
      }
      j <- j+1
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
Obs.cons <- function(x,cor.mat=mean_cor){
  T.PC <- x$T.PC
  T.PCA <- x$T.PCA
  x = c(x$P,x$Alpha,x$Sigma2)
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
  j <- 1
  for(entry in seq(1,length(ExpX)-3,3)){
    k<-1
    for(a in seq(entry+3,length(ExpX),3)){
      Var.mat[entry,a] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry,a+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry,a+2] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+1,a] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+1,a+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+1,a+2] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+2,a] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+2,a+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[entry+2,a+2] <- cor.mat[j,j+k] * var[j]
      
      Var.mat[a,entry] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+1,entry] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+2,entry] <- cor.mat[j,j+k] * var[j]
      Var.mat[a,entry+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+1,entry+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+2,entry+1] <- cor.mat[j,j+k] * var[j]
      Var.mat[a,entry+2] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+1,entry+2] <- cor.mat[j,j+k] * var[j]
      Var.mat[a+2,entry+2] <- cor.mat[j,j+k] * var[j]
      
      k <- k+1
    }
    j <- j+1
  }
  if(any(!is.numeric(Var.mat))){
    return(NA)}
  else if(any(!is.finite(Var.mat))){
    return (NA)
  }
  else if( !is.positive.semi.definite(Var.mat)){return (NA)}
  else{
    retval <- dmvnorm(x=ExpX,mean=ExpX, sigma=Var.mat,log=1)  
    if (!is.finite(retval))
    {return (NA)
    }else { return(c(ExpX,retval)) }
  }
  
}
#clear screen and time my code
cat("\014")  
t1 <- Sys.time()
eps = 0
indices = 1:400
sig=0.05
sig2 <- 0.05
# Generating Simulated Data
X = matrix(0,nrow=length(indices),ncol=3)
TimesPC = rep(0,length(indices))
TimesPCA = rep(0,length(indices))
Sim.Alphas = matrix(0,nrow=length(indices),ncol=6)
Sim.Sigmas2 = rep(0,length(indices))
Sim.ThetaP = rep(0,length(indices))
Sim.ThetaC = rep(0,length(indices))
Sim.ThetaA = rep(0,length(indices))
PC_ind = c(rep(c(1,3),times=6) + rep(3*(0:5),each=2))
Thetamin = -3
Thetamax = 3
Alphamax = 12
Sigma2max = 12
load("Mean_Cor_Real_Data.Rdata")

no_cores = detectCores(TRUE,TRUE)-1
cl=makeCluster(no_cores,type="FORK")

#Simulating Neochild Data
clusterSetRNGStream(cl,iseed=NULL)
Neochild.Sims <- clusterApply(cl,1:2, function(i){
  while(1){
    TimePC <- runif(1,0.001,1.15)
    TimePCA <- runif(1,0.1,15)
    Alphas = 10**runif(6,-Alphamax,Alphamax)
    Sigma2 = 10**runif(1,-Sigma2max,Sigma2max)
    while(1){
      #T1
      Theta1 = 10**runif(1,Thetamin,Thetamax)
      #T2
      Theta2 = 10**runif(1,Thetamin,Thetamax)
      if(abs(Theta1-Theta2)>eps){break}
    }
    sim.data = list("P" = Theta1,"C"=Theta2,"Alphas"=Alphas,"Sigma2"=Sigma2,"T.PC"=TimePC,"T.PCA"=TimePCA)
    x = Obs.neochild(sim.data)
    if(all(is.finite(x))){
      neoc.logll = x[length(x)]
      x = x[-length(x)]
      Obs = list("X_carcass"=x[1:3],"X_fHead"=x[4:6],"X_ovary"=x[7:9],"X_mHead"=x[10:12],"X_testis"=x[13:15],"X_accg"=x[16:18])
      Obs.data = list("Obs"=Obs,"T.PC"=TimePC,"T.PCA"=TimePCA)
      cons.x = c(Theta1,Alphas,Sigma2)
      neop.x = c(Theta2,Theta1,Alphas,Sigma2)
      cons.logll = logll.cons(x=cons.x,other.data=Obs.data)
      neop.logll = logll.neoparent(x=neop.x,other.data=Obs.data)
      if(all(is.finite(c(cons.logll,neop.logll)))){
        LR_neoc_v_cons = 2*(neoc.logll - cons.logll)
        LR_neop_v_cons = 2*(neop.logll - cons.logll)
        P_neoc_v_cons = pchisq(LR_neoc_v_cons,df=1,lower.tail=F)
        P_neop_v_cons = pchisq(LR_neop_v_cons,df=1,lower.tail=F)
        if(P_neoc_v_cons<sig/2 & P_neoc_v_cons>0 & P_neop_v_cons>sig/2 ) break
      }

    }
  }
  return(c("X"=x,"TPC"=TimePC,"TPCA"=TimePCA,"ThetaP"=Theta1,"ThetaC"=Theta2,"ThetaA"=Theta1,"Alphas"=Alphas,"Sigma2"=Sigma2,"loglikelihood"=neoc.logll))
})
Neoc.Sims = t(sapply(1:length(Neochild.Sims), function(x) Neochild.Sims[[x]]))
t2 <- Sys.time()
cat("Neochild Simulation done in ")
print(t2-t1)
t3 <- Sys.time()

# #Simulating Spec data
clusterSetRNGStream(cl,iseed=NULL)
Specialization.Sims <- clusterApply(cl,1:1, function(i){
  while(1){
    TimePC <- runif(1,0.001,1.15)
    TimePCA <- runif(1,0.1,15)
    Alphas = 10**runif(6,0,Alphamax)
    Sigma2 = 10**runif(1,0,Sigma2max)
    while(1){
      Theta1 = 10**runif(1,Thetamin,Thetamax)
      Theta2 = 10**runif(1,Thetamin,Thetamax)
      Theta3 = 10**runif(1,Thetamin,Thetamax)
      if(abs(Theta1-Theta2)>eps & abs(Theta2-Theta3)>eps & abs(Theta1 - Theta3)>eps){break}
    }
    sim.data = list("P" = Theta2,"C"=Theta3,"A"=Theta1,"Alphas"=Alphas,"Sigma2"=Sigma2,"T.PC"=TimePC,"T.PCA"=TimePCA)
    x = Obs.spec(sim.data)
    if(all(is.finite(x))){
      spec.logll = x[length(x)]
      x=x[-length(x)]
      Obs = list("X_carcass"=x[1:3],"X_fHead"=x[4:6],"X_ovary"=x[7:9],"X_mHead"=x[10:12],"X_testis"=x[13:15],"X_accg"=x[16:18])
      Obs.data = list("Obs"=Obs,"T.PC"=TimePC,"T.PCA"=TimePCA)
      cons.x = c(Theta1,Alphas,Sigma2)
      neoc.x = c(Theta1,Theta3,Alphas,Sigma2)
      neop.x = c(Theta2,Theta1,Alphas,Sigma2)
      sub.x = c(Theta2,Theta3,Alphas,Sigma2)
      cons.logll = logll.cons(x=cons.x,other.data=Obs.data)
      neoc.logll = logll.neochild(x=neoc.x,other.data=Obs.data)
      neop.logll = logll.neoparent(x=neop.x,other.data=Obs.data)
      sub.logll = logll.sub(x=sub.x,other.data=Obs.data)
      if(all(is.finite(c(cons.logll,neoc.logll,neop.logll,sub.logll)))){
        LR_spec_v_neoc = 2*(spec.logll - neoc.logll)
        LR_spec_v_neop = 2*(spec.logll - neop.logll)
        LR_spec_v_cons = 2*(spec.logll - cons.logll)
        LR_spec_v_sub = 2*(spec.logll - sub.logll)
        P_spec_v_neoc = pchisq(LR_spec_v_neoc,df=1,lower.tail=F)
        P_spec_v_neop = pchisq(LR_spec_v_neop,df=1,lower.tail=F)
        P_spec_v_cons = pchisq(LR_spec_v_cons,df=2,lower.tail=F)
        P_spec_v_sub = pchisq(LR_spec_v_sub,df=1,lower.tail=F)
        if(P_spec_v_neoc<sig/3 & P_spec_v_neop<sig/3 & P_spec_v_cons<sig/3 & P_spec_v_sub<sig2 &P_spec_v_neoc>0 &P_spec_v_neop>0 & P_spec_v_cons>0 & P_spec_v_sub>0 ) {
          break}
      }
    }
  }
  return(c("X"=x,"TPC"=TimePC,"TPCA"=TimePCA,"ThetaP"=Theta2,"ThetaC"=Theta3,"ThetaA"=Theta1,"Alphas"=Alphas,"Sigma2"=Sigma2,"loglikelihood"=spec.logll))
})
Spec.Sims = t(sapply(1:length(Specialization.Sims), function(x) Specialization.Sims[[x]]))
t4 <- Sys.time()
cat("Specialization Simulation done in ")
print(t4-t3)

# #Simulating Subfunc data
# clusterSetRNGStream(cl,iseed=NULL)
Subfunctionalization.Sims <- clusterApply(cl,1:1, function(i){
  while(1){
    TimePC <- runif(1,0.001,1.15)
    TimePCA <- runif(1,0.1,15)
    Alphas = 10**runif(6,-Alphamax,Alphamax)
    Sigma2 = 10**runif(1,-Sigma2max,Sigma2max)
    while(1){
      Theta1 = 10**runif(1,Thetamin,Thetamax)
      Theta2 = 10**runif(1,Thetamin,Thetamax)
      if(abs(Theta1-Theta2)>eps & (Theta1+Theta2)<=10**Thetamax){break}
    }
    sim.data = list("P" = Theta1,"C"=Theta2,"A"=Theta1+Theta2,"Alphas"=Alphas,"Sigma2"=Sigma2,"T.PC"=TimePC,"T.PCA"=TimePCA)
    x = Obs.spec(sim.data)
    if(all(is.finite(x))){
      sub.logll = x[length(x)]
      x=x[-length(x)]
      Obs = list("X_carcass"=x[1:3],"X_fHead"=x[4:6],"X_ovary"=x[7:9],"X_mHead"=x[10:12],"X_testis"=x[13:15],"X_accg"=x[16:18])
      Obs.data = list("Obs"=Obs,"T.PC"=TimePC,"T.PCA"=TimePCA)
      get.theta.data = list("Obs"=x,"T.PC"=TimePC,"T.PCA"=TimePCA,"params"=c(Theta1,Theta2,Alphas,Sigma2))
      cons.x = c(Theta1+Theta2,Alphas,Sigma2)
      neoc.x = c(Theta1+Theta2,Theta2,Alphas,Sigma2)
      neop.x = c(Theta1,Theta1+Theta2,Alphas,Sigma2)
      cons.logll = logll.cons(x=cons.x,other.data=Obs.data)
      neoc.logll = logll.neochild(x=neoc.x,other.data=Obs.data)
      neop.logll = logll.neoparent(x=neop.x,other.data=Obs.data)

      spec.runs <- lapply(1:100, function(j){
        while(1){
          spec.init <- 10**runif(1,Thetamin,Thetamax)
          # spec.init <- runif(1,sum(get.theta.data$params[1:2])-1,sum(get.theta.data$params[1:2])+1)
          if (!is.na(get.spec.theta(spec.init,get.theta.data))){break}
        }
        optim(spec.init, fn=get.spec.theta,other.data=get.theta.data, method="Brent",lower=10**Thetamin,upper=10**Thetamax,control=list("fnscale"=-1))
      })
      spec.runs=matrix(sapply(1:length(spec.runs), function(x) unlist(spec.runs[[x]][1:2])),ncol=length(spec.runs[[1]][1:2]), byrow=T)
      spec.ind = which.max(spec.runs[,ncol(spec.runs)])
      spec.logll = spec.runs[spec.ind,ncol(spec.runs)]
      if(all(is.finite(c(cons.logll,neoc.logll,neop.logll,spec.logll)))){
        LR_spec_v_neoc = 2*(spec.logll - neoc.logll)
        LR_spec_v_neop = 2*(spec.logll - neop.logll)
        LR_spec_v_cons = 2*(spec.logll - cons.logll)
        LR_spec_v_sub = 2*(spec.logll - sub.logll)
        P_spec_v_neoc = pchisq(LR_spec_v_neoc,df=1,lower.tail=F)
        P_spec_v_neop = pchisq(LR_spec_v_neop,df=1,lower.tail=F)
        P_spec_v_cons = pchisq(LR_spec_v_cons,df=2,lower.tail=F)
        P_spec_v_sub = pchisq(LR_spec_v_sub,df=1,lower.tail=F)
        if(P_spec_v_neoc<sig/3 & P_spec_v_neop<sig/3 & P_spec_v_cons<sig/3 & P_spec_v_sub>sig2 &P_spec_v_neoc>0 &P_spec_v_neop>0 & P_spec_v_cons>0 ) {
          break}
      }
    }
  }
  return(c("X"=x,"TPC"=TimePC,"TPCA"=TimePCA,"ThetaP"=Theta1,"ThetaC"=Theta2,"ThetaA"=spec.runs[spec.ind,1],"Alphas"=Alphas,"Sigma2"=Sigma2,"log.likelihood"=sub.logll))
})
Sub.Sims = t(sapply(1:length(Subfunctionalization.Sims), function(x) Subfunctionalization.Sims[[x]]))
t5 <- Sys.time()
# cat("Subfunctionalization Simulation done in ")
# print(t5-t4)

#Simulating Conservation data
clusterSetRNGStream(cl,iseed=NULL)
Conservation.Sims <- clusterApply(cl,1:1, function(i){
  while(1){
    TimePC <- runif(1,0.001,1.15)
    TimePCA <- runif(1,0.1,1.15)
    Alphas = 10**runif(6,-Alphamax,Alphamax)
    Sigma2 = 10**runif(1,-Sigma2max,Sigma2max)
    Theta1 = 10**runif(1,Thetamin,Thetamax)
    sim.data = list("P" = Theta1,"Alphas"=Alphas,"Sigma2"=Sigma2,"T.PC"=TimePC,"T.PCA"=TimePCA)
    x = Obs.cons(sim.data)
    if(all(is.finite(x))){
      cons.logll = x[length(x)]
      x = x[-length(x)]
      break
    }
  }
  return(c("X"=x,"TPC"=TimePC,"TPCA"=TimePCA,"ThetaP"=Theta1,"ThetaC"=Theta1,"ThetaA"=Theta1,"Alphas"=Alphas,"Sigma2"=Sigma2,"log.likelihood"=cons.logll))
})
Cons.Sims = t(sapply(1:length(Conservation.Sims), function(x) Conservation.Sims[[x]]))
t6 <- Sys.time()
cat("Conservation Simulation done in ")
print(t6-t5)
print("Done")
X = rbind(Neoc.Sims,Spec.Sims,Sub.Sims,Cons.Sims)
TimesPC = X[,19]
TimesPCA = X[,20]
Sim.ThetaP = X[,21]
Sim.ThetaC = X[,22]
Sim.ThetaA = X[,23]
Sim.Alpha = X[,24:29]
Sim.Sigma2 = X[,30]
Sim.logll = X[,31]

X_all <- list("carcass"=X[,1:3],"fHead"=X[,4:6],"ovary"=X[,7:9],"mHead"=X[,10:12],"testis"=X[,13:15],"accg"=X[,16:18])

Sim.Class = c(rep("Neochild",100),rep("Spec",100),rep("Sub",100),rep("Cons",100))
Sim.Class = factor(Sim.Class,levels=c("Cons","Neochild","Neoparent","Sub","Spec"))
Simulated.Models = list("ThetaP"=Sim.ThetaP,"ThetaC"=Sim.ThetaC,"ThetaA"=Sim.ThetaA,"TimesPC"=TimesPC,"TimesPCA"=TimesPCA,"Alpha"=Sim.Alpha,"Sigma2"=Sim.Sigma2,"X_all"=X_all,"Sim.Class"=Sim.Class,"logll"=Sim.logll)
stopCluster(cl)

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
cons.lower <- c(0,rep(eps,2))
neochild.lower <- c(0,0,rep(eps,2))
neoparent.lower <- c(0,0,rep(eps,2))
sub.lower <- c(0,0,rep(eps,2))
spec.lower <- c(0,0,0,rep(eps,2))
up.a <- 10**12
up.theta <- log10(range(Simulated.Models$X_all)[2]+1)
low.theta <- log10(range(Simulated.Models$X_all)[1]+1)

cons.upper <- c(up.a,rep(up.a,7))
neochild.upper <- c(up.a,up.a,rep(up.a,7))
neoparent.upper <- c(up.a,up.a,rep(up.a,7))
neoparent.upper <- c(up.a,up.a,rep(up.a,7))
spec.upper <- c(up.a,up.a,up.a,rep(up.a,7))
n <- 100
Pred.Class = rep("Cons",length(indices))
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
# varlist2 =  c("i","neochild.choice","neoparent.choice","cons.choice","spec.choice")
varlist2 =  c("i","neochild.choice","cons.choice","spec.choice","neoparent.choice","sub.choice")
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
    Pred.Class[i]= "Neochild"}
  if(Neochild.Model$Pvalues_cons[i]>=sig & Neoparent.Model$Pvalues_cons[i]<sig & Neoparent.Model$Pvalues_cons[i]>0 & Spec.Model$Pvalues_neop[i]>=sig){
    Pred.Class[i]= "Neoparent"}
  if(Spec.Model$Pvalues_neoc[i]<sig & Spec.Model$Pvalues_neop[i]<sig & Spec.Model$Pvalues_cons[i]<sig & Spec.Model$Pvalues_neoc[i]>0 & Spec.Model$Pvalues_neop[i]>0 & Spec.Model$Pvalues_cons[i]>0){
    Pred.Class[i]= "Spec"
    
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
    if(Pval.spec_v_sub>sig2){Pred.Class[i]= "Sub"}
  }
  print(i)
}
stopCluster(cl)
gc()

# library(gridExtra)
Pred.Class = factor(Pred.Class,levels=c("Cons","Neochild","Neoparent","Sub","Spec"))
print(table(Pred.Class,Sim.Class))
save(Pred.Class,Simulated.Models,Neochild.Model,Neoparent.Model,Cons.Model, Spec.Model,Sub.Model, file="Sim_Models_Aallcor.Rdata")
t2 <- Sys.time()
print(t2 - t1)
cat("Exited Successfully!","\n")