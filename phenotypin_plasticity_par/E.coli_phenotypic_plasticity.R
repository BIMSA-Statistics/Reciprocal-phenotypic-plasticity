rm(list=ls())
setwd("R:\\Rcode\\上传脚本\\phenotypin_plasticity_par")
#-----------------------------------requiredPackages----------------------------
requiredPackages = c("mvtnorm","parallel","pbapply")
for(packages in requiredPackages){
  if(!require(packages,character.only = TRUE)) install.packages(packages)
  require(packages,character.only = TRUE)
}
#-----------------------------------inputdata-----------------------------------
E_mo <- read.csv("E_mo.csv",row.names = 1)
S_mo <- read.csv("S_mo.csv",row.names = 1)
ES_E <- read.csv("ES_E_co.csv",row.names = 1)
ES_S <- read.csv("ES_S_co.csv",row.names = 1)
S_SNP <- read.table("S-SNP.txt",row.names = 1)
E_SNP <- read.table("E-SNP.txt",row.names = 1)
E_mo <- log(E_mo)
S_mo <- log(S_mo)
ES_E <- log(ES_E)
ES_S <- log(ES_S)
X_E <-E_mo-ES_E
#-----------------------------------function------------------------------------
get_miu<-function(par, t, options=list())
{
  y <- par[1]/(1+par[2]*exp(-par[3]*t))-par[4]/(1+par[5]*exp(-par[6]*t))
  return (y);
}
#-----------------------------------SAD-----------------------------------------
SAD1_get_matrix = function(par, times =Time, options=list()) {   
  n <- ifelse (is.vector(times), length(times), NCOL(times) )   
  phi<- par[1]   
  v2 <- par[2]
  sigma <- array(0,dim = c(n,n))
  sigma1 <- array( 0, dim = c(n,n))
  
  
  for(i in 1:n)   
  {
    sigma1[i:n,i] <- phi^( times[i:n] - times[i] )     
  }
  
  sigma <- sigma1%*%t(sigma1)*abs(v2)
  return(sigma = sigma); 
}
#-----------------------------------H0------------------------------------------
t=c(0.5,1,1.5,2,4,6,8,10,12,16,20)
par=c(26.4691709,102.7073464,1.0095179,27.8153760,70.2178921,0.9077182,0.8961068,2.7064985)
H0 = function(yt,t,par){
  miu = get_miu(par=c(par[1],par[2],par[3],par[4],par[5],par[6]),t)
  sigma = SAD1_get_matrix(par = c(par[7],par[8]), times = t )
  L0 = c()
  L0 = sum(dmvnorm(yt,miu,sigma,log = T))
  return(-L0)
}
H0_par<-optim(par=par,H0,yt=X_E,t=t,method="Nelder-Mead",control=list(maxit=200000))
H0_par
sigma= SAD1_get_matrix(par = c(H0_par$par[7:8]), times = t )
#-----------------------------------H1------------------------------------------
H1_lr <- function(s,eall){
  tem_SNP <- t(sapply(1:nrow(eall),function(c)paste0(eall[c,],s)))
  rownames(tem_SNP) <- rownames(eall)
  colnames(tem_SNP) <- colnames(eall)
  fn<-function(k){ 
    x11<-which(k=="11")
    x10<-which(k=="10")
    x00<-which(k=="00")
    x01<-which(k=="01")
    H1 = function(par,t){
      miu11= get_miu(par=c(par[1],par[2],par[3],par[4],par[5],par[6]),t)
      sigma = sigma
      L11=c()
      L11 = sum(dmvnorm(X_E[x11,],miu11,sigma,log = T))
      return(-L11)
    }
    H2 = function(par,t){
      miu10= get_miu(par=c(par[1],par[2],par[3],par[4],par[5],par[6]),t)
      sigma= sigma
      L10=c()
      L10 = sum(dmvnorm(X_E[x10,],miu10,sigma,log = T))
      return(-L10)
    }
    H3 = function(par,t){
      miu00= get_miu(par=c(par[1],par[2],par[3],par[4],par[5],par[6]),t)
      sigma= sigma
      L00=c()
      L00 = sum(dmvnorm(X_E[x00,],miu00,sigma,log = T))
      return(-L00)
    }
    H4 = function(par,t){
      miu01= get_miu(par=c(par[1],par[2],par[3],par[4],par[5],par[6]),t)
      sigma= sigma
      L01=c()
      L01 = sum(dmvnorm(X_E[x01,],miu01,sigma,log = T))
      return(-L01)
    }
    
    a1<-optim(par=c(H0_par$par[1:6]),t=t,H1,method="BFGS",control=list(maxit=1000))
    a2<-optim(par=c(H0_par$par[1:6]),t=t,H2,method="BFGS",control=list(maxit=1000))
    a3<-optim(par=c(H0_par$par[1:6]),t=t,H3,method="BFGS",control=list(maxit=1000))
    a4<-optim(par=c(H0_par$par[1:6]),t=t,H4,method="BFGS",control=list(maxit=1000))
    
    Lr<-(-2*((a1$value+a2$value+a3$value+a4$value)-H0_par$value))
    
    par <- c(Lr,a1$par,a2$par,a3$par,a4$par)
    return(par)
  }
  
  core.number <- detectCores()
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterExport(cl, c("dmvnorm","t","get_miu","X_E","sigma","H0_par"),envir = environment())
  aaa <- pbapply(tem_SNP,1,fn,cl=cl)
  stopCluster(cl)
  return(t(aaa))
}
#-----------------------------------one time combination------------------------
system.time(cc<- H1_lr(S_SNP[1,],E_SNP))
#-----------------------------------All S_SNP with E_SNP combination------------
all.cc <- list()
for(z in 1:dim(S_SNP)[1]){
  all.cc[[z]] <- H1_lr(S_SNP[z,],E_SNP)
}


