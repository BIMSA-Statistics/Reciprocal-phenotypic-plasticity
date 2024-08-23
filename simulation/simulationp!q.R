library(mvtnorm)
get_simiu_data <- function(par=matrix(c(1,1,0.13,0.08,0.06,0.1,0.03,0.05,0.2,0.15,0.5,0.5),
                                      nrow = 2),n=100,p=0.6,q=0.3){
  
  
  miu_e <- par[1,1]
  a_e_e <- par[1,2]
  a_e_s <- par[1,3]
  i_e <- par[1,4]
  sigma_e <- par[1,5]
  phi <- par[1,6]
  
  miu_s <- par[2,1]
  a_s_e <- par[2,2]
  a_s_s <- par[2,3]
  i_s <- par[2,4]
  sigma_s <- par[2,5]
  phi <- par[2,6]
  
  
  sigma <- matrix(c(sigma_e,phi*sqrt(sigma_e*sigma_s),phi*sqrt(sigma_e*sigma_s),sigma_s),2,2)
  
  miu11 <- c(miu_e+a_e_e+a_e_s+i_e,miu_s+a_s_e+a_s_s+i_s)
  miu12 <- c(miu_e+a_e_e-a_e_s-i_e,miu_s+a_s_e-a_s_s-i_s)
  miu21 <- c(miu_e-a_e_e+a_e_s-i_e,miu_s-a_s_e+a_s_s-i_s)
  miu22 <- c(miu_e-a_e_e-a_e_s+i_e,miu_s-a_s_e-a_s_s+i_s)
  simiu_data <- matrix(NA,nrow=n,ncol=4)
  # X=sample(x=c("11","12","21","22"),size=n,replace = T,prob=c(p*q,p*(1-q),(1-p)*q,(1-p)*(1-q)))
  # for (i in 1:n) {
  #  
  #   if( X[i]=="11"){
  #     simiu_data[i,] <- c(rmvnorm(1,miu11,sigma),1,1)
  #   }else if(X[i]=="12"){
  #     simiu_data[i,] <- c(rmvnorm(1,miu12,sigma),1,2)
  #   }else if(X[i]=="21"){
  #     simiu_data[i,] <- c(rmvnorm(1,miu21,sigma),2,1)
  #   }else{
  #     simiu_data[i,] <- c(rmvnorm(1,miu22,sigma),2,2)
  #   }}
  X=sample(x=c("11","12","21","22"),size=n,replace = T,prob=c(p*q,p*(1-q),(1-p)*q,(1-p)*(1-q)))
  for (i in 1:n) {
     
    if( 1<=i&i<=(n*p*q)){
      simiu_data[i,] <- c(rmvnorm(1,miu11,sigma),1,1)
    }else if((n*p*q+1)<=i&i<=(n*p*q+n*p*(1-q))){
      simiu_data[i,] <- c(rmvnorm(1,miu12,sigma),1,2)
    }else if((n*p*q+n*p*(1-q)+1)<=i&i<=(n*p*q+n*p*(1-q)+n*(1-p)*q)){
      simiu_data[i,] <- c(rmvnorm(1,miu21,sigma),2,1)
    }else{
      simiu_data[i,] <- c(rmvnorm(1,miu22,sigma),2,2)
    }}
    
 
       
  colnames(simiu_data) <- c("miu_e","miu_s","geno_e","geno_s")
  return((simiu_data))
}


simiu_data <- get_simiu_data(n=100)

length(which(simiu_data[,3]==1&simiu_data[,4]==2))




res <- matrix(NA,nrow = 1000,ncol = 8)
for (i in 1:1000) {
  simiu_data <- get_simiu_data()
  res[i,] <- get_par_res(simiu_data)
  
  
}


rres <-    rbind(colMeans(res),
                 sapply(1:8,function(c)sd(res[,c])))

write.csv(rres,file = "模拟不同频率optim_n100.csv")


res <- matrix(NA,nrow = 1000,ncol = 8)
for (i in 1:1000) {
  simiu_data <- get_simiu_data(n=500)
  res[i,] <- get_par_res(simiu_data)
  
  
}


rres <-    rbind(colMeans(res),
                 sapply(1:8,function(c)sd(res[,c])))

write.csv(rres,file = "模拟不同频率optim_n500.csv")

get_power <- function(par=matrix(c(1,1,0.13,0.08,0.06,0.1,0.03,0.05,0.2,0.15,0.5,0.5),
                                 nrow = 2),n=100,p=0.6,q=0.3){
  res <- c()
  
  
  simiu_data <- get_simiu_data(par,n,p=0.6,q=0.3)
  LR <- get_LR(simiu_data)
  LR_test <- c()
  for (j in 1:100) {
    simiu_data[,3:4] <- simiu_data[sample(1:nrow(simiu_data),nrow(simiu_data)),3:4]
    LR_test[j] <- get_LR(simiu_data)
  }
  
  
  thre <- sort(LR_test,decreasing = T)[5]
  res[1] <- LR
  res[2] <- thre
  return(res)
}


get_power2 <- function(par=matrix(c(1,1,0.13,0.08,0.06,0.1,0.03,0.05,0.2,0.15,0.5,0.5),
                                  nrow = 2),n=100,p=0.6,q=0.3){
  res <- c()
  
  
  simiu_data <- get_simiu_data(par,n,p=0.6,q=0.3)
  LR <- get_LR2(simiu_data)
  LR_test <- c()
  for (j in 1:100) {
    simiu_data[,3:4] <- simiu_data[sample(1:nrow(simiu_data),nrow(simiu_data)),3:4]
    LR_test[j] <- get_LR2(simiu_data)
  }
  
  
  thre <- sort(LR_test,decreasing = T)[5]
  res[1] <- LR
  res[2] <- thre
  return(res)
}


library(parallel)
library(pbapply)
cl.cores <- 10
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("get_simiu_data",'get_LR','get_power','get_LR2','get_power2',"LL2","LL1","LL0","get_covmatrix"))
clusterEvalQ(cl,{library(mvtnorm)})
res <- pbsapply(1:1000,function(c)get_power(par=matrix(c(1,1,0.13,0.08,0.06,0.1,0.03,0.05,0.2,0.15,0.5,0.5),
                                                      nrow = 2),n=100,p=0.6,q=0.3),cl=cl)

stopCluster(cl)

power2 <- length(which(res[1,]>res[2,]))/1000



library(parallel)
library(pbapply)
cl.cores <- 10
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("get_simiu_data",'get_LR','get_power','get_LR2','get_power2',"LL2","LL1","LL0","get_covmatrix"))
clusterEvalQ(cl,{library(mvtnorm)})
res <- pbsapply(1:1000,function(c)get_power(par=matrix(c(1,1,0.13,0.08,0.06,0.1,0.03,0.05,0.2,0.15,0.5,0.5),
                                                      nrow = 2),n=500,p=0.6,q=0.3),cl=cl)

stopCluster(cl)

power2_n500 <- length(which(res[1,]>res[2,]))/1000

library(parallel)
library(pbapply)
cl.cores <- 10
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("get_simiu_data",'get_LR','get_power','get_LR2','get_power2',"LL2","LL1","LL0","get_covmatrix"))
clusterEvalQ(cl,{library(mvtnorm)})
res <- pbsapply(1:1000,function(c)get_power2(par=matrix(c(1,1,0.13,0.08,0.06,0.1,0.03,0.05,0.2,0.15,0.5,0.5),
                                                       nrow = 2),n=100,p=0.6,q=0.3),cl=cl)

stopCluster(cl)

tpower2 <- length(which(res[1,]>res[2,]))/1000

library(parallel)
library(pbapply)
cl.cores <- 10
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("get_simiu_data",'get_LR','get_power','get_LR2','get_power2',"LL2","LL1","LL0","get_covmatrix"))
clusterEvalQ(cl,{library(mvtnorm)})
res <- pbsapply(1:1000,function(c)get_power2(par=matrix(c(1,1,0.13,0.08,0.06,0.1,0.03,0.05,0.2,0.15,0.5,0.5),
                                                       nrow = 2),n=500,p=0.6,q=0.3),cl=cl)

stopCluster(cl)

tpower2_n500 <- length(which(res[1,]>res[2,]))/1000


library(parallel)
library(pbapply)
cl.cores <- 10
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("get_simiu_data",'get_LR','get_power','get_LR2','get_power2',"LL2","LL1","LL0","get_covmatrix"))
clusterEvalQ(cl,{library(mvtnorm)})
res <- pbsapply(1:1000,function(c)get_power(par=matrix(c(1,1,0,0,0,0,0,0,0.2,0.15,0.5,0.5),
                                                      nrow = 2),n=100,p=0.6,q=0.3),cl=cl)

stopCluster(cl)

false2 <- length(which(res[1,]>res[2,]))/1000


library(parallel)
library(pbapply)
cl.cores <- 10
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("get_simiu_data",'get_LR','get_power','get_LR2','get_power2',"LL2","LL1","LL0","get_covmatrix"))
clusterEvalQ(cl,{library(mvtnorm)})
res <- pbsapply(1:1000,function(c)get_power(par=matrix(c(1,1,0,0,0,0,0,0,0.2,0.15,0.5,0.5),
                                                      nrow = 2),n=500,p=0.6,q=0.3),cl=cl)

stopCluster(cl)

false2_n500 <- length(which(res[1,]>res[2,]))/1000


library(parallel)
library(pbapply)
cl.cores <- 10
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("get_simiu_data",'get_LR','get_power','get_LR2','get_power2',"LL2","LL1","LL0","get_covmatrix"))
clusterEvalQ(cl,{library(mvtnorm)})
res <- pbsapply(1:1000,function(c)get_power2(par=matrix(c(1,1,0.13,0,0,0.1,0,0,0.2,0.15,0.5,0.5),
                                                           nrow = 2),n=100,p=0.6,q=0.3),cl=cl)

stopCluster(cl)

tfalse2 <- length(which(res[1,]>res[2,]))/1000


library(parallel)
library(pbapply)
cl.cores <- 10
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("get_simiu_data",'get_LR','get_power','get_LR2','get_power2',"LL2","LL1","LL0","get_covmatrix"))
clusterEvalQ(cl,{library(mvtnorm)})
res <- pbsapply(1:1000,function(c)get_power2(par=matrix(c(1,1,0.13,0,0,0.1,0,0,0.2,0.15,0.5,0.5),
                                                       nrow = 2),n=500,p=0.6,q=0.3),cl=cl)

stopCluster(cl)

tfalse2_n500 <- length(which(res[1,]>res[2,]))/1000

cl.cores <- 10
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("get_simiu_data",'get_LR','get_power','get_LR2','get_power2',"LL2","LL1","LL0","get_covmatrix"))
clusterEvalQ(cl,{library(mvtnorm)})
res <- pbsapply(1:1000,function(c)get_power2(par=matrix(c(1,1,0.13,0,0,0.1,0,0,0.2,0.15,0.5,0.5),
                                                       nrow = 2),n=1000,p=0.6,q=0.3),cl=cl)

tfalse2_n5002 <- length(which(res[1,]>res[2,]))/1000
cl.cores <- 10
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("get_simiu_data",'get_LR','get_power','get_LR2','get_power2',"LL2","LL1","LL0","get_covmatrix"))
clusterEvalQ(cl,{library(mvtnorm)})
res <- pbsapply(1:1000,function(c)get_power2(par=matrix(c(1,1,0.13,0,0,0.1,0,0,0.2,0.15,0.5,0.5),
                                                       nrow = 2),n=1000,p=0.6,q=0.3),cl=cl)

tfalse2_n5003 <- length(which(res[1,]>res[2,]))/1000




