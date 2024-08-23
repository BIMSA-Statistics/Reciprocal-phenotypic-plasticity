library(mvtnorm)
get_simiu_data <- function(par=matrix(c(1,1,0.13,0.08,0.06,0.1,0.03,0.05),
                                      nrow = 2),n=100,he=0.05,hs=0.05){
  
  
  miu_e <- par[1,1]
  a_e_e <- par[1,2]
  a_e_s <- par[1,3]
  i_e <- par[1,4]
  
  
  miu_s <- par[2,1]
  a_s_e <- par[2,2]
  a_s_s <- par[2,3]
  i_s <- par[2,4]
  
  
  
  
  miu11 <- c(miu_e+a_e_e+a_e_s+i_e,miu_s+a_s_e+a_s_s+i_s)
  miu12 <- c(miu_e+a_e_e-a_e_s-i_e,miu_s+a_s_e-a_s_s-i_s)
  miu21 <- c(miu_e-a_e_e+a_e_s-i_e,miu_s-a_s_e+a_s_s-i_s)
  miu22 <- c(miu_e-a_e_e-a_e_s+i_e,miu_s-a_s_e-a_s_s+i_s)
  
  
  sigma_ge <- var(c(miu11[1],miu12[1],miu21[1],miu22[1]))
  
  sigma_gs <- var(c(miu11[2],miu12[2],miu21[2],miu22[2]))
  
  
  sigma_e <- (sigma_ge-sigma_ge*he)/he
  
  sigma_s <- (sigma_gs-sigma_gs*hs)/hs
  
  sigma <- matrix(c(sigma_e,0.5*sqrt(sigma_e*sigma_s),0.5*sqrt(sigma_e*sigma_s),sigma_s),2,2)
  
  simiu_data <- matrix(NA,nrow=n,ncol=4)
  for (i in 1:(n/4)) {
    
    simiu_data[(i-1)*4+1,] <- c(rmvnorm(1,miu11,sigma),1,1)
    
    simiu_data[(i-1)*4+2,] <- c(rmvnorm(1,miu12,sigma),1,2)
    
    simiu_data[(i-1)*4+3,] <- c(rmvnorm(1,miu21,sigma),2,1)
    
    simiu_data[(i-1)*4+4,] <- c(rmvnorm(1,miu22,sigma),2,2)}
  
  colnames(simiu_data) <- c("miu_e","miu_s","geno_e","geno_s")
  return((simiu_data))
}

simiu_data <- get_simiu_data()

get_par_res <- function(simiu_data ){
  
  par11 <-colMeans( simiu_data[which(simiu_data[,3]==1&simiu_data[,4]==1),1:2])
  par12 <-colMeans( simiu_data[which(simiu_data[,3]==1&simiu_data[,4]==2),1:2])
  par21 <-colMeans( simiu_data[which(simiu_data[,3]==2&simiu_data[,4]==1),1:2])
  par22 <-colMeans( simiu_data[which(simiu_data[,3]==2&simiu_data[,4]==2),1:2])
  
  
  
  LL <- optim(par = c(0.1,0.08,0.5,par11,par12,par21,par22),fn = LL2,
              yt11=simiu_data[which(simiu_data[,3]==1&simiu_data[,4]==1),1:2],
              yt12=simiu_data[which(simiu_data[,3]==1&simiu_data[,4]==2),1:2],
              yt21=simiu_data[which(simiu_data[,3]==2&simiu_data[,4]==1),1:2],
              yt22=simiu_data[which(simiu_data[,3]==2&simiu_data[,4]==2),1:2],
              method = "BFGS"#lower = rep(0,11),upper = rep(Inf,11),method = "L-BFGS-B"
  )
  
  
  miu11 <- LL$par[4:5]
  miu12 <- LL$par[6:7]
  miu21 <- LL$par[8:9]
  miu22 <- LL$par[10:11]
  
  a_e_e <- 1/4*(miu11[1]+miu12[1]-miu21[1]-miu22[1])
  a_s_s <- 1/4*(miu11[2]+miu21[2]-miu12[2]-miu22[2])
  
  a_s_e <- 1/4*(miu11[2]+miu12[2]-miu21[2]-miu22[2])
  
  a_e_s <- 1/4*(miu11[1]+miu21[1]-miu12[1]-miu22[1])
  
  i_e <- 1/4*(miu11[1]+miu22[1]-miu12[1]-miu21[1]) 
  i_s <- 1/4*(miu11[2]+miu22[2]-miu12[2]-miu21[2]) 
  miu_e <- 1/4*(miu11[1]+miu12[1]+miu21[1]+miu22[1])
  miu_s <- 1/4*(miu11[2]+miu12[2]+miu21[2]+miu22[2])
  
  res <- c(miu_e,a_e_e, a_e_s,i_e, miu_s,a_s_e, a_s_s,i_s)
  names(res) <- c("miu_e", "a_e_e", "a_e_s","i_e", "miu_s", "a_s_e", "a_s_s","i_s")
  return(res)
}




res <- matrix(NA,nrow = 1000,ncol = 8)
for (i in 1:1000) {
  simiu_data <- get_simiu_data()
  res[i,] <- get_par_res(simiu_data)
  
  
}


rres <-    rbind(colMeans(res),
                 sapply(1:8,function(c)sd(res[,c])))

write.csv(rres,file = "同遗传力0.05optim_n100.csv")


res <- matrix(NA,nrow = 1000,ncol = 8)
for (i in 1:1000) {
  simiu_data <- get_simiu_data(n=500)
  res[i,] <- get_par_res(simiu_data)
  
  
}


rres <-    rbind(colMeans(res),
                 sapply(1:8,function(c)sd(res[,c])))

write.csv(rres,file = "同遗传力0.05optim_n500.csv")

get_LR <- function(simiu_data ){
  par11 <-colMeans( simiu_data[which(simiu_data[,3]==1&simiu_data[,4]==1),1:2])
  par12 <-colMeans( simiu_data[which(simiu_data[,3]==1&simiu_data[,4]==2),1:2])
  par21 <-colMeans( simiu_data[which(simiu_data[,3]==2&simiu_data[,4]==1),1:2])
  par22 <-colMeans( simiu_data[which(simiu_data[,3]==2&simiu_data[,4]==2),1:2])
  
  
  LL <- optim(par = c(0.2,0.15,0.5,par11,par12,par21,par22),fn = LL2,
              yt11=simiu_data[which(simiu_data[,3]==1&simiu_data[,4]==1),1:2],
              yt12=simiu_data[which(simiu_data[,3]==1&simiu_data[,4]==2),1:2],
              yt21=simiu_data[which(simiu_data[,3]==2&simiu_data[,4]==1),1:2],
              yt22=simiu_data[which(simiu_data[,3]==2&simiu_data[,4]==2),1:2],
              method = "BFGS"#lower = rep(0,11),upper = rep(Inf,11),method = "L-BFGS-B"
  )
  
  L1=-LL$value
  
  par <- colMeans( simiu_data[,1:2])
  
  ll <- optim(par =c(0.2,0.15,0.5,par),fn = LL0,
              yt=simiu_data[,1:2],
              method = "BFGS"#lower = rep(0,11),upper = rep(Inf,11),method = "L-BFGS-B"
  )
  
  
  L0 <- -ll$value
  LR=  -2*(L0-L1)
  
  
  return(LR)
}
get_LR2 <- function(simiu_data ){
  ######H1 
  par11 <-colMeans( simiu_data[which(simiu_data[,3]==1&simiu_data[,4]==1),1:2])
  par12 <-colMeans( simiu_data[which(simiu_data[,3]==1&simiu_data[,4]==2),1:2])
  par21 <-colMeans( simiu_data[which(simiu_data[,3]==2&simiu_data[,4]==1),1:2])
  par22 <-colMeans( simiu_data[which(simiu_data[,3]==2&simiu_data[,4]==2),1:2])
  
  
  
  LL <- optim(par = c(0.2,0.15,0.5,par11,par12,par21,par22),fn = LL2,
              yt11=simiu_data[which(simiu_data[,3]==1&simiu_data[,4]==1),1:2],
              yt12=simiu_data[which(simiu_data[,3]==1&simiu_data[,4]==2),1:2],
              yt21=simiu_data[which(simiu_data[,3]==2&simiu_data[,4]==1),1:2],
              yt22=simiu_data[which(simiu_data[,3]==2&simiu_data[,4]==2),1:2],
              method = "BFGS"#lower = rep(0,11),upper = rep(Inf,11),method = "L-BFGS-B"
  )
  
  L1=-LL$value
  ###########H0
  miu1 <- c(mean(simiu_data[which(simiu_data[,3]==1),1]),
            mean(simiu_data[which(simiu_data[,4]==1),2]))
  
  miu2 <- c(mean(simiu_data[which(simiu_data[,3]==2),1]),
            mean(simiu_data[which(simiu_data[,4]==2),2]))
  
  ll <-  optim(par = c(0.2,0.15,0.5,miu1,miu2),fn = LL1,
               yt11=simiu_data[which(simiu_data[,3]==1&simiu_data[,4]==1),1:2],
               yt12=simiu_data[which(simiu_data[,3]==1&simiu_data[,4]==2),1:2],
               yt21=simiu_data[which(simiu_data[,3]==2&simiu_data[,4]==1),1:2],
               yt22=simiu_data[which(simiu_data[,3]==2&simiu_data[,4]==2),1:2],
               method = "BFGS"#lower = rep(0,11),upper = rep(Inf,11),method = "L-BFGS-B"
  )
  
  
  L0 <- -ll$value
  LR=  -2*(L0-L1)
  
  
  return(LR)
}

get_power <- function(par=matrix(c(1,1,0.13,0.08,0.06,0.1,0.03,0.05),
                                 nrow = 2),n=100,he=0.05,hs=0.05){
  res <- c()
  
  
  simiu_data <- get_simiu_data(par,n,he=0.05,hs=0.05)
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


get_power2 <- function(par=matrix(c(1,1,0.13,0.08,0.06,0.1,0.03,0.05),
                                  nrow = 2),n=100,he=0.05,hs=0.05){
  res <- c()
  
  
  simiu_data <- get_simiu_data(par,n,he=0.05,hs=0.05)
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
res <- pbsapply(1:1000,function(c)get_power(par=matrix(c(1,1,0.13,0.08,0.06,0.1,0.03,0.05),
                                                      nrow = 2),n=100,he=0.05,hs=0.05),cl=cl)

stopCluster(cl)

power3 <- length(which(res[1,]>res[2,]))/1000



library(parallel)
library(pbapply)
cl.cores <- 10
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("get_simiu_data",'get_LR','get_power','get_LR2','get_power2',"LL2","LL1","LL0","get_covmatrix"))
clusterEvalQ(cl,{library(mvtnorm)})
res <- pbsapply(1:1000,function(c)get_power(par=matrix(c(1,1,0.13,0.08,0.06,0.1,0.03,0.05),
                                                      nrow = 2),n=500,he=0.05,hs=0.05),cl=cl)

stopCluster(cl)

power3_n500 <- length(which(res[1,]>res[2,]))/1000

library(parallel)
library(pbapply)
cl.cores <- 10
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("get_simiu_data",'get_LR','get_power','get_LR2','get_power2',"LL2","LL1","LL0","get_covmatrix"))
clusterEvalQ(cl,{library(mvtnorm)})
res <- pbsapply(1:1000,function(c)get_power2(par=matrix(c(1,1,0.13,0.08,0.06,0.1,0.03,0.05),
                                                       nrow = 2),n=100,he=0.05,hs=0.05),cl=cl)

stopCluster(cl)

tpower3 <- length(which(res[1,]>res[2,]))/1000

library(parallel)
library(pbapply)
cl.cores <- 10
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("get_simiu_data",'get_LR','get_power','get_LR2','get_power2',"LL2","LL1","LL0","get_covmatrix"))
clusterEvalQ(cl,{library(mvtnorm)})
res <- pbsapply(1:1000,function(c)get_power2(par=matrix(c(1,1,0.13,0.08,0.06,0.1,0.03,0.05),
                                                       nrow = 2),n=500,he=0.05,hs=0.05),cl=cl)

stopCluster(cl)

tpower3_n500 <- length(which(res[1,]>res[2,]))/1000

# 
# library(parallel)
# library(pbapply)
# cl.cores <- 10
# cl <- makeCluster(getOption("cl.cores", cl.cores ))
# clusterExport(cl,c("get_simiu_data",'get_LR','get_power','get_LR2','get_power2',"LL2","LL1","LL0","get_covmatrix"))
# clusterEvalQ(cl,{library(mvtnorm)})
# res <- pbsapply(1:100,function(c)get_power(par=matrix(c(1,1,0,0,0,0,0,0),
#                                                       nrow = 2),n=100,he=0.1,hs=0.05),cl=cl)
# 
# stopCluster(cl)
# 
# false3 <- length(which(res[1,]>res[2,]))/100
# 
# 
# library(parallel)
# library(pbapply)
# cl.cores <- 10
# cl <- makeCluster(getOption("cl.cores", cl.cores ))
# clusterExport(cl,c("get_simiu_data",'get_LR','get_power','get_LR2','get_power2',"LL2","LL1","LL0","get_covmatrix"))
# clusterEvalQ(cl,{library(mvtnorm)})
# res <- pbsapply(1:100,function(c)get_power(par=matrix(c(1,1,0,0,0,0,0,0),
#                                                       nrow = 2),n=500,he=0.1,hs=0.05),cl=cl)
# 
# stopCluster(cl)
# 
# false3_n500 <- length(which(res[1,]>res[2,]))/100


library(parallel)
library(pbapply)
cl.cores <- 10
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("get_simiu_data",'get_LR','get_power','get_LR2','get_power2',"LL2","LL1","LL0","get_covmatrix"))
clusterEvalQ(cl,{library(mvtnorm)})
res <- pbsapply(1:1000,function(c)get_power2(par=matrix(c(1,1,0.13,0,0,0.1,0,0),
                                                       nrow = 2),n=100,he=0.05,hs=0.05),cl=cl)

stopCluster(cl)

tfalse3 <- length(which(res[1,]>res[2,]))/1000


library(parallel)
library(pbapply)
cl.cores <- 10
cl <- makeCluster(getOption("cl.cores", cl.cores ))
clusterExport(cl,c("get_simiu_data",'get_LR','get_power','get_LR2','get_power2',"LL2","LL1","LL0","get_covmatrix"))
clusterEvalQ(cl,{library(mvtnorm)})
res <- pbsapply(1:1000,function(c)get_power2(par=matrix(c(1,1,0.13,0,0,0.1,0,0),
                                                       nrow = 2),n=500,he=0.05,hs=0.05),cl=cl)

stopCluster(cl)

tfalse3_n500 <- length(which(res[1,]>res[2,]))/1000








