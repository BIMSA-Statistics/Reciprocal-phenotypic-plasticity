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
LC_init_par_fn<-function(pheno_data, times){
  m <- dim(pheno_data)[2]
  um <- mean(pheno_data[,m], na.rm=T)
  um.1 <- mean(pheno_data[,m-1], na.rm=T)
  u2 <- mean(pheno_data[,2], na.rm=T)
  u1 <- mean(pheno_data[,1], na.rm=T)
  if (u2>u1) 
    u0 <- u1*3/4 
  else
    u0 <- u1*4/3
  par.a <- um * um / um.1
  par.b <- par.a / u0 -1
  if(par.b==Inf){
    par.b <- u2
  }
  par.r <- try( -log((par.a/um-1)/par.b)/m )
  if (is.na(par.r)|par.r==Inf|par.r==-Inf) par.r <- 0.6
  if(class(par.r)=="try-error")
    par.r <- 1
  par.A <- par.a
  par.R <- par.A*par.r/4
  par.lamda <- (par.A*(log(par.b)-2))/4/par.R
  return(c(par.A, par.R, par.lamda))
  # return(c(par.a, par.b, par.r))
}
SAD_init_par_fn<-function(pheno_data){
  m <- dim(pheno_data)[2]
  dat.t0<-pheno_data[ , 1] 
  par.s2 <- sd(dat.t0, na.rm=T)
  
  dat.t1<-pheno_data[ , 2 ]
  par.rho <- cor(dat.t1, dat.t0, use="complete.obs")
  par.phi <- par.rho/par.s2^2
  
  cat("cov.1",par.phi,par.s2,"\n")
  return( c( par.phi, par.s2))
}
get_miu<-function(par, t, options=list())
{
  y <- par[1]/(1+par[2]*exp(-par[3]*t))
  return (y);
}
