#-----------------------------------requiredPackages----------------------------
requiredPackages = c("mvtnorm","parallel","RColorBrewer","mvtnorm","scales","reshape2","pbapply","plyr","patchwork","ggplot2","deSolve","Rmisc","glmnet","gridExtra","tidyverse","ggpubr")
for(packages in requiredPackages){
  if(!require(packages,character.only = TRUE)) install.packages(packages)
  require(packages,character.only = TRUE)
}
#-----------------------------------inputdata-----------------------------------
t=c(0.5,1,1.5,2,4,6,8,10,12,16,20)
E_mo <- read.csv("E_mo.csv",row.names = 1,header=T)
S_mo <- read.csv("S_mo.csv",row.names = 1,header=T)
ES_E <- read.csv("ES_E_co.csv",row.names = 1,header=T)
ES_S <- read.csv("ES_s_co.csv",row.names = 1,header=T)
S_SNP <- read.table("S-SNP.txt",row.names = 1,header=T)
E_SNP <- read.table("E-SNP.txt",row.names = 1,header=T)
E_mo <- log10(E_mo)
S_mo <- log10(S_mo)
ES_E<-log10(ES_E)
ES_S<-log10(ES_S)
X_E <-E_mo-ES_E
X_S <-S_mo-ES_S
#-----------------------------------function------------------------------------
plot.abundance<-function(mo){
  mo_t<-rbind(t,mo)
  colnames(mo_t)<-t
  mo_t<-as.data.frame(t(mo_t))
  mo_t<-melt(mo_t, id="V1")
  mo.abundance<-ggplot((mo_t),aes())
  mo.abundance<-mo.abundance+geom_line(data = mo_t,aes(x=V1,y=value,group=variable),color="#7fc3bd")+
    labs(x=NULL)+
    labs(y=NULL)+
    theme_bw()+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
  theme(panel.grid=element_blank())
  return(mo.abundance)
}
options (warn = -1)
#-----------------------------------Fig1.A--------------------------------------
#-----------------------------------1st.pic.E.mo.abundance----------------------
#8086a8
plot.abundance<-function(mo){
  mo_t<-rbind(t,mo)
  colnames(mo_t)<-t
  mo_t<-as.data.frame(t(mo_t))
  mo_t<-melt(mo_t, id="V1")
  mo.abundance<-ggplot((mo_t),aes())
  mo.abundance<-mo.abundance+geom_line(data=mo_t,aes(x=V1,y=value,group=variable),color="#527dad",alpha=0.4)+
    labs(x=NULL)+
    labs(y=NULL)+scale_y_continuous(labels=c(0,expression(10^2.5),expression(10^5),expression(10^7.5),expression(10^10)),
                                    limits = c(0,10))+
    theme_bw()+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))+
    theme(panel.grid=element_blank())
  return(mo.abundance)
}
E.mo.abundance<-plot.abundance(E_mo)
E.mo.abundance
#-----------------------------------1st.pic.E.mo.abundance's Curve fitting------
H0 <- function(mo){
  source("function.R")
  par=c(18,2,0.2,1,0.5)
  H0 = function(yt,t,par){
    miu = get_miu(par=c(par[1],par[2],par[3]),t)
    sigma = SAD1_get_matrix(par = c(par[4],par[5]), times = t )
    L0 = c()
    L0 = sum(dmvnorm(yt,miu,sigma,log = T))
    return(-L0)
  }
  H0_par<-optim(par=par,H0,yt=mo,t=t,method="BFGS",control=list(maxit=5000))
  H0_par<-optim(par=H0_par$par,H0,yt=mo,t=t,method="BFGS",control=list(maxit=5000))
  H0_frame <- as.data.frame(cbind(get_miu(H0_par$par[1:3],t=c(0.5:20,20)),t=c(0.5:20,20)))
  colnames(H0_frame)<-c("miu","time")
  b <- list()
  b[[1]] <- H0_frame
  b[[2]] <- H0_par$par[1:5]
  return(b)
}
abcuvre<- function(mo){
  mo_H0<- H0(mo)
  cat(mo_H0[[2]],log(mo_H0[[2]][2])/mo_H0[[2]][3],get_miu(par=c(mo_H0[[2]][1:3]),log(mo_H0[[2]][2])/mo_H0[[2]][3]))
  a <- c(log(mo_H0[[2]][2])/mo_H0[[2]][3],get_miu(par=c(mo_H0[[2]][1:3]),log(mo_H0[[2]][2])/mo_H0[[2]][3]))
  E.mo.abundance<-E.mo.abundance+geom_line(data=mo_H0[[1]],aes(x=time,y=miu),size=2,color="#41648B")
  #geom_segment(aes(x =a[1], y = 0, xend = a[1], yend = a[2]),colour="#41648B")
  return(E.mo.abundance)
}
E.mo_abline<- abcuvre(E_mo)
E.mo_abline
#-----------------------------------2nd.pic.ES_e.abundance----------------------
#527dad 
plot.abundance<-function(mo){
  mo_t<-rbind(t,mo)
  colnames(mo_t)<-t
  mo_t<-as.data.frame(t(mo_t))
  mo_t<-melt(mo_t, id="V1")
  mo.abundance<-ggplot((mo_t),aes())
  mo.abundance<-mo.abundance+geom_line(data = mo_t,aes(x=V1,y=value,group=variable),color="#527dad",alpha=0.4)+
    labs(x=NULL)+
    labs(y=NULL)+scale_y_continuous(labels=c(0,expression(10^2.5),expression(10^5),expression(10^7.5),expression(10^10)),
                                    limits = c(0,10))+
    theme_bw()+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))+
    theme(panel.grid=element_blank())
  return(mo.abundance)
}
ES_E.abundance<-plot.abundance(ES_E)
ES_E.abundance
#-----------------------------------2nd.pic.ES_E.abundance's Curve fitting------
H0 <- function(mo){
  source("R:/Rcode/肠道微生物/画三张大图需要的函数.R")
  par=c(8,1,3,0.08,5)
  H0 = function(yt,t,par){
    miu = get_miu(par=c(par[1],par[2],par[3]),t)
    sigma = SAD1_get_matrix(par = c(par[4],par[5]), times = t )
    L0 = c()
    L0 = sum(dmvnorm(yt,miu,sigma,log = T))
    return(-L0)
  }
  H0_par<-optim(par=par,H0,yt=mo,t=t,method="BFGS",control=list(maxit=5000))
  H0_par<-optim(par=H0_par$par,H0,yt=mo,t=t,method="BFGS",control=list(maxit=5000))
  H0_frame <- as.data.frame(cbind(get_miu(H0_par$par[1:3],t=c(0.5:20,20)),t=c(0.5:20,20)))
  colnames(H0_frame)<-c("miu","time")
  b <- list()
  b[[1]] <- H0_frame
  b[[2]] <- H0_par$par[1:5]
  return(b)
}
abcuvre<- function(mo){
  mo_H0<- H0(mo)
  cat(mo_H0[[2]],log(mo_H0[[2]][2])/mo_H0[[2]][3],get_miu(par=c(mo_H0[[2]][1:3]),log(mo_H0[[2]][2])/mo_H0[[2]][3]))
  a <- c(log(mo_H0[[2]][2])/mo_H0[[2]][3],get_miu(par=c(mo_H0[[2]][1:3]),log(mo_H0[[2]][2])/mo_H0[[2]][3]))
  ES_E.abundance<-ES_E.abundance+geom_line(data=mo_H0[[1]],aes(x=time,y=miu),size=2,color="#41648B")
  #geom_segment(aes(x =a[1],y = 0,xend = a[1],yend = a[2]),colour="#41648B")
  return(ES_E.abundance)
}
ES_E_abline<- abcuvre(ES_E)
ES_E_abline
#-----------------------------------3rd.pic.ES_s.abundance----------------------
#7fc3bd  "#5CB2AB"
plot.abundance<-function(mo){
  mo_t<-rbind(t,mo)
  colnames(mo_t)<-t
  mo_t<-as.data.frame(t(mo_t))
  mo_t<-melt(mo_t, id="V1")
  mo.abundance<-ggplot((mo_t),aes())
  mo.abundance<-mo.abundance+geom_line(data = mo_t,aes(x=V1,y=value,group=variable),color="#a67e6c",alpha=0.4)+
    labs(x=NULL)+
    labs(y=NULL)+scale_y_continuous(labels=c(0,expression(10^2.5),expression(10^5),expression(10^7.5),expression(10^10)),
                                    limits = c(0,10))+
    theme_bw()+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))+
    theme(panel.grid=element_blank())
  return(mo.abundance)
}
ES_S.abundance<-plot.abundance(ES_S)
ES_S.abundance
#-----------------------------------3rd.pic.ES_s.abundance's Curve fitting------
H0 <- function(mo){
  source("R:/Rcode/肠道微生物/画三张大图需要的函数.R")
  par=c(10,1,3,0.08,5)
  H0 = function(yt,t,par){
    miu = get_miu(par=c(par[1],par[2],par[3]),t)
    sigma = SAD1_get_matrix(par = c(par[4],par[5]), times = t )
    L0 = c()
    L0 = sum(dmvnorm(yt,miu,sigma,log = T))
    return(-L0)
  }
  H0_par<-optim(par=par,H0,yt=mo,t=t,method="BFGS",control=list(maxit=5000))
  H0_par<-optim(par=H0_par$par,H0,yt=mo,t=t,method="BFGS",control=list(maxit=5000))
  H0_frame <- as.data.frame(cbind(get_miu(H0_par$par[1:3],t=c(0.5:20,20)),t=c(0.5:20,20)))
  colnames(H0_frame)<-c("miu","time")
  b <- list()
  b[[1]] <- H0_frame
  b[[2]] <- H0_par$par[1:5]
  return(b)
}
abcuvre<- function(mo){
  mo_H0<- H0(mo)
  cat(mo_H0[[2]],log(mo_H0[[2]][2])/mo_H0[[2]][3],get_miu(par=c(mo_H0[[2]][1:3]),log(mo_H0[[2]][2])/mo_H0[[2]][3]))
  a <- c(log(mo_H0[[2]][2])/mo_H0[[2]][3],get_miu(par=c(mo_H0[[2]][1:3]),log(mo_H0[[2]][2])/mo_H0[[2]][3]))
  ES_S.abundance<-ES_S.abundance+geom_line(data=mo_H0[[1]],aes(x=time,y=miu),size=2,color="#8C6554")
  return(ES_S.abundance)
}
ES_S_abline<- abcuvre(ES_S)
ES_S_abline
#-----------------------------------4th.pic.S.mo.abundance----------------------
#a67e6c "#8C6554"
plot.abundance<-function(mo){
  mo_t<-rbind(t,mo)
  colnames(mo_t)<-t
  mo_t<-as.data.frame(t(mo_t))
  mo_t<-melt(mo_t, id="V1")
  mo.abundance<-ggplot((mo_t),aes())
  mo.abundance<-mo.abundance+geom_line(data = mo_t,aes(x=V1,y=value,group=variable),color="#a67e6c",alpha=0.4)+
    labs(x=NULL)+
    labs(y=NULL)+scale_y_continuous(labels=c(0,expression(10^2.5),expression(10^5),expression(10^7.5),expression(10^10)),
                                    limits = c(0,10))+
    theme_bw()+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))+
    theme(panel.grid=element_blank())
  return(mo.abundance)
}
S.mo.abundance<-plot.abundance(S_mo)
S.mo.abundance
#-----------------------------------4th.pic.S.mo.abundance----------------------
H0 <- function(mo){
  source("R:/Rcode/肠道微生物/画三张大图需要的函数.R")
  par=c(8,2.5,2,0.08,5)
  H0 = function(yt,t,par){
    miu = get_miu(par=c(par[1],par[2],par[3]),t)
    sigma = SAD1_get_matrix(par = c(par[4],par[5]), times = t )
    L0 = c()
    L0 = sum(dmvnorm(yt,miu,sigma,log = T))
    return(-L0)
  }
  H0_par<-optim(par=par,H0,yt=mo,t=t,method="BFGS",control=list(maxit=5000))
  H0_par<-optim(par=H0_par$par,H0,yt=mo,t=t,method="BFGS",control=list(maxit=5000))
  H0_frame <- as.data.frame(cbind(get_miu(H0_par$par[1:3],t=c(0.5:20,20)),t=c(0.5:20,20)))
  colnames(H0_frame)<-c("miu","time")
  b <- list()
  b[[1]] <- H0_frame
  b[[2]] <- H0_par$par[1:5]
  return(b)
}
abcuvre<- function(mo){
  mo_H0<- H0(mo)
  cat(mo_H0[[2]],log(mo_H0[[2]][2])/mo_H0[[2]][3],get_miu(par=c(mo_H0[[2]][1:3]),log(mo_H0[[2]][2])/mo_H0[[2]][3]))
  a <- c(log(mo_H0[[2]][2])/mo_H0[[2]][3],get_miu(par=c(mo_H0[[2]][1:3]),log(mo_H0[[2]][2])/mo_H0[[2]][3]))
  S.mo.abundance<-S.mo.abundance+geom_line(data=mo_H0[[1]],aes(x=time,y=miu),size=2,color="#8C6554")
  return(S.mo.abundance)
}
S.mo_abline<- abcuvre(S_mo)
S.mo_abline
abundance_abline<-cowplot::plot_grid(E.mo_abline,ES_E_abline,ES_S_abline,S.mo_abline,nrow = 1)
#-----------------------------------Fig1.B--------------------------------------
#-----------------------------------function------------------------------------
plot.plasticity<-function(mo){
  mo_t<-rbind(t,mo)
  colnames(mo_t)<-t
  mo_t<-as.data.frame(t(mo_t))
  mo_t<-melt(mo_t, id="V1")
  mo.abundance<-ggplot((mo_t),aes())
  mo.abundance<-mo.abundance+geom_line(data = mo_t,aes(x=V1,y=value,group=variable),color="#5f867e")+labs(x=NULL)+labs(y=NULL)+theme_bw()+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))+theme(panel.grid=element_blank())
  return(mo.abundance)
}
#-----------------------------------1st.X.E.plasticity--------------------------
#a47d96
plot.plasticity<-function(mo){
  mo_t<-rbind(t,mo)
  colnames(mo_t)<-t
  mo_t<-as.data.frame(t(mo_t))
  mo_t<-melt(mo_t, id="V1")
  mo.abundance<-ggplot((mo_t),aes())
  mo.abundance<-mo.abundance+geom_line(data = mo_t,aes(x=V1,y=value,group=variable),color="#6388B6")+
    labs(x=NULL)+labs(y=NULL)+geom_hline(aes(yintercept=0),size=0.8,linetype="dashed",color="#585D7E")+
    scale_y_continuous(labels=c(expression(-10^5),expression(-10^2.5),0,expression(10^2.5),expression(10^5)),limits = c(-6.5,5.5))+
    theme_bw()+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))+theme(panel.grid=element_blank())
  return(mo.abundance)
}
X.E.plasticity<-plot.plasticity(X_E)
X.E.plasticity
#-----------------------------------2nd.ES.E.plasticity-------------------------
#6ba8a2
zhuanhuan<-function(mo)
{mo_t<-rbind(t,mo)
colnames(mo_t)<-t
mo_t<-as.data.frame(t(mo_t))
mo_t<-melt(mo_t, id="V1")}
co.es.e.mat<-(cbind(zhuanhuan(ES_E),zhuanhuan(E_mo))[,-(4:5)])
colnames(co.es.e.mat)<-c("V1","var","co.e","e")
X.ES.E.plasticity<-ggplot(NULL)
X.ES.E.plasticity<-X.ES.E.plasticity+
  geom_point(data=co.es.e.mat,aes(x=e,y=co.e,group=var),color="#6388B6")+
  scale_y_continuous(labels=(math_format(10^.x)),limits = c(0.5,10))+
  scale_x_continuous(labels=(math_format(10^.x)),limits = c(0.5,10))+
  theme_bw()+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))+labs(x=NULL)+labs(y="Co-culture")+
  theme(axis.title.y=element_text(size=20))+
  theme(panel.grid=element_blank())+
  geom_abline(slope=1,intercept=(0))
X.ES.E.plasticity
#-----------------------------------3rd.ES.S.plasticity-------------------------
#a4a57d
co.es.s.mat<-(cbind(zhuanhuan(ES_S),zhuanhuan(S_mo))[,-(4:5)])
colnames(co.es.s.mat)<-c("V1","var","co.s","s")
X.ES.S.plasticity<-ggplot(NULL)
X.ES.S.plasticity<-X.ES.S.plasticity+geom_point(data=co.es.s.mat,aes(x=s,y=co.s,group=var),color="#BA9C8C")+
  scale_y_continuous(labels=(math_format(10^.x)),limits = c(0.5,10))+
  scale_x_continuous(labels=(math_format(10^.x)),limits = c(0.5,10))+
  theme_bw()+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))+labs(x=NULL)+labs(y=NULL)+
  theme(panel.grid=element_blank())+
  geom_abline(slope=1,intercept=(0))
X.ES.S.plasticity
#-----------------------------------4th.X.S.plasticity--------------------------
#5f867e
plot.plasticity<-function(mo){
  mo_t<-rbind(t,mo)
  colnames(mo_t)<-t
  mo_t<-as.data.frame(t(mo_t))
  mo_t<-melt(mo_t, id="V1")
  mo.abundance<-ggplot((mo_t),aes())
  mo.abundance<-mo.abundance+geom_line(data = mo_t,aes(x=V1,y=value,group=variable),color="#BA9C8C")+
    labs(x=NULL)+labs(y=NULL)+theme_bw()+theme(axis.text.x=element_text(size=14))+
    scale_y_continuous(labels=c(expression(-10^5),expression(-10^2.5),0,expression(10^2.5),expression(10^5)),limits = c(-6.5,5.5))+
    theme(axis.text.y=element_text(size=14))+theme(panel.grid=element_blank())+
    labs(y="Phenotypic plasticity")+theme(axis.title.y = (element_text(size = 20)))+
    geom_hline(aes(yintercept=0),size=0.8,linetype="dashed",color="#004E75")
  return(mo.abundance)
}
X.S.plasticity<-plot.plasticity(X_S)
X.S.plasticity
#-----------------------------------output--------------------------------------
plasticity<-cowplot::plot_grid(X.E.plasticity,X.ES.E.plasticity,X.ES.S.plasticity,X.S.plasticity,nrow = 1)
Fig1abline<-cowplot::plot_grid(abundance_abline,plasticity,nrow = 2)
Fig1abline
ggsave("Fig1abline",Fig1abline,width = 16,height = 8)
#-----------------------------------Resampling----------------------------------
#-----------------------------------E.mo----------------------------------------
E.mo.var <- matrix(NA,100,3)
for(i in 1:100){
  set.seed(i)
  roww <- sample(c(1:100),100,replace = T)
  E_mos <- E_mo[roww,]
  H0 <- function(mo){
    source("function.R")
    par=c(8,1,0.2,1,0.5)
    H0 = function(yt,t,par){
      miu = get_miu(par=c(par[1],par[2],par[3]),t)
      sigma = SAD1_get_matrix(par = c(par[4],par[5]), times = t )
      L0 = c()
      L0 = sum(dmvnorm(yt,miu,sigma,log = T))
      return(-L0)
    }
    H0_par<-optim(par=par,H0,yt=mo,t=t,method="BFGS",control=list(maxit=5000))
    H0_par<-optim(par=H0_par$par,H0,yt=mo,t=t,method="BFGS",control=list(maxit=5000))
    H0_frame <- as.data.frame(cbind(get_miu(H0_par$par[1:3],t=c(0.5:20,20)),t=c(0.5:20,20)))
    colnames(H0_frame)<-c("miu","time")
    b <- list()
    b[[1]] <- H0_frame
    b[[2]] <- H0_par$par[1:3]
    return(b)
  }
  mo_H0<- H0(E_mos)
  E.mo.var[i,]<- mo_H0[[2]]
}
E.mo.var
sd(E.mo.var[,1])
sd(E.mo.var[,2])
sd(E.mo.var[,3])
#-----------------------------------ES_E----------------------------------------
ES_E.var <- matrix(NA,100,3)
for(i in 1:100){
  set.seed(i)
  roww <- sample(c(1:100),100,replace = T)
  ES_Es <- ES_E[roww,]
  H0 <- function(mo){
    source("R:/Rcode/肠道微生物/画三张大图需要的函数.R")
    par=c(8,1,0.5,1,0.3)
    H0 = function(yt,t,par){
      miu = get_miu(par=c(par[1],par[2],par[3]),t)
      sigma = SAD1_get_matrix(par = c(par[4],par[5]), times = t )
      L0 = c()
      L0 = sum(dmvnorm(yt,miu,sigma,log = T))
      return(-L0)
    }
    H0_par<-optim(par=par,H0,yt=mo,t=t,method="BFGS",control=list(maxit=5000))
    H0_par<-optim(par=H0_par$par,H0,yt=mo,t=t,method="BFGS",control=list(maxit=5000))
    H0_frame <- as.data.frame(cbind(get_miu(H0_par$par[1:3],t=c(0.5:20,20)),t=c(0.5:20,20)))
    colnames(H0_frame)<-c("miu","time")
    b <- list()
    b[[1]] <- H0_frame
    b[[2]] <- H0_par$par[1:3]
    return(b)
  }
  mo_H0<- H0(ES_Es)
  ES_E.var[i,]<- mo_H0[[2]]
}
ES_E.var
sd(ES_E.var[,1])
sd(ES_E.var[,2])
sd(ES_E.var[,3])
#-----------------------------------ES_S的标准差--------------------------------
for(i in 1:100){
  set.seed(i)
  roww <- sample(c(1:100),100,replace = T)
  ES_Ss <- ES_S[roww,]
  H0 <- function(mo){
    source("R:/Rcode/肠道微生物/画三张大图需要的函数.R")
    par=c(7,1.5,0.2,0.9,0.25)
    H0 = function(yt,t,par){
      miu = get_miu(par=c(par[1],par[2],par[3]),t)
      sigma = SAD1_get_matrix(par = c(par[4],par[5]), times = t )
      L0 = c()
      L0 = sum(dmvnorm(yt,miu,sigma,log = T))
      return(-L0)
    }
    H0_par<-optim(par=par,H0,yt=mo,t=t,method="BFGS",control=list(maxit=5000))
    H0_par<-optim(par=H0_par$par,H0,yt=mo,t=t,method="BFGS",control=list(maxit=5000))
    H0_frame <- as.data.frame(cbind(get_miu(H0_par$par[1:3],t=c(0.5:20,20)),t=c(0.5:20,20)))
    colnames(H0_frame)<-c("miu","time")
    b <- list()
    b[[1]] <- H0_frame
    b[[2]] <- H0_par$par[1:3]
    return(b)
  }
  mo_H0<- H0(ES_Ss )
  ES_S.var[i,]<- mo_H0[[2]]
}
ES_S.var
sd(ES_S.var[,1])
sd(ES_S.var[,2])
sd(ES_S.var[,3])
#-----------------------------------s.mo的标准差--------------------------------
s.mo.var <- matrix(NA,100,3)
for(i in 1:100){
  set.seed(i)
  roww <- sample(c(1:100),100,replace = T)
  S_mos <- S_mo[roww,]
  H0 <- function(mo){
    source("R:/Rcode/肠道微生物/画三张大图需要的函数.R")
    par=c(8.45,3.15,0.28,0.8,-0.4)
    H0 = function(yt,t,par){
      miu = get_miu(par=c(par[1],par[2],par[3]),t)
      sigma = SAD1_get_matrix(par = c(par[4],par[5]), times = t )
      L0 = c()
      L0 = sum(dmvnorm(yt,miu,sigma,log = T))
      return(-L0)
    }
    H0_par<-optim(par=par,H0,yt=mo,t=t,method="BFGS",control=list(maxit=5000))
    H0_par<-optim(par=H0_par$par,H0,yt=mo,t=t,method="BFGS",control=list(maxit=5000))
    H0_frame <- as.data.frame(cbind(get_miu(H0_par$par[1:3],t=c(0.5:20,20)),t=c(0.5:20,20)))
    colnames(H0_frame)<-c("miu","time")
    b <- list()
    b[[1]] <- H0_frame
    b[[2]] <- H0_par$par[1:3]
    return(b)
  }
  mo_H0<- H0(S_mos)
  s.mo.var[i,]<- mo_H0[[2]]
}
s.mo.var
sd(s.mo.var[,1])
sd(s.mo.var[,2])
sd(s.mo.var[,3])