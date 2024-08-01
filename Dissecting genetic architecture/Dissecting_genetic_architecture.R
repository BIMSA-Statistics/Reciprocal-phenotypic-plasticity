#---------------------requiredPackages----------------------------
requiredPackages = c("mvtnorm","parallel","scales","pbapply","plyr","scales","patchwork","ggplot2","deSolve","Rmisc","glmnet","gridExtra","tidyverse","ggpubr")
for(packages in requiredPackages){
  if(!require(packages,character.only = TRUE)) install.packages(packages)
  require(packages,character.only = TRUE)
}
#---------------------inputdata-----------------------------------
all.s.sig<-read.csv("all.s.sig.csv",row.names=1)
all.s.sig.loci<-read.csv("loci.s.0.01.csv",row.names=1)
S_SNP <- read.table("S-SNP.txt",row.names = 1,header=T)
E_SNP <- read.table("E-SNP.txt",row.names = 1,header=T)
all.e.sig<-read.csv("all.e.sig.csv",row.names=1)
all.e.sig.loci<-read.csv("loci.ee.0.01.csv",row.names=1)
#---------------------Data preprocessing--------------------------
t=c(0.5,1,1.5,2,4,6,8,10,12,16,20)
get_miu<-function(par, t, options=list())
{
  y <- par[1]/(1+par[2]*exp(-par[3]*t))-par[4]/(1+par[5]*exp(-par[6]*t))
  return (y);
}
es.co <- rbind(all.e.sig.loci[,c(1,2)],all.s.sig.loci[,c(1,2)])
se.co <- rbind(all.s.sig.loci[,c(1,2)],all.e.sig.loci[,c(1,2)])
es.130 <- all.s.sig.loci[which(duplicated(es.co))-132535,]#S表型下的
se.130 <- all.e.sig.loci[which(duplicated(se.co))-3543,]#E表型下的
ba <- c()
for(i in 1:130){
  a <- which(all.s.sig.loci[which(all.s.sig.loci[,2]==es.130[,2][i]),][,1]==es.130[,1][i])
  b1 <- which(all.s.sig[,1]==rownames(S_SNP[es.130[,2][i],]))[a]
  ba <- c(ba,b1)
}

be <- c()
for(i in 1:130){
  a <- which(all.e.sig.loci[which(all.e.sig.loci[,1]==es.130[,1][i]),][,2]==es.130[,2][i])#13
  b1 <- which(all.e.sig[,1]==rownames(E_SNP[es.130[,1][i],]))[a]
  be <- c(be,b1)
}
#---------------------ALL.S.Significant.location.average_value--------------------------
s.all.effect<-list()
s.all.aee.list <- list()
s.all.aes.list <- list()
s.all.ies.list <- list()
for (z in 1:dim(all.s.sig)[1]){
  s.all.miu11<-c(as.numeric(all.s.sig[z,3:8]))
  s.all.miu10<-c(as.numeric(all.s.sig[z,9:14]))
  s.all.miu00<-c(as.numeric(all.s.sig[z,15:20]))
  s.all.miu01<-c(as.numeric(all.s.sig[z,21:26]))
  
  s.all.mmiu11 <- get_miu(s.all.miu11,t)
  s.all.mmiu10 <- get_miu(s.all.miu10,t)
  s.all.mmiu00 <- get_miu(s.all.miu00,t)
  s.all.mmiu01 <- get_miu(s.all.miu01,t)
  
  s.all.aee <- c()
  s.all.aes <- c()
  s.all.ies <- c()
  s.all.aee<-(s.all.mmiu11+s.all.mmiu10-s.all.mmiu01-s.all.mmiu00)/4
  s.all.aes<-(s.all.mmiu11+s.all.mmiu01-s.all.mmiu10-s.all.mmiu00)/4
  s.all.ies<-(s.all.mmiu11+s.all.mmiu00-s.all.mmiu10-s.all.mmiu01)/4
  
  s.all.aee.list[[z]]<-s.all.aee
  s.all.aes.list[[z]]<-s.all.aes
  s.all.ies.list[[z]]<-s.all.ies
}
s.all.aee.ave<-Reduce('+',s.all.aee.list)/3543
s.all.aes.ave<-Reduce('+',s.all.aes.list)/3543
s.all.ies.ave<-Reduce('+',s.all.ies.list)/3543
mat<-data.frame(t,s.all.aee.ave,s.all.aes.ave,s.all.ies.ave)
colnames(mat)<-c("t","aee","aes","ies")
mat[,c(2:4)] <- log10(exp(mat[,c(2:4)]))
s.all.gg<- ggplot()
s.all.gg<-s.all.gg+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
s.all.gg<-s.all.gg+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
s.all.gg<-s.all.gg+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
s.all.gg<-s.all.gg+ theme_bw()+scale_y_continuous(labels=c(expression(-10^0.25),0,expression(10^0.25),
                                                           expression(10^0.5),expression(10^0.75)),limits = c(-0.3,0.75),position = "right")+
  theme(panel.grid=element_blank())+
  theme(plot.margin = unit(c(1,0,0,0),"lines"))+scale_x_continuous(breaks = c(5,10,15,20))+
  geom_hline(aes(yintercept=0),linetype="dotted")+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
s.all.effect.ave<-s.all.gg
s.all.effect.ave

#---------------------130.S.Significant.location.average_value------------------
s.130.aee.ave<-Reduce('+',s.all.aee.list[ba])/130
s.130.aes.ave<-Reduce('+',s.all.aes.list[ba])/130
s.130.ies.ave<-Reduce('+',s.all.ies.list[ba])/130
mat<-data.frame(t,s.130.aee.ave,s.130.aes.ave,s.130.ies.ave)
colnames(mat)<-c("t","aee","aes","ies")
mat[,c(2:4)] <- log10(exp(mat[,c(2:4)]))
s.130.gg<- ggplot()
s.130.gg<-s.130.gg+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
s.130.gg<-s.130.gg+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
s.130.gg<-s.130.gg+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
s.130.gg<-s.130.gg+ theme_bw()+scale_y_continuous(labels=c(expression(-10^1),expression(-10^0.5),0,
                                                           expression(10^0.5),expression(10^0.75)),
                                                  limits = c(-1,0.75),position = "right")+
  theme(panel.grid=element_blank(),panel.background = element_rect(color = 'black', fill = "#FFEBEE"))+
  theme(plot.margin = unit(c(1,0,0,0),"lines"))+scale_x_continuous(breaks = c(5,10,15,20))+
  geom_hline(aes(yintercept=0),linetype="dotted")+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
s.130.effect.ave<-s.130.gg
s.130.effect.ave
#---------------------RM130.S.Significant.location.average_value----------------
s.rm130.aee.ave<-Reduce('+',s.all.aee.list[-ba])/(dim(all.s.sig)[1]-130)
s.rm130.aes.ave<-Reduce('+',s.all.aes.list[-ba])/(dim(all.s.sig)[1]-130)
s.rm130.ies.ave<-Reduce('+',s.all.ies.list[-ba])/(dim(all.s.sig)[1]-130)
mat<-data.frame(t,s.rm130.aee.ave,s.rm130.aes.ave,s.rm130.ies.ave)
colnames(mat)<-c("t","aee","aes","ies")
mat[,c(2:4)] <- log10(exp(mat[,c(2:4)]))
s.rm130.gg<- ggplot()
s.rm130.gg<-s.rm130.gg+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
s.rm130.gg<-s.rm130.gg+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
s.rm130.gg<-s.rm130.gg+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
s.rm130.gg<-s.rm130.gg+ theme_bw()+scale_y_continuous(labels=c(expression(-10^0.25),0,
                                                               expression(10^0.25),expression(10^0.5),expression(10^0.75)),
                                                      limits = c(-0.3,0.75),position = "right")+
  theme(panel.grid=element_blank(),panel.background = element_rect(color = 'black', fill = "#F0FBFF"))+
  theme(plot.margin = unit(c(1,0,0,0),"lines"))+scale_x_continuous(breaks = c(5,10,15,20))+
  geom_hline(aes(yintercept=0),linetype="dotted")+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
s.rm130.effect.ave<-s.rm130.gg
s.rm130.effect.ave
#---------------------ALL.S.Significant.location.variance--------------------------
s.all.aee.var.list <- list()
s.all.aes.var.list <- list()
s.all.ies.var.list <- list()
for(i in 1:dim(all.s.sig)[1]){
  p1 <- length(which(S_SNP[all.s.sig.loci[,2],][i,]==1))/100
  q1 <- length(which(S_SNP[all.s.sig.loci[,2],][i,]==0))/100
  p2 <- length(which(E_SNP[all.s.sig.loci[,1],][i,]==1))/100
  q2 <- length(which(E_SNP[all.s.sig.loci[,1],][i,]==0))/100
  
  s.all.aee.var.list[[i]]<-(1-(p1-q1)^2)*(s.all.aee.list[[i]]+(p2-q2)*s.all.ies.list[[i]])^2
  s.all.aes.var.list[[i]]<-(1-(p2-q2)^2)*(s.all.aes.list[[i]]+(p1-q1)*s.all.ies.list[[i]])^2
  s.all.ies.var.list[[i]]<-(1-(p1-q1)^2-(p2-q2)^2+(p1-q1)^2*(p2-q2)^2)*(s.all.ies.list[[i]])^2
  
}
s.all.aee.ave.var <-Reduce('+',s.all.aee.var.list)/dim(all.s.sig)[1]
s.all.aes.ave.var <-Reduce('+',s.all.aes.var.list)/dim(all.s.sig)[1]
s.all.ies.ave.var <-Reduce('+',s.all.ies.var.list)/dim(all.s.sig)[1]
mat<-data.frame(t,s.all.aee.ave.var,s.all.aes.ave.var,s.all.ies.ave.var)
mat[,c(2:4)] <- log10(exp(mat[,c(2:4)]))*dim(all.s.sig)[1]
colnames(mat)<-c("t","aee","aes","ies")
s.all.gg.v<- ggplot()
s.all.gg.v<-s.all.gg.v+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
s.all.gg.v<-s.all.gg.v+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
s.all.gg.v<-s.all.gg.v+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
s.all.gg.v<-s.all.gg.v+ theme_bw()+
  scale_y_continuous(labels=(math_format(10^.x)),position = "right")+
  theme(panel.grid=element_blank())+scale_x_continuous(breaks = c(5,10,15,20))+
  theme(plot.margin = unit(c(1,0,0,0),"lines"))+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
s.all.effect.ave.var<-s.all.gg.v
s.all.effect.ave.var
#---------------------130.S.Significant.location.variance--------------------------
s.130.aee.ave.var <-Reduce('+',s.all.aee.var.list[ba])/130
s.130.aes.ave.var <-Reduce('+',s.all.aes.var.list[ba])/130
s.130.ies.ave.var <-Reduce('+',s.all.ies.var.list[ba])/130
mat<-data.frame(t,s.130.aee.ave.var,s.130.aes.ave.var,s.130.ies.ave.var)
mat[,c(2:4)] <- log10(exp(mat[,c(2:4)]))*130
colnames(mat)<-c("t","aee","aes","ies")
s.130.gg.v<- ggplot()
s.130.gg.v<-s.130.gg.v+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
s.130.gg.v<-s.130.gg.v+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
s.130.gg.v<-s.130.gg.v+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
s.130.gg.v<-s.130.gg.v+ theme_bw()+
  scale_y_continuous(labels=(math_format(10^.x)),position = "right")+
  theme(panel.grid=element_blank(),panel.background = element_rect(color = 'black', fill = "#FFEBEE"))+scale_x_continuous(breaks = c(5,10,15,20))+
  theme(plot.margin = unit(c(1,0,0,0),"lines"))+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
s.130.effect.ave.var<-s.130.gg.v
s.130.effect.ave.var

#---------------------RM130.S.Significant.location.variance------------------------
s.rm130.aee.ave.var <-Reduce('+',s.all.aee.var.list[-ba])/(dim(all.s.sig)[1]-130)
s.rm130.aes.ave.var <-Reduce('+',s.all.aes.var.list[-ba])/(dim(all.s.sig)[1]-130)
s.rm130.ies.ave.var <-Reduce('+',s.all.ies.var.list[-ba])/(dim(all.s.sig)[1]-130)
mat<-data.frame(t,s.rm130.aee.ave.var,s.rm130.aes.ave.var,s.rm130.ies.ave.var)
mat[,c(2:4)] <- log10(exp(mat[,c(2:4)]))*(dim(all.s.sig)[1]-130)
colnames(mat)<-c("t","aee","aes","ies")
s.rm130.gg.v<- ggplot()
s.rm130.gg.v<-s.rm130.gg.v+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
s.rm130.gg.v<-s.rm130.gg.v+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
s.rm130.gg.v<-s.rm130.gg.v+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
s.rm130.gg.v<-s.rm130.gg.v+ theme_bw()+
  scale_y_continuous(labels=(math_format(10^.x)),position = "right")+
  theme(panel.grid=element_blank(),panel.background = element_rect(color = 'black', fill = "#F0FBFF"))+scale_x_continuous(breaks = c(5,10,15,20))+
  theme(plot.margin = unit(c(1,0,0,0),"lines"))+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
s.rm130.effect.ave.var<-s.rm130.gg.v
s.rm130.effect.ave.var


#---------------------ALL.S.Significant.location.Proportion--------------------------
s.all.aee.pro.list <- list()
s.all.aes.pro.list <- list()
s.all.ies.pro.list <- list()
s.all.pro.list <- list()

for(i in 1:dim(all.s.sig)[1]){
  s.all.pro.list[[i]] <- s.all.aee.var.list[[i]]+s.all.aes.var.list[[i]]+s.all.ies.var.list[[i]]
  s.all.aee.pro.list[[i]]<-s.all.aee.var.list[[i]]/s.all.pro.list[[i]]
  s.all.aes.pro.list[[i]]<-s.all.aes.var.list[[i]]/s.all.pro.list[[i]]
  s.all.ies.pro.list[[i]]<-s.all.ies.var.list[[i]]/s.all.pro.list[[i]]
}

s.all.aee.ave.pro <-Reduce('+',s.all.aee.pro.list)/dim(all.s.sig)[1]
s.all.aes.ave.pro <-Reduce('+',s.all.aes.pro.list)/dim(all.s.sig)[1]
s.all.ies.ave.pro <-Reduce('+',s.all.ies.pro.list)/dim(all.s.sig)[1]

mat<-data.frame(t,s.all.aee.ave.pro,
                s.all.aes.ave.pro,
                s.all.ies.ave.pro)
colnames(mat)<-c("t","aee","aes","ies")
s.all.gg.p<- ggplot()
s.all.gg.p<-s.all.gg.p+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
s.all.gg.p<-s.all.gg.p+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
s.all.gg.p<-s.all.gg.p+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
s.all.gg.p<-s.all.gg.p+ theme_bw()+scale_y_continuous(position = "right",limits = c(0,1))+
  theme(panel.grid=element_blank())+scale_x_continuous(breaks = c(5,10,15,20))+
  theme(plot.margin = unit(c(1,0,0,0),"lines"))+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
s.all.effect.ave.pro<-s.all.gg.p
s.all.effect.ave.pro
#---------------------130.S.Significant.location.Proportion --------
s.130.aee.ave.pro <-Reduce('+',s.all.aee.pro.list[ba])/130
s.130.aes.ave.pro <-Reduce('+',s.all.aes.pro.list[ba])/130
s.130.ies.ave.pro <-Reduce('+',s.all.ies.pro.list[ba])/130

mat<-data.frame(t,s.130.aee.ave.pro,
                s.130.aes.ave.pro,
                s.130.ies.ave.pro)
colnames(mat)<-c("t","aee","aes","ies")
s.130.gg.p<- ggplot()
s.130.gg.p<-s.130.gg.p+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
s.130.gg.p<-s.130.gg.p+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
s.130.gg.p<-s.130.gg.p+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
s.130.gg.p<-s.130.gg.p+ theme_bw()+scale_y_continuous(position = "right",limits = c(0,1))+
  theme(panel.grid=element_blank(),panel.background = element_rect(color = 'black', fill = "#FFEBEE"))+scale_x_continuous(breaks = c(5,10,15,20))+
  theme(plot.margin = unit(c(1,0,0,0),"lines"))+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
s.130.effect.ave.pro<-s.130.gg.p
s.130.effect.ave.pro
#---------------------RM130.S.Significant.location.Proportion------------------------
s.rm130.aee.ave.pro <-Reduce('+',s.all.aee.pro.list[-ba])/(dim(all.s.sig)[1]-130)
s.rm130.aes.ave.pro <-Reduce('+',s.all.aes.pro.list[-ba])/(dim(all.s.sig)[1]-130)
s.rm130.ies.ave.pro <-Reduce('+',s.all.ies.pro.list[-ba])/(dim(all.s.sig)[1]-130)
mat<-data.frame(t,s.rm130.aee.ave.pro,
                s.rm130.aes.ave.pro,
                s.rm130.ies.ave.pro)
colnames(mat)<-c("t","aee","aes","ies")
s.rm130.gg.p<- ggplot()
s.rm130.gg.p<-s.rm130.gg.p+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
s.rm130.gg.p<-s.rm130.gg.p+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
s.rm130.gg.p<-s.rm130.gg.p+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
s.rm130.gg.p<-s.rm130.gg.p+ theme_bw()+scale_y_continuous(position = "right",limits = c(0,1))+
  theme(panel.grid=element_blank(),panel.background = element_rect(color = 'black', fill = "#F0FBFF"))+scale_x_continuous(breaks = c(5,10,15,20))+
  theme(plot.margin = unit(c(1,0,0,0),"lines"))+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
s.rm130.effect.ave.pro<-s.rm130.gg.p
s.rm130.effect.ave.pro
#---------------------ALL.E.Significant.location.average_value------------------
e.all.effect<-list()
e.all.aee.list <- list()
e.all.aes.list <- list()
e.all.ies.list <- list()
for (z in 1:dim(all.e.sig)[1]){
  e.all.miu11<-c(as.numeric(all.e.sig[z,3:8]))
  e.all.miu10<-c(as.numeric(all.e.sig[z,9:14]))
  e.all.miu00<-c(as.numeric(all.e.sig[z,15:20]))
  e.all.miu01<-c(as.numeric(all.e.sig[z,21:26]))
  
  e.all.mmiu11 <- get_miu(e.all.miu11,t)
  e.all.mmiu10 <- get_miu(e.all.miu10,t)
  e.all.mmiu00 <- get_miu(e.all.miu00,t)
  e.all.mmiu01 <- get_miu(e.all.miu01,t)
  
  e.all.aee.list[[z]]<-(e.all.mmiu11+e.all.mmiu10-e.all.mmiu01-e.all.mmiu00)/4
  e.all.aes.list[[z]]<-(e.all.mmiu11+e.all.mmiu01-e.all.mmiu10-e.all.mmiu00)/4
  e.all.ies.list[[z]]<-(e.all.mmiu11+e.all.mmiu00-e.all.mmiu10-e.all.mmiu01)/4
  
}

e.all.aee.ave<-Reduce('+',e.all.aee.list)/dim(all.e.sig)[1]
e.all.aes.ave<-Reduce('+',e.all.aes.list)/dim(all.e.sig)[1]
e.all.ies.ave<-Reduce('+',e.all.ies.list)/dim(all.e.sig)[1]

mat<-data.frame(t,e.all.aee.ave,e.all.aes.ave,e.all.ies.ave)
colnames(mat)<-c("t","aee","aes","ies")
mat[,c(2:4)] <- log10(exp(mat[,c(2:4)]))
e.all.gg<- ggplot()
e.all.gg<-e.all.gg+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
e.all.gg<-e.all.gg+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
e.all.gg<-e.all.gg+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
e.all.gg<-e.all.gg+ theme_bw()+scale_y_continuous(labels=c(expression(-10^0.25),0,expression(10^0.25),
                                                           expression(10^0.5),expression(10^0.75)),limits = c(-0.3,0.75),position = "left")+
  theme(panel.grid=element_blank())+
  theme(plot.margin = unit(c(1,0,0,0),"lines"))+
  geom_hline(aes(yintercept=0),linetype="dotted")+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
e.all.effect.ave<-e.all.gg
e.all.effect.ave
#---------------------130.E.Significant.location.average_value--------------------------
e.130.aee.ave<-Reduce('+',e.all.aee.list[be])/130
e.130.aes.ave<-Reduce('+',e.all.aes.list[be])/130
e.130.ies.ave<-Reduce('+',e.all.ies.list[be])/130
mat<-data.frame(t,e.130.aee.ave,e.130.aes.ave,e.130.ies.ave)
colnames(mat)<-c("t","aee","aes","ies")
mat[,c(2:4)] <- log10(exp(mat[,c(2:4)]))
e.130.gg<- ggplot()
e.130.gg<-e.130.gg+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
e.130.gg<-e.130.gg+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
e.130.gg<-e.130.gg+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
e.130.gg<-e.130.gg+ theme_bw()+scale_y_continuous(labels=c(expression(-10^1),expression(-10^0.5),0,
                                                           expression(10^0.5),expression(10^0.75)),
                                                  limits = c(-1,0.75),position = "left")+
  theme(panel.grid=element_blank(),panel.background = element_rect(color = 'black', fill = "#FFEBEE"))+
  theme(plot.margin = unit(c(1,0,0,0),"lines"))+
  geom_hline(aes(yintercept=0),linetype="dotted")+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
e.130.effect.ave<-e.130.gg
e.130.effect.ave
#---------------------RM130.E.Significant.location.average_value----------------
e.rm130.aee.ave<-Reduce('+',e.all.aee.list[-be])/(dim(all.e.sig)[1]-130)
e.rm130.aes.ave<-Reduce('+',e.all.aes.list[-be])/(dim(all.e.sig)[1]-130)
e.rm130.ies.ave<-Reduce('+',e.all.ies.list[-be])/(dim(all.e.sig)[1]-130)
mat<-data.frame(t,e.rm130.aee.ave,e.rm130.aes.ave,e.rm130.ies.ave)
colnames(mat)<-c("t","aee","aes","ies")
mat[,c(2:4)] <- log10(exp(mat[,c(2:4)]))
e.rm130.gg<- ggplot()
e.rm130.gg<-e.rm130.gg+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
e.rm130.gg<-e.rm130.gg+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
e.rm130.gg<-e.rm130.gg+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
e.rm130.gg<-e.rm130.gg+ theme_bw()+scale_y_continuous(labels=c(expression(-10^0.25),0,
                                                               expression(10^0.25),expression(10^0.5),expression(10^0.75)),
                                                      limits = c(-0.3,0.75),position = "left")+
  theme(panel.grid=element_blank(),panel.background = element_rect(color = 'black', fill = "#F0FBFF"))+
  theme(plot.margin = unit(c(1,0,0,0),"lines"))+
  geom_hline(aes(yintercept=0),linetype="dotted")+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
e.rm130.effect.ave<-e.rm130.gg
e.rm130.effect.ave



#---------------------ALL.E.Significant.location.variance-----------------------
load("e.var.fun.RData")
e.all.aee.var.list <- list()
e.all.aes.var.list <- list()
e.all.ies.var.list <- list()
e.all.aee.var <- c()
e.all.aes.var <- c()
e.all.ies.var <- c()
for(i in 1:dim(all.e.sig)[1]){
  e.all.aee.var.list[i] <- e.all.var.list[[i]][1]
  e.all.aes.var.list[i] <- e.all.var.list[[i]][2]
  e.all.ies.var.list[i] <- e.all.var.list[[i]][3]
}
e.all.aee.ave.var <-Reduce('+',e.all.aee.var.list)/dim(all.e.sig)[1]
e.all.aes.ave.var <-Reduce('+',e.all.aes.var.list)/dim(all.e.sig)[1]
e.all.ies.ave.var <-Reduce('+',e.all.ies.var.list)/dim(all.e.sig)[1]

mat<-data.frame(t,e.all.aee.ave.var,e.all.aes.ave.var,e.all.ies.ave.var)
mat[,c(2:4)] <- log10(exp(mat[,c(2:4)]))*dim(all.e.sig)[1]
colnames(mat)<-c("t","aee","aes","ies")
e.all.gg.v<- ggplot()
e.all.gg.v<-e.all.gg.v+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
e.all.gg.v<-e.all.gg.v+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
e.all.gg.v<-e.all.gg.v+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
e.all.gg.v<-e.all.gg.v+ theme_bw()+
  theme(panel.grid=element_blank())+
  scale_y_continuous(labels=(math_format(10^.x)),position = "left")+
  theme(plot.margin = unit(c(1,0,0,0),"lines"))+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
e.all.effect.ave.var<-e.all.gg.v
e.all.effect.ave.var
#---------------------130.E.Significant.location.variance-----------------------
e.130.aee.ave.var <-Reduce('+',e.all.aee.var.list[be])/130
e.130.aes.ave.var <-Reduce('+',e.all.aes.var.list[be])/130
e.130.ies.ave.var <-Reduce('+',e.all.ies.var.list[be])/130
mat<-data.frame(t,e.130.aee.ave.var,e.130.aes.ave.var,e.130.ies.ave.var)
mat[,c(2:4)] <- log10(exp(mat[,c(2:4)]))*130
colnames(mat)<-c("t","aee","aes","ies")
e.130.gg.v<- ggplot()
e.130.gg.v<-e.130.gg.v+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
e.130.gg.v<-e.130.gg.v+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
e.130.gg.v<-e.130.gg.v+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
e.130.gg.v<-e.130.gg.v+ theme_bw()+
  scale_y_continuous(labels=(math_format(10^.x)),position = "left")+
  theme(panel.grid=element_blank(),panel.background = element_rect(color = 'black', fill = "#FFEBEE"))+
  theme(plot.margin = unit(c(1,0,0,0),"lines"))+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
e.130.effect.ave.var<-e.130.gg.v
e.130.effect.ave.var

#---------------------RM130.E.Significant.location.variance---------------------
e.rm130.aee.ave.var <-Reduce('+',e.all.aee.var.list[-be])/(dim(all.e.sig)[1]-130)
e.rm130.aes.ave.var <-Reduce('+',e.all.aes.var.list[-be])/(dim(all.e.sig)[1]-130)
e.rm130.ies.ave.var <-Reduce('+',e.all.ies.var.list[-be])/(dim(all.e.sig)[1]-130)
mat<-data.frame(t,e.rm130.aee.ave.var,e.rm130.aes.ave.var,e.rm130.ies.ave.var)
mat[,c(2:4)] <- log10(exp(mat[,c(2:4)]))*(dim(all.e.sig)[1]-130)
colnames(mat)<-c("t","aee","aes","ies")
e.rm130.gg.v<- ggplot()
e.rm130.gg.v<-e.rm130.gg.v+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
e.rm130.gg.v<-e.rm130.gg.v+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
e.rm130.gg.v<-e.rm130.gg.v+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
e.rm130.gg.v<-e.rm130.gg.v+ theme_bw()+
  scale_y_continuous(labels=(math_format(10^.x)),position = "left")+
  theme(panel.grid=element_blank(),panel.background = element_rect(color = 'black', fill = "#F0FBFF"))+
  theme(plot.margin = unit(c(1,0,0,0),"lines"))+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
e.rm130.effect.ave.var<-e.rm130.gg.v
e.rm130.effect.ave.var


#---------------------ALL.E.Significant.location.Proportion---------------------
e.all.aee.pro.list <- list()
e.all.aes.pro.list <- list()
e.all.ies.pro.list <- list()
e.all.pro.list <- list()

for(i in 1:dim(all.e.sig)[1]){
  e.all.pro.list[[i]] <- e.all.aee.var.list[[i]]+e.all.aes.var.list[[i]]+e.all.ies.var.list[[i]]
  e.all.aee.pro.list[[i]]<-e.all.aee.var.list[[i]]/e.all.pro.list[[i]]
  e.all.aes.pro.list[[i]]<-e.all.aes.var.list[[i]]/e.all.pro.list[[i]]
  e.all.ies.pro.list[[i]]<-e.all.ies.var.list[[i]]/e.all.pro.list[[i]]
}

e.all.aee.ave.pro <-Reduce('+',e.all.aee.pro.list)/dim(all.e.sig)[1]
e.all.aes.ave.pro <-Reduce('+',e.all.aes.pro.list)/dim(all.e.sig)[1]
e.all.ies.ave.pro <-Reduce('+',e.all.ies.pro.list)/dim(all.e.sig)[1]

mat<-data.frame(t,e.all.aee.ave.pro,
                e.all.aes.ave.pro,
                e.all.ies.ave.pro)
colnames(mat)<-c("t","aee","aes","ies")
e.all.gg.p<- ggplot()
e.all.gg.p<-e.all.gg.p+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
e.all.gg.p<-e.all.gg.p+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
e.all.gg.p<-e.all.gg.p+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
e.all.gg.p<-e.all.gg.p+ theme_bw()+scale_y_continuous(position = "left",limits = c(0,1))+
  theme(panel.grid=element_blank())+
  theme(plot.margin = unit(c(1,0,0,0),"lines"))+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
e.all.effect.ave.pro<-e.all.gg.p
e.all.effect.ave.pro


#---------------------130.E.Significant.location.Proportion---------------------
e.130.aee.ave.pro <-Reduce('+',e.all.aee.pro.list[be])/130
e.130.aes.ave.pro <-Reduce('+',e.all.aes.pro.list[be])/130
e.130.ies.ave.pro <-Reduce('+',e.all.ies.pro.list[be])/130

mat<-data.frame(t,e.130.aee.ave.pro,
                e.130.aes.ave.pro,
                e.130.ies.ave.pro)
colnames(mat)<-c("t","aee","aes","ies")
e.130.gg.p<- ggplot()
e.130.gg.p<-e.130.gg.p+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
e.130.gg.p<-e.130.gg.p+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
e.130.gg.p<-e.130.gg.p+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
e.130.gg.p<-e.130.gg.p+ theme_bw()+scale_y_continuous(position = "left",limits = c(0,1))+
  theme(panel.grid=element_blank(),panel.background = element_rect(color = 'black', fill = "#FFEBEE"))+
  theme(plot.margin = unit(c(1,0,0,0),"lines"))+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
e.130.effect.ave.pro<-e.130.gg.p
e.130.effect.ave.pro
#---------------------RM130.E.Significant.location.Proportion-------------------
e.rm130.aee.ave.pro <-Reduce('+',e.all.aee.pro.list[-be])/(dim(all.e.sig)[1]-130)
e.rm130.aes.ave.pro <-Reduce('+',e.all.aes.pro.list[-be])/(dim(all.e.sig)[1]-130)
e.rm130.ies.ave.pro <-Reduce('+',e.all.ies.pro.list[-be])/(dim(all.e.sig)[1]-130)

mat<-data.frame(t,e.rm130.aee.ave.pro,
                e.rm130.aes.ave.pro,
                e.rm130.ies.ave.pro)
colnames(mat)<-c("t","aee","aes","ies")
e.rm130.gg.p<- ggplot()
e.rm130.gg.p<-e.rm130.gg.p+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
e.rm130.gg.p<-e.rm130.gg.p+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
e.rm130.gg.p<-e.rm130.gg.p+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
e.rm130.gg.p<-e.rm130.gg.p+ theme_bw()+scale_y_continuous(position = "left",limits = c(0,1))+
  theme(panel.grid=element_blank(),panel.background = element_rect(color = 'black', fill = "#F0FBFF"))+
  theme(plot.margin = unit(c(1,0,0,0),"lines"))+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
e.rm130.effect.ave.pro<-e.rm130.gg.p
e.rm130.effect.ave.pro

#---------------------S.four.SNP------------------------------------------------
#四个位点 21681/4400  57116/6194  37949/7326  1278/9104  
View(t(S_SNP[c(4400,6194,7326,9104),]))
View(t(E_SNP[c(21681,57116,37949,1278),]))

cbind(ba,be)
length(unique(be))


e.sig <- paste0(all.e.sig.loci[,1],".",all.e.sig.loci[,2])
s.sig <- paste0(all.s.sig.loci[,1],".",all.s.sig.loci[,2])
es.split <- t(as.data.frame(strsplit(intersect(e.sig,s.sig),"[.]")))
rownames(es.split) <- c(1:130)
colnames(es.split) <- c("E","S")
e.split <- t(as.data.frame(strsplit(intersect(e.sig,s.sig),"[.]")))[,1]
s.split <- t(as.data.frame(strsplit(intersect(e.sig,s.sig),"[.]")))[,2]


which(all.s.sig.loci[which(all.s.sig.loci[,2]==4400),][,1]==21681)#143
which(all.s.sig.loci[which(all.s.sig.loci[,2]==6194),][,1]==57116)#8
which(all.s.sig.loci[which(all.s.sig.loci[,2]==7326),][,1]==37949)#13
which(all.s.sig.loci[which(all.s.sig.loci[,2]==9104),][,1]==1278)#7

which(all.s.sig[,1]==rownames(S_SNP[4400,]))[143]#1206
which(all.s.sig[,1]==rownames(S_SNP[6194,]))[8]#1939
which(all.s.sig[,1]==rownames(S_SNP[7326,]))[13]#2128
which(all.s.sig[,1]==rownames(S_SNP[9104,]))[7]#2624


s.effect<-list()
for (z in c(1206,1939,2128,2624)){
  mat<-data.frame(t,s.all.aee.list[[z]],s.all.aes.list[[z]],s.all.ies.list[[z]])
  mat[,c(2:4)] <- log10(exp(mat[,c(2:4)]))
  colnames(mat)<-c("t","aee","aes","ies")
  s.gg<- ggplot(data = mat,aes(x=t,y=aee))
  s.gg<-s.gg+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
  s.gg<-s.gg+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
  s.gg<-s.gg+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
  s.gg<-s.gg+ theme_bw()+scale_y_continuous(labels =c(expression(-10^1),expression(-10^0.5),0,expression(10^0.5),expression(10^1)),
                                            limits=c(-1.1,1),position = "right")+
    theme(panel.grid=element_blank())+
    theme(plot.margin = unit(c(1,0,0,0),"lines"))+
    geom_hline(aes(yintercept=0),linetype="dotted")+scale_x_continuous(breaks = c(5,10,15,20))+
    theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
  s.effect[[z]]<-s.gg
}
effect.s<-cowplot::plot_grid(s.effect[[1206]],s.effect[[1939]],s.effect[[2128]],
                             s.effect[[2624]],s.130.effect.ave,s.rm130.effect.ave,nrow = 6)
effect.s

#---------------------S.four.SNP.variance---------------------------------------
s.variance <- list()
for(z in c(1206,1939,2128,2624)){
  mat<-data.frame(t,s.all.aee.var.list[[z]],s.all.aes.var.list[[z]],s.all.ies.var.list[[z]])
  mat[,c(2:4)] <- log10(exp(mat[,c(2:4)]))
  colnames(mat)<-c("t","aee","aes","ies")
  s.gg.v<- ggplot()
  s.gg.v<-s.gg.v+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
  s.gg.v<-s.gg.v+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
  s.gg.v<-s.gg.v+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
  s.gg.v<-s.gg.v+ 
    theme_bw()+scale_y_continuous(labels=(math_format(10^.x)),limits =c(0,1),position = "right")+
    theme(panel.grid=element_blank())+scale_x_continuous(breaks = c(5,10,15,20))+
    theme(plot.margin = unit(c(1,0,0,0),"lines"))+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
  s.variance[[z]]<-s.gg.v}

s.var.gg<-cowplot::plot_grid(s.variance[[1206]],s.variance[[1939]],
                             s.variance[[2128]],s.variance[[2624]],
                             s.130.effect.ave.var,s.rm130.effect.ave.var,nrow = 6)
s.var.gg
#---------------------S.four.SNP.Proportion-------------------------------------
s.pro <- list()
for(z in c(1206,1939,2128,2624)){
  mat<-data.frame(t,s.all.aee.pro.list[[z]],s.all.aes.pro.list[[z]],s.all.ies.pro.list[[z]])
  colnames(mat)<-c("t","aee","aes","ies")
  s.gg.p<- ggplot(data = mat,aes(x=t,y=aee))
  s.gg.p<-s.gg.p+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
  s.gg.p<-s.gg.p+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
  s.gg.p<-s.gg.p+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
  s.gg.p<-s.gg.p+ theme_bw()+scale_y_continuous(limits = c(0,1),position = "right")+
    theme(panel.grid=element_blank())+scale_x_continuous(breaks = c(5,10,15,20))+
    theme(plot.margin = unit(c(1,0,0,0),"lines"))+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
  s.pro[[z]]<-s.gg.p
}

s.pro.gg<-cowplot::plot_grid(s.pro[[1206]],s.pro[[1939]],s.pro[[2128]],
                             s.pro[[2624]],s.130.effect.ave.pro,
                             s.rm130.effect.ave.pro,nrow = 6)
s.pro.gg
#---------------------E.four.SNP------------------------------------------------
#四个位点 21681/4400 12.55752  57116/6194 12.36908 37949/7326 11.59264  1278/9104  12.03646
S_SNP.4 <- t(S_SNP[c(4400,6194,7326,9104),])
E_SNP.4 <- t(E_SNP[c(21681,57116,37949,1278),])
write.csv(S_SNP.4,"SNP4/S_SNP.4.csv")
write.csv(E_SNP.4,"SNP4/E_SNP.4.csv")
write.csv(es.130,"SNP4/es.130.csv")



which(all.e.sig.loci[which(all.e.sig.loci[,1]==21681),][,2]==4400)#114
which(all.e.sig.loci[which(all.e.sig.loci[,1]==57116),][,2]==6194)#8
which(all.e.sig.loci[which(all.e.sig.loci[,1]==37949),][,2]==7326)#17
which(all.e.sig.loci[which(all.e.sig.loci[,1]==1278),][,2]==9104)#9

which(all.e.sig[,1]==rownames(E_SNP[21681,]))[114]#40209
which(all.e.sig[,1]==rownames(E_SNP[57116,]))[8]#60277
which(all.e.sig[,1]==rownames(E_SNP[37949,]))[17]#76542
which(all.e.sig[,1]==rownames(E_SNP[1278,]))[9]#93515


e.effect<-list()
for (z in c(40209,60277,76542,93515)){
  mat<-data.frame(t,e.all.aee.list[[z]],e.all.aes.list[[z]],e.all.ies.list[[z]])
  mat[,c(2:4)] <- log10(exp(mat[,c(2:4)]))
  colnames(mat)<-c("t","aee","aes","ies")
  e.gg<- ggplot()
  e.gg<-e.gg+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
  e.gg<-e.gg+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
  e.gg<-e.gg+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
  e.gg<-e.gg+ theme_bw()+theme(panel.grid=element_blank())+
    theme(plot.margin = unit(c(1,0,0,0),"lines"))+geom_hline(aes(yintercept=0),linetype="dotted")+
    scale_y_continuous(labels=c(expression(-10^1),expression(-10^0.5),0,expression(10^0.5),expression(10^1)),
                       limits=c(-1.1,1),position = "left")+
    theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
  e.effect[[z]]<-e.gg
}
effect.e<-cowplot::plot_grid(e.effect[[40209]],e.effect[[60277]],e.effect[[76542]],
                             e.effect[[93515]],e.130.effect.ave,e.rm130.effect.ave,nrow = 6)
effect.e
#---------------------E.four.SNP.variance---------------------------------------
e.variance <- list()
for(z in c(40209,60277,76542,93515)){
  mat<-data.frame(t,e.all.aee.var.list[[z]],e.all.aes.var.list[[z]],e.all.ies.var.list[[z]])
  mat[,c(2:4)] <- log10(exp(mat[,c(2:4)]))
  colnames(mat)<-c("t","aee","aes","ies")
  e.gg.v<- ggplot(data = mat,aes(x=t,y=aee))
  e.gg.v<-e.gg.v+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
  e.gg.v<-e.gg.v+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
  e.gg.v<-e.gg.v+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
  e.gg.v<-e.gg.v+ theme_bw()+theme(panel.grid=element_blank())+
    theme(plot.margin = unit(c(1,0,0,0),"lines"))+
    scale_y_continuous(labels=(math_format(10^.x)),limits =c(0,1),position = "left")+
    theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
  e.variance[[z]]<-e.gg.v}
e.var.gg<-cowplot::plot_grid(e.variance[[40209]],e.variance[[60277]],e.variance[[76542]],
                             e.variance[[93515]],e.130.effect.ave.var,e.rm130.effect.ave.var,nrow = 6)
e.var.gg
#---------------------E.four.SNP.Proportion-------------------------------------
e.pro<-list()
for(z in c(40209,60277,76542,93515)){
  mat<-data.frame(t,e.all.aee.pro.list[[z]],e.all.aes.pro.list[[z]],e.all.ies.pro.list[[z]])
  colnames(mat)<-c("t","aee","aes","ies")
  e.gg.p<- ggplot(data = mat,aes(x=t,y=aee))
  e.gg.p<-e.gg.p+ geom_line(data = mat,aes(x=t,y=aee),colour="blue",size=1)+labs(y=NULL)
  e.gg.p<-e.gg.p+ geom_line(data = mat,aes(x=t,y=aes),colour="red",size=1)+labs(x=NULL)
  e.gg.p<-e.gg.p+ geom_line(data = mat,aes(x=t,y=ies),colour="green",size=1)
  e.gg.p<-e.gg.p+ theme_bw()+
    theme(panel.grid=element_blank())+
    theme(plot.margin = unit(c(1,0,0,0),"lines"))+ylim(0,1)+theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))
  e.pro[[z]]<-e.gg.p}
e.pro.gg<-cowplot::plot_grid(e.pro[[40209]],e.pro[[60277]],e.pro[[76542]],e.pro[[93515]],
                             e.130.effect.ave.pro,e.rm130.effect.ave.pro,nrow = 6)
e.pro.gg

fig4<-cowplot::plot_grid(effect.e,effect.s,NULL,
                          e.var.gg,s.var.gg,NULL,
                          e.pro.gg,s.pro.gg,nrow = 1,
                          rel_widths = c(1,1,0.1,1,1,0.1,1,1))  
fig4
ggsave("Fig4.pdf",fig4,width = 20,height = 14)
