#-----------------------------------requiredPackages----------------------------
requiredPackages = c("mvtnorm","ggnewscale","VennDiagram","dplyr","parallel","grid","pbapply","futile.logger","plyr","scales","patchwork","ggplot2","deSolve","Rmisc","glmnet","gridExtra","tidyverse","ggpubr")
for(packages in requiredPackages){
  if(!require(packages,character.only = TRUE)) install.packages(packages)
  require(packages,character.only = TRUE)
}
#-----------------------------------inputdata-----------------------------------
all.e.sig<-read.csv("all.e.sig.csv",row.names=1)
all.e.sig.loci<-read.csv("loci.ee.0.01.csv",row.names=1)
all.s.sig<-read.csv("all.s.sig.csv",row.names=1)
all.s.sig.loci<-read.csv("loci.s.0.01.csv",row.names=1)
S_SNP <- read.table("S-SNP.txt",row.names = 1,header=T)
E_SNP <- read.table("E-SNP.txt",row.names = 1,header=T)
#-----------------------------------Data preprocessing-------------------------- 
es.co <- rbind(all.e.sig.loci[,c(1,2)],all.s.sig.loci[,c(1,2)])
se.co <- rbind(all.s.sig.loci[,c(1,2)],all.e.sig.loci[,c(1,2)])
es.130 <- all.s.sig.loci[which(duplicated(es.co))-132535,]#S表型下的
se.130 <- all.e.sig.loci[which(duplicated(se.co))-3543,]#E表型下的



e <- read.csv("loci.ee.0.01.csv",row.names = 1,header=T)
s <- read.csv("loci.s.0.01.csv",row.names = 1,header=T)
colnames(e) <- c("E_SNP","S_SNP","P")
colnames(s) <- c("E_SNP","S_SNP","P")
e[,1] <- as.numeric(rownames(E_SNP)[e[,1]])
e[,2] <- as.numeric(rownames(S_SNP)[e[,2]])

s[,1] <- as.numeric(rownames(E_SNP)[s[,1]])
s[,2] <- as.numeric(rownames(S_SNP)[s[,2]])

se.130.loci <- se.130
colnames(se.130.loci) <- c("E_SNP","S_SNP","P")
se.130.loci[,1] <- as.numeric(rownames(E_SNP)[se.130.loci[,1]])
se.130.loci[,2] <- as.numeric(rownames(S_SNP)[se.130.loci[,2]])
which(duplicated(rbind(se.130.loci[,c(1:2)],e[,c(1:2)])))

#-----------------------------------Fig3.A--------------------------------------
g <- ggplot(e,aes(y=E_SNP,x=S_SNP,colour = P))+
  geom_point(shape=21,alpha=0.4)+
  theme_bw()+scale_color_gradient(limits=c(10,22.5),low = "#1E90FF", high = "#000080")+
  labs(colour=NULL)+ggtitle("E") +
  theme(plot.title = element_text(hjust = 0.5,size=18))+
  theme(panel.grid=element_blank(),
        axis.text=element_text(color='black'))+
  labs(x=NULL)+labs(y=NULL)+theme(legend.position="none")+
  scale_x_continuous(expand = c(0.02,0.02),breaks=(c(500000,1000000,1500000,2000000)),labels = c(500,1000,1500,2000))+
  scale_y_continuous(expand = c(0.02,0.02),breaks=(c(0,1000000,2000000,3000000,4000000)),labels = c(0,1000,2000,3000,4000))+
  theme(axis.text.x=element_text(size=14))+theme(axis.text.y=element_text(size=14))

g <- g+new_scale_colour()+geom_point(data=se.130.loci,aes(y=E_SNP,x=S_SNP,colour = P),shape=21,alpha=1,size=2)+
  scale_color_gradient(limits=c(10,16),low = "red", high = "red")

g  

es.130.loci <- es.130
colnames(es.130.loci) <- c("E_SNP","S_SNP","P")
es.130.loci[,1] <- as.numeric(rownames(E_SNP)[es.130.loci[,1]])
es.130.loci[,2] <- as.numeric(rownames(S_SNP)[es.130.loci[,2]])

gb <- ggplot(s,aes(y=E_SNP,x=S_SNP,colour=P))+
  geom_point(shape=21,alpha=0.4)+
  theme_bw()+
  scale_color_gradient(limits=c(10,22.5),low = "#1E90FF", high = "#000080")+
  labs(colour="-log(P)")+ggtitle("S")+
  theme(plot.title = element_text(hjust = 0.5,size=18))+
  labs(x=NULL)+labs(y=NULL)+
  scale_x_continuous(expand = c(0.02,0.02),breaks=(c(500000,1000000,1500000,2000000)),labels = c(500,1000,1500,2000))+
  scale_y_continuous(expand = c(0.02,0.02),breaks=(c(0,1000000,2000000,3000000,4000000)),labels = c(0,1000,2000,3000,4000))+
  theme(axis.text.x=element_text(size=14,colour = "black"))+theme(axis.text.y=element_text(size=14,colour = "black"))
gb <- gb+new_scale_colour()+geom_point(data=es.130.loci,aes(y=E_SNP,x=S_SNP,colour = P),shape=21,alpha=1,size=2)+
  scale_color_gradient(limits=c(10,16),low = "red", high = "red")



gb

genome<-cowplot::plot_grid(g, gb, nrow = 1) #labels = LETTERS[1:2]
genome <- annotate_figure(genome,
                          left = text_grob("E. coil Genome (Mb)",color = "black",size = 20,rot = 90),
                          bottom = text_grob("S. aureus Genome (Mb)", color = "black",size = 20))

genome
#ggsave("AB.pdf",genome,width = 9,height = 6)

#-----------------------------------Fig3.B--------------------------------------
loci.ee.0.01<-read.csv("loci.ee.0.01.csv",header = T,row.names = 1)
loci.s.0.01<-read.csv("loci.s.0.01.csv",header = T,row.names = 1)
dim(inner_join(loci.ee.0.01[,1:2],loci.s.0.01[,1:2]))
dim(loci.ee.0.01)[1]
dim(loci.s.0.01)[1]
ven <- draw.pairwise.venn(132535,3543,130,c("E", "S"),
                          scaled = FALSE,
                          lty = "blank",  
                          fill = c("#d8986e","#d53f89"),filename=NULL,
                          cex = 2,cat.cex = 2)
#-----------------------------------Fig3.C--------------------------------------
ven2 <- draw.pairwise.venn(2251, 285, 255,
                           scaled = FALSE,
                           lty = "blank",  
                           fill = c("#41C6F6","#7141F6"),filename=E,
                           cex = 2,cat.cex = 2)

ven3 <- draw.pairwise.venn(1173, 322, 294,
                           scaled = FALSE,
                           lty = "blank",  
                           fill = c("#f6416c","#00b8a9"),filename=E,
                           cex = 2,cat.cex = 2)
grid.newpage()  # new page
pushViewport(viewport(layout = grid.layout(2,1)))
print(genome, vp=viewport(layout.pos.row=1,layout.pos.col=1))
pushViewport(viewport(layout.pos.row=2,layout.pos.col=1))
grid.draw(ven)#14*12

grid.newpage()
pushViewport(viewport(layout = grid.layout(1,1)))
pushViewport(viewport(layout.pos.row=1,layout.pos.col=1))
grid.draw(ven2)#3*5

grid.newpage()
pushViewport(viewport(layout = grid.layout(1,1)))
pushViewport(viewport(layout.pos.row=1,layout.pos.col=1))
grid.draw(ven3)
