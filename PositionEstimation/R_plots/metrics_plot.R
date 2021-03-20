if(!require(plotly)){
  install.packages('plotly')
  library(plotly)
}
if(!require(ggplot2)){
  install.packages('ggplot2')
  library(ggplot2) 
}
if(!require(dplyr)){
  install.packages('dplyr')
  library(dplyr)
}
if(!require(plyr)){
  install.packages('plyr')
  library(plyr)
}
if(!require(R.matlab)){
  install.packages('R.matlab')
  library(R.matlab)
}
if(!require(yarrr)){
  install.packages('yarrr')
  library(yarrr)
}
if(!require(grid)){
  install.packages('grid')
  library(grid)
}
if(!require(reshape2)){
  install.packages('reshape2')
  library(reshape2)
}
if(!require(pracma)){
  install.packages('pracma')
  library(pracma)
}
if(!require(MASS)){
  install.packages('MASS')
  library(MASS)
}
if(!require(PRROC)){
  install.packages('PRROC')
  library(PRROC)
}
if(!require(emoa)){
  install.packages('emoa')
  library(emoa)
}

if(!require(readr)){
  install.packages("readr")
  library(readr)
}
if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}
if(!require(multcomp)){
  install.packages("multcomp")
  library(multcomp)
}

if(!require(car)){
  install.packages("car")
  library(car)
}
if(!require(knitr)){
  install.packages("knitr")
  library(knitr)
}
if(!require(agricolae)){
  install.packages("agricolae")
  library(agricolae)
}

if(!require(plotly)){
  install.packages('plotly')
  library(plotly)
}
if(!require(ggplot2)){
  install.packages('ggplot2')
  library(ggplot2)
}
if(!require(dplyr)){
  install.packages('dplyr')
  library(dplyr)
}
if(!require(R.matlab)){
  install.packages('R.matlab')
  library(R.matlab)
}

if(!require(grid)){
  install.packages('grid')
  library(grid)
}
if(!require(reshape2)){
  install.packages('reshape2')
  library(reshape2)
}
if(!require(PRROC)){
  install.packages('PRROC')
  library(PRROC)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

getwd()


# pasta onde ficam alocados os resultados obtidos da simulacao
pasta = 'resultados1'

# Coleta de erros
#ERRO_caso_1 <- data.frame(readMat(paste(pasta,'/ERR_TOTAL_IMG_Rotacionada_spre.mat', sep="")))
#colnames(ERRO_caso_1)<-c("ERRO_caso_1")

#ERRO_caso_2 <- data.frame(readMat(paste(pasta,'/ERR_TOTAL_IMG_Rotacionada_cpre.mat', sep="")))
#colnames(ERRO_caso_2)<-c("ERRO_caso_2")

#ERRO_caso_3 <- data.frame(readMat(paste(pasta,'/ERR_TOTAL_IMG_Rotacionada_spre_bnorm.mat', sep="")))
#colnames(ERRO_caso_3)<-c("ERRO_caso_3")

#ERRO_caso_4 <- data.frame(readMat(paste(pasta,'/ERR_TOTAL_IMG_Rotacionada_cpre_badap.mat', sep="")))
#colnames(ERRO_caso_4)<-c("ERRO_caso_4")

#ERRO_caso_5 <- data.frame(readMat(paste(pasta,'/ERR_TOTAL_IMG_ransac_matlab_bnormal.mat', sep="")))
#colnames(ERRO_caso_5)<-c("ERRO_caso_5")

#ERRO_caso_6 <- data.frame(readMat(paste(pasta,'/ERR_TOTAL_IMG_ransac_matlab_badap.mat', sep="")))
#colnames(ERRO_caso_6)<-c("ERRO_caso_6")

ERRO_caso_7 <- data.frame(readMat(paste(pasta,'/ERR_TOTAL_IMG_ransac_opencv_badap.mat', sep="")))
colnames(ERRO_caso_7)<-c("ERRO_caso_7")

x <- c(1:10)

#ERROR<- data.frame(x,ERRO_caso_1,ERRO_caso_2,ERRO_caso_3,ERRO_caso_4,ERRO_caso_5,ERRO_caso_6,ERRO_caso_7)
ERROR<- data.frame(x,ERRO_caso_7)
# colnames(ERROR) <- c('x','Erro_caso_1', 'Erro_caso_2', 'Erro_caso_3', 'Erro_caso_4', 'Erro_caso_5','Erro_caso_6','Erro_caso_7')

# calculo de la media e del desvio padrao SD
ERROR_1 <- data.frame(ERRO_caso_7)
med_ERROR_1 <- data.frame(apply(ERROR_1,2,mean,na.rm=T))
sd_ERROR_1 <- data.frame(apply(ERROR_1,2,sd,na.rm=T))

#plot de una rota individual
p1<-ggplot(ERROR, aes(x=x,y=ERRO_caso_7, fill=ERRO_caso_7))+
  geom_point(aes(color = ERRO_caso_7)) +
  # geom_line()+
  geom_smooth(alpha=.7)+
  geom_bar(position=position_dodge(), stat="identity",alpha=.4) +
  # geom_linerange() +
  xlab("Images") + ylab("EPE (meters)") + ggtitle("EPE for Case 1")
p1


##### ajuste de los datos para plotar todos en un solo plot
ERROR_11<-data.frame(ERRO_caso_1,ERRO_caso_2,ERRO_caso_3,ERRO_caso_4,ERRO_caso_5,ERRO_caso_6,ERRO_caso_7)
colnames(ERROR_11) <- c('EPE case 1', 'EPE case 2', 'EPE case 3', 'EPE case 4', 'EPE case 5',
                        'EPE case 6','EPE case 7')
ERROR_2 <- melt(ERROR_11)
xx=(rep(1:1593,7))
ERROR_22 <- data.frame(ERROR_2,xx)
ERROR_22a<-(ERROR_22[1:6372,])
p22a<-ggplot(ERROR_22a, aes(x=xx,y=value, fill=factor(variable)))+
  geom_point(aes(color = variable),alpha=.3) +
  geom_smooth(alpha=.5)+
  # geom_bar(position=position_dodge(), stat="identity",alpha=.4) +
  xlab("Images") + ylab("EPE (meters)") + ggtitle("EPE for Cases 1 - 4")
p22a
p22a + geom_smooth(alpha=.9)+facet_wrap(~variable)+geom_bar(position=position_dodge(), 
                                                            stat="identity",alpha=.3)

ERROR_22b<-(ERROR_22[4780:length(ERROR_22[,1]),])
p22b<-ggplot(ERROR_22b, aes(x=xx,y=value, fill=factor(variable)))+
  geom_point(aes(color = variable),alpha=.3) +
  geom_smooth(alpha=.5)+
  # geom_bar(position=position_dodge(), stat="identity",alpha=.4) +
  xlab("Images") + ylab("EPE (meters)") + ggtitle("EPE for Cases 4 - 7")
p22b
p22b + geom_smooth(alpha=.9)+facet_wrap(~variable)+geom_bar(position=position_dodge(), 
                                                            stat="identity",alpha=.3)


p22c<-ggplot(ERROR_22, aes(x=xx,y=value, fill=factor(variable)))+
  geom_point(aes(color = variable),alpha=.3) +
  geom_smooth(alpha=.5)+
  # geom_bar(position=position_dodge(), stat="identity",alpha=.4) +
  xlab("Images") + ylab("EPE (meters)") + ggtitle("EPE for all Cases")+
  labs(color='Cases') +labs(fill='Cases')
p22c
p22c + geom_smooth(alpha=.9)+facet_wrap(~variable)+geom_bar(position=position_dodge(), 
                                                            stat="identity",alpha=.3)+
  theme(legend.justification=c(1,0), legend.position=c(1,0), legend.box = "horizontal",
        plot.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(1.1)),
        axis.text.x=element_text(hjust=1,vjust = 1),
        axis.title.y = element_text(size = rel(1.1)),
        axis.title.x = element_text(size = rel(1.1)),
        legend.title=element_text(size = rel(1.2)),
        legend.text=element_text(size = rel(1)),
        legend.key.size = unit(1.2, 'lines'))+
  labs(color='Cases') +labs(fill='Cases') 




#plot de densidades en los erros
ERROR_11<-data.frame(ERRO_caso_1,ERRO_caso_2,ERRO_caso_3,ERRO_caso_4,ERRO_caso_5,ERRO_caso_6,ERRO_caso_7)
colnames(ERROR_11) <- c('EPE case 1', 'EPE case 2', 'EPE case 3', 'EPE case 4', 'EPE case 5',
                        'EPE case 6','EPE case 7')

ERROR_2 <- melt(ERROR_11)
densi <- ggplot(ERROR_2, aes(x= value,fill=factor(variable)))+
  geom_density(alpha=0.55)+
  # geom_vline(data=med_ERROR_11, aes(xintercept=Means, color=Cases),
  #            linetype="dashed", size=1)+
  labs(title="Density EPE plot for all cases", 
       # subtitle="City Mileage Grouped by Number of cylinders",
       x="EPE (meters)",
       y="Density",
       fill="# Cases")
densi


med_ERROR_11 <- data.frame(c('EPE case 1', 'EPE case 2', 'EPE case 3', 'EPE case 4', 'EPE case 5',
                             'EPE case 6','EPE case 7'),med_ERROR_1)
colnames(med_ERROR_11) <- c('Cases','Means')
med_ERROR_11 <- data.frame(med_ERROR_11)
densi + geom_density(alpha=0.55)+
  # geom_vline(aes(xintercept=Means, color=factor(Cases)),med_ERROR_11,linetype="dashed", size=1)+
  facet_wrap(~variable)+
  geom_vline(data=filter(ERROR_22, variable=="EPE case 1"), aes(xintercept=med_ERROR_11$Means[1]), colour="black",
             linetype="dashed", size=1) +
  geom_text(data=filter(ERROR_22, variable=="EPE case 1"), aes(x=med_ERROR_11$Means[1], y=0.045),label=round((med_ERROR_11$Means[1]),digits=2),hjust=-0.3, size=6)+
  geom_vline(data=filter(ERROR_22, variable=="EPE case 2"), aes(xintercept=med_ERROR_11$Means[2]), colour="black",
             linetype="dashed", size=1) +
  geom_text(data=filter(ERROR_22, variable=="EPE case 2"), aes(x=med_ERROR_11$Means[2], y=0.045),label=round((med_ERROR_11$Means[2]),digits=2),hjust=-0.3, size=6)+
  geom_vline(data=filter(ERROR_22, variable=="EPE case 3"), aes(xintercept=med_ERROR_11$Means[3]), colour="black",
             linetype="dashed", size=1) +
  geom_text(data=filter(ERROR_22, variable=="EPE case 3"), aes(x=med_ERROR_11$Means[3], y=0.045),label=round((med_ERROR_11$Means[3]),digits=2),hjust=-0.3, size=6)+
  geom_vline(data=filter(ERROR_22, variable=="EPE case 4"), aes(xintercept=med_ERROR_11$Means[4]), colour="black",
             linetype="dashed", size=1) +
  geom_text(data=filter(ERROR_22, variable=="EPE case 4"), aes(x=med_ERROR_11$Means[4], y=0.045),label=round((med_ERROR_11$Means[4]),digits=2),hjust=-0.3, size=6)+
  geom_vline(data=filter(ERROR_22, variable=="EPE case 5"), aes(xintercept=med_ERROR_11$Means[5]), colour="black",
             linetype="dashed", size=1) +
  geom_text(data=filter(ERROR_22, variable=="EPE case 5"), aes(x=med_ERROR_11$Means[5], y=0.045),label=round((med_ERROR_11$Means[5]),digits=2),hjust=-0.3, size=6)+
  geom_vline(data=filter(ERROR_22, variable=="EPE case 6"), aes(xintercept=med_ERROR_11$Means[6]), colour="black",
             linetype="dashed", size=1) +
  geom_text(data=filter(ERROR_22, variable=="EPE case 6"), aes(x=med_ERROR_11$Means[6], y=0.045),label=round((med_ERROR_11$Means[6]),digits=2),hjust=-0.3, size=6)+
  geom_vline(data=filter(ERROR_22, variable=="EPE case 7"), aes(xintercept=med_ERROR_11$Means[7]), colour="black",
             linetype="dashed", size=1) +
  geom_text(data=filter(ERROR_22, variable=="EPE case 7"), aes(x=med_ERROR_11$Means[7], y=0.045),label=round((med_ERROR_11$Means[7]),digits=2),hjust=-0.3, size=6)+
  theme(legend.justification=c(1,0), legend.position=c(1,0), 
        legend.box = "horizontal",
        plot.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(1.1)),
        axis.text.x=element_text(hjust=1,vjust = 1),
        axis.title.y = element_text(size = rel(1.1)),
        axis.title.x = element_text(size = rel(1.1)),
        legend.title=element_text(size = rel(1.2)),
        legend.text=element_text(size = rel(1)),
        legend.key.size = unit(1.2, 'lines'))+
  labs(color='Cases') +labs(fill='Cases') 
