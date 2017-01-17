##############################################################
#######    Visualizing DIRAC Output & Interpretation     #####
########       Alison G Paquette, January 11 2017  ##########
##############################################################

#Part 1: Read in Data
Data<-read.csv("~/Documents/Project 1. Preterm Birth/PaperPreprocessingJan/DIRAC/GeneOntology_DIRAC_1_11_2017.csv")

#Format Data
rownames(Data)<-Data$Pathway
Data<-Data[order(Data$BH.Q.Value),] #Sort by P value

Data<-Data[,c(5:6)]
#Create Barplot of top 100 Pathways
Data<-Data[1:100,]
Data<-Data[order(Data$CV.ACC,decreasing=T),]
par(mar=c(5,20,4,2))

mycolours <- as.character(cut(Data$BH.Q.Value,c(1,0.05,0.01,0.005,0), right=F, labels=c("azure3","yellow","orange","red"), include.lowest=T))

barplot(Data$CV.ACC, main="Top 100 Reactome Pathways", horiz=TRUE,col=mycolours,
        names.arg=rownames(Data),cex.names=0.3,las=2,
        xlim=c(0,1)
)

abline(v=0.7,col="black",lty=2,lwd=3)
legend("topright",legend=c("Not significant","<0.05","<0.01","<0.005"),pch=16,col=c("azure3","yellow","orange","red"))
