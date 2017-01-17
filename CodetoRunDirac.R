---
  title: "DIRAC Code"
author: "Alison Paquette"
date: "11/16/2016"
output: html_document
---
  ##################################Finalized Code to run DIRAC####################################
#Originally written by Vineet Sangar, https://github.com/vinisan/DIRAC_Rcode/blob/master/Dirac.R#

#Load Dependencies####

library(dplyr)
library(ggplot2)
library(doParallel)
library(pheatmap)


#Load Functions####

#Note: Read pathway is set to read all pathways with between 5-100 genes.  >100 Genes greatly slows down #computational time!
readPathway<-function(path.File, exp.data){
  pathways<-readLines(path.File)
  path<-strsplit(pathways, "\t")

  return.pathways = list()
  for ( i in 1:length(path)){
    namesDIR<-as.matrix(path[[i]])

    if (length(path[[i]]) > 5 && (length(path[[i]])< 100) && length(which(rownames(exp.data) %in% namesDIR[,1])) > 5 && length(which(rownames(exp.data) %in% namesDIR[,1]))<100){
      return.pathways <-append(return.pathways, list(path[[i]]))
    }
  }

  names(return.pathways)<-sapply(return.pathways,'[',1)
  final.pathways = lapply(return.pathways,'[',-(1:2))
  return(final.pathways)
}
doPairWise<-function(pathwayOrderMatrix, gene.pairs){

  pathCond.cols<-nrow(pathwayOrderMatrix)
  class.Matrix<-matrix(, nrow = gene.pairs)

  for ( sample.Count in 1:ncol(pathwayOrderMatrix)){
    sample<-vector('numeric', gene.pairs)
    pathCond<-as.matrix(pathwayOrderMatrix[,sample.Count])
    sample<-rankVector(pathCond, pathCond.cols, gene.pairs)
    class.Matrix<-cbind(class.Matrix,sample)
  }

  class.Matrix<-class.Matrix[,-1]
  mode.Vector<-vector('numeric',nrow(class.Matrix))
  nSamples <- ncol(class.Matrix)

  for ( gene.path in 1:nrow(class.Matrix) ){
    mode.gene = 0
    mode.gene = sum(class.Matrix[gene.path,])/nSamples

    if(mode.gene < 0.5){
      mode.Vector[gene.path] = 0
    }

    else{
      mode.Vector[gene.path] = 1
    }

  }
  return(list(class.Matrix, mode.Vector))
}
rankVector<-function(pathCond, pathCond.cols, gene.pairs){
  counter= 0

  rank.compared<-vector('numeric', gene.pairs)
  for(i in 1:pathCond.cols){

    for (j in i:pathCond.cols){

      if ( i != j){
        counter= counter+1

        if(pathCond[i,1] >= pathCond[j,1]){
          rank.compared[counter] = 1

        }

        else{
          rank.compared[counter] = 0
        }

      }
    }
  }

  return(rank.compared)

}
calcRankMatching<-function(expr.Matrix,template){

  pairs<-nrow(expr.Matrix)
  rank.conservation<-vector("numeric", ncol(expr.Matrix))

  for ( columns in 1:ncol(expr.Matrix)){

    different<-which(expr.Matrix[,columns] != template)
    percentage.different = 1-(length(different)/pairs)
    rank.conservation[columns] = percentage.different

  }
  return(rank.conservation)
}
calculateAccuracy<-function(rank.difference,nCond1,nCond2){

  true.positives = 0
  ties = 0
  specificity = 0
  sensitivity = 0
  accuracy = 0
  true.positives <- length(which(rank.difference[1:nCond1] >= 0))
  sensitivity = (true.positives)/nCond1
  true.negatives = 0
  true.negatives <- length(which(rank.difference[(nCond1 +1):(nCond1 + nCond2)] < 0))
  ties = length(which(rank.difference[(nCond1 +1):(nCond1 + nCond2)] == 0))
  specificity = (true.negatives )/nCond2
  accuracy = (specificity * 0.5) + (sensitivity * 0.5)

  return(accuracy)
}
runDirac<-function(data.Cond1, data.Cond2, pathway.List, pathway.accuracy){


  for (pathwayN in 1:length(pathway.List)){
    namesDIR<-as.matrix(pathway.List[[pathwayN]])

    pathwayNdata.Cond1<-data.Cond1[which(rownames(data.Cond1) %in% namesDIR[,1]),]
    pathwayNdata.Cond2<-data.Cond2[which(rownames(data.Cond2) %in% namesDIR[,1]),]

    pathwayNDataCond1.order = apply(pathwayNdata.Cond1,2, rank)
    pathwayNDataCond2.order = apply(pathwayNdata.Cond2,2, rank)

    pathCond.cols<-nrow(pathwayNDataCond1.order)
    gene.pairs <-(pathCond.cols*(pathCond.cols-1))/2

    cond1.matrix<-matrix(, nrow = gene.pairs)
    cond2.matrix<-matrix(, nrow = gene.pairs)

    cond1.list<-doPairWise(pathwayNDataCond1.order,gene.pairs)
    cond2.list<-doPairWise(pathwayNDataCond2.order,gene.pairs)

    total.Matrix<-cbind(cond1.list[[1]],cond2.list[[1]])
    rank.matching1<-calcRankMatching(total.Matrix,cond1.list[[2]])
    rank.matching2<-calcRankMatching(total.Matrix,cond2.list[[2]])
    rank.difference = rank.matching1 - rank.matching2

    nCond1 = ncol(data.Cond1)
    nCond2 = ncol(data.Cond2)
    pathway.accuracy[pathwayN] = calculateAccuracy(rank.difference,nCond1,nCond2)
  }
  return(pathway.accuracy)
}
permuteDirac<-function(dirac.data,pathway.List, nCond1, nCond2, pathway.Difference, calculated.Accuracy){

  for (run in 1:permutations){

    print (run)
    permuted.Data<-dirac.data[,sample(ncol(dirac.data), replace = TRUE)]
    data.Cond1<-permuted.Data[,1:nCond1]
    data.Cond2<-permuted.Data[,(nCond1 +1):(nCond1 + nCond2)]

    pathway.Accuracy<-vector('numeric',length(pathway.List))
    permuted.Difference <- runDirac(data.Cond1, data.Cond2,pathway.List, pathway.Accuracy)
    more.Difference<-which(permuted.Difference >= pathway.Difference)

    if (length(more.Difference) > 0)
    {
      calculated.Accuracy[more.Difference]<-(calculated.Accuracy[more.Difference] + 1)
    }
  }
  return(calculated.Accuracy)
}
#Note: do Permutations is corrected from vineets code
doPermutations<-function(expr.data, pathway.list, nCond1, nCond2,  PathwayAccuracy, calculated.Accuracy, CVaccuracy, permutations, cores){

  calculated.Accuracy<-vector('numeric',length(pathway.list))
  registerDoParallel(cores)

  parallelPermuted <- foreach(i=1:cores, .combine='cbind') %dopar% permuteDirac(expr.data, pathway.list, nCond1, nCond2,  PathwayAccuracy, calculated.Accuracy)
  PermutedAccuracy =apply(parallelPermuted, 1, sum)
  Pvalue<-PermutedAccuracy/(permutations*cores)
  q<-p.adjust(Pvalue,method = "BH", n = length(Pvalue))


  #Storing the values in a data matrix
  results<-cbind(names(pathway.list), PathwayAccuracy)
  results<-cbind(results, as.numeric(Pvalue))
  results<-cbind(results, as.numeric(q))
  results<-cbind(results, as.numeric(CVaccuracy))

  return(results)
}
#Note: CV accuracy is corrected
doCrossValidation<-function(cond1.data,cond2.data,CV.pathwayList){

  CVaccuracy<-vector('numeric',length(CV.pathwayList) )
  index<-vector('character',length = (ncol(cond1.data) + ncol(cond2.data)))
  index[1:ncol(cond1.data)] = "1"
  index[(ncol(cond1.data)+1):(ncol(cond1.data) + ncol(cond2.data))] = "0"
  #print(index)

  CVaccuracy[1:length(CV.pathwayList)] = 0
  rank.difference<-vector('numeric',length(CV.pathwayList))
  dirac.data<-cbind(cond1.data,cond2.data)
  #dim(dirac.data)
  for (samples in 1:ncol(dirac.data)){

    print(samples)
    pathway.accuracy<-vector('numeric',length(CV.pathwayList))
    cvData<-dirac.data[,-samples]
    cvIndex<-index[-samples]
    holdData<-as.matrix(dirac.data[,samples])
    holdIndex<-index[samples]

    rownames(holdData)<-rownames(dirac.data)
    #print(dim(holdData))
    #print(holdIndex)
    data1<-which(cvIndex == "1")
    data2<-which(cvIndex == "0")
    #print(data1)
    data.Cond1<-cvData[,data1]
    data.Cond2<-cvData[,data2]
    cvPathway.Accuracy<-vector('numeric',length(CV.pathwayList))

    nCond1<-ncol(data.Cond1)
    nCond2<-ncol(data.Cond2)
    #print(dim(dirac.data))
    for (pathwayN in 1:length(CV.pathwayList)){

      namesDIR<-as.matrix(CV.pathwayList[[pathwayN]])

      pathwayNdata.Cond1<-data.Cond1[which(rownames(data.Cond1) %in% namesDIR[,1]),]
      pathwayNdata.Cond2<-data.Cond2[which(rownames(data.Cond2) %in% namesDIR[,1]),]
      #print(dim(pathwayNdata.Cond1))
      #print(which(rownames(holdData) %in% namesDIR[,1]))
      holdDataGenes<-holdData[which(rownames(holdData) %in% namesDIR[,1]),]
      #print(holdDataGenes)
      pathwayNDataCond1.order = apply(pathwayNdata.Cond1,2, rank)
      pathwayNDataCond2.order = apply(pathwayNdata.Cond2,2, rank)
      holdDataGenes.order     = as.matrix(rank(holdDataGenes))

      pathCond.cols<-nrow(pathwayNDataCond1.order)
      #print(pathCond.cols)

      gene.pairs <-(pathCond.cols*(pathCond.cols-1))/2



      cond1.matrix<-matrix(, nrow = gene.pairs)
      cond2.matrix<-matrix(, nrow = gene.pairs)
      hold.matrix<-matrix(, nrow = gene.pairs)

      #pathCond<-as.matrix(holdDataGenes.order[,sample.Count])
      cond1.list<-doPairWise(pathwayNDataCond1.order,gene.pairs)
      cond2.list<-doPairWise(pathwayNDataCond2.order,gene.pairs)
      #print(holdDataGenes.order)
      hold.list<-rankVector(holdDataGenes.order, nrow(holdDataGenes.order), gene.pairs)

      #print(hold.list)
      total.Matrix<-cbind(cond1.list[[1]],cond2.list[[1]])
      template1<-cond1.list[[2]]
      template2<-cond2.list[[2]]
      #print(which(hold.list != template1))
      #print(which(hold.list != template2))

      different1<-which(hold.list != template1)
      percentage.different1 = 1-(length(different1)/gene.pairs)
      rank.conservation1 = percentage.different1

      different2<-which(hold.list != template2)
      percentage.different2= 1-(length(different2)/gene.pairs)
      rank.conservation2 = percentage.different2
      print(rank.conservation1)
      print(rank.conservation2)
      rank.difference = rank.conservation1 - rank.conservation2
      print(rank.difference)
      if(rank.difference > 0 & holdIndex == 1)
      {
        CVaccuracy[pathwayN] = CVaccuracy[pathwayN] +1
      }
      else if(rank.difference < 0 & holdIndex == 0)

      {
        CVaccuracy[pathwayN] = CVaccuracy[pathwayN] +1
      }
      else if(rank.difference == 0 )

      {
        CVaccuracy[pathwayN] = CVaccuracy[pathwayN] + 0.5
      }
      #rank.difference[pathwayN] = rank.conservation1 - rank.conservation2
      #print(rank.difference)
      #pathway.accuracy[pathwayN] = calculateAccuracy(rank.difference,nCond1,nCond2)
      #print(CVaccuracy)

    }
    #print(pathway.accuracy)
  }
  print((CVaccuracy/(ncol(dirac.data))))
  return((CVaccuracy/(ncol(dirac.data))))
}


####RUN DIRAC ON GENE EXPRESSION MATRIX####

#Load Normalized Gene expression file#
#note: first row is binary vector of phenotype of interest

load("~/DATA/282015Option4forDIRAC.RData")

#Check # of Cases and Controls
TindexDir<-DATA[1,]
table(TindexDir)

Texpr.data<-DATA[-1,]

#Format Data to be run
cond1.DIR<-which(TindexDir == 0)
cond2.DIR<-which(TindexDir == 1)

Tdata.Cond1<-Texpr.data[,cond1.DIR]
Tdata.Cond2<-Texpr.data[,cond2.DIR]
Tdirac.data<-cbind(Tdata.Cond1,Tdata.Cond2)

nCond1<-ncol(Tdata.Cond1)
nCond2<-ncol(Tdata.Cond2)


#Create Network List####
#The networks are different pathways archived for GSEA downloaded from http://software.broadinstitute.org/gsea/msigdb

path.File<-"~/GSEA msdb/GSEAentrez2kegg.gmt"
Tpathway.list<-readPathway(path.File, Tdirac.data)
print(length(Tpathway.list)) #how many pathways are you counting?


#Perform DIRAC & LOOCV###
#Note: This is a very computationally intensive and time consumng step that is NOT parellelized!


Tpathway.Accuracy<-vector('numeric',length(Tpathway.list))
TPathwayAccuracy<-runDirac(Tdata.Cond1, Tdata.Cond2, Tpathway.list,  Tpathway.Accuracy)
CVaccuracy<-doCrossValidation(Tdata.Cond1, Tdata.Cond2,Tpathway.list)

#Perform Permutations to Obtain P Value

set.seed(32)#random number, to gt reproducible results

#Need to do at least 1000 Permutations.
permutations = 50
cores= 20
storeResultsT<-doPermutations(Tdirac.data, Tpathway.list, nCond1, nCond2,  TPathwayAccuracy, Tcalculated.Accuracy, CVaccuracy, permutations, cores)

colnames(storeResultsT)<-c("Pathway","PathwayAccuracy","PValue","BH Q Value","CV ACC")

#save as either CSV or RData file
save(storeResultsT,file="~/Results/PTBREACTOME1_3_2017.Rdata")




