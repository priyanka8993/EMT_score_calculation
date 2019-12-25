
library(readxl)
library(matlabr)

EMT76GS = function(finalGSEMat,gseID){
  cat("Calculating EMT Score by using 76 gene signatures ....\n")

  remIdx = which(apply(finalGSEMat,1,function(x) any(x == NaN | x == -Inf)) == TRUE)
  if(length(remIdx) > 0 ) finalGSEMat = finalGSEMat[-remIdx, ]

  finalGSEMat[is.na(finalGSEMat)] = 0

  sampleNum = ncol(finalGSEMat)
  genes = finalGSEMat[,1]
  exp = apply(finalGSEMat[ ,3:sampleNum],2,as.numeric)

  EMTSignature = data.frame(read_excel("./Gene_Signatures/76GS/EMT_signature_final.xlsx",col_names = TRUE))
  EMTIdx = match(unique(na.omit(EMTSignature[,2])),genes)
  geneFound = length(na.omit(EMTIdx))
  cat(paste(geneFound ,"gene's expression values found \n",sep = " "))

  EMTMat = exp[na.omit(EMTIdx),]
  row.names(EMTMat) = genes[na.omit(EMTIdx)]

  ## get the weights for each gene
  ecadhExp = grep("^CDH1$",row.names(EMTMat))
  if(length(ecadhExp) == 0 ){
    cat("CDH1 gene not found- 76 GS EMT score cannot be calculated\n") 
    EMTScoreStd = rep(0, ncol(exp))
  } else{
    ecadhExpVal = EMTMat[ecadhExp, ]
    weightVal = apply(EMTMat,1,function(x) cor(x,ecadhExpVal,use = "complete.obs"))
    EMTMatWt = weightVal*EMTMat
    EMTScore = apply(EMTMatWt,2,sum)
    EMTScoreMean = mean(EMTScore)
    EMTScoreStd = EMTScore-EMTScoreMean

  }  
  outFile = paste("./",gseID,"/",paste(gseID, "_EMT_Score_76GS.txt",sep=""),sep = "")
  if(!file.exists(outFile)) file.create(outFile)
  #write(gseID,file =outFile,append = TRUE)
  emtWrite = cbind(colnames(exp),EMTScoreStd)
  #write(colnames(exp),file = outFile,sep = '\t',ncolumns= length(colnames(exp)),append = TRUE)
  write.table(emtWrite,file = outFile,sep = '\t',row.names = F,quote = F)

  return(EMTScoreStd)

}


KSScore = function(expMat,gseID,outDirectory){

  cat("Calculating EMT Score by KS score method ....")
  remIdx = which(apply(expMat,1,function(x) any(x == NaN | x == -Inf)) == TRUE)
  if(length(remIdx) > 0 ) expMat = expMat[-remIdx, ]
  sampleNum = ncol(expMat)
  genes = expMat[,1]
  exp = apply(expMat[ ,3:sampleNum],2,as.numeric)
  
  EMTSignature = data.frame(read_excel("./Gene_Signatures/KS/EM_gene_signature_cellLine.xlsx",col_names = FALSE))
  commonSig = intersect(EMTSignature[,1],genes)
  EMTExpIdx = match(commonSig,genes)
  EMTExp = exp[EMTExpIdx, ]
  EMTGrpIdx = match(commonSig,EMTSignature[,1])
  geneCat = EMTSignature[EMTGrpIdx,2]
  epiIdx = which(geneCat == "Epi")
  mesIdx = which(geneCat == "Mes")

    ## Perform KS test
    sampleScore2 = matrix(0,nrow=ncol(EMTExp),ncol=6)
    rownames(sampleScore2) = colnames(EMTExp)
    for(i in 1:ncol(EMTExp)){
        ## Two sided test
        ksTwoSided =  ks.test(EMTExp[mesIdx,i],EMTExp[epiIdx,i])
        ## One sided test: ecdf(Mes) > ecdf(Epi)
        ksResGrt = ks.test(EMTExp[mesIdx,i],EMTExp[epiIdx,i],alternative = "greater")
        ## One sided test: ecdf(Epi) > ecdf(Mes)
        ksResLess = ks.test(EMTExp[epiIdx,i],EMTExp[mesIdx,i],alternative = "greater")
        sampleScore2[i, ] = c(ksTwoSided$statistic,ksTwoSided$p.value,
                  ksResGrt$statistic,ksResGrt$p.value,
                  ksResLess$statistic,ksResLess$p.value)
    }

    ## Assign signs to EMT score of sample based on test statistic and pvalue
    finalScore = matrix(0,nrow = nrow(sampleScore2),ncol = 1)
    for (i in 1:nrow(sampleScore2)) {
        if(sampleScore2[i,4] < 0.05){
          finalScore[i, ] = c(-1 * sampleScore2[i,3])
        } else if (sampleScore2[i,6] < 0.05){
          finalScore[i, ] = sampleScore2[i,5]
        } else {

          if(sampleScore2[i,5] == max(c(sampleScore2[i,3],sampleScore2[i,5]))){
            finalScore[i, ] = max(c(sampleScore2[i,3],sampleScore2[i,5]))
          } else {
            finalScore[i, ] = (-1 * max(c(sampleScore2[i,3],sampleScore2[i,5])))
          }
          
      }
  }

  outputFile = paste("./",gseID,"/",paste(gseID,"_EMT_KS_Score_GSE.txt",sep = ''),sep = "")
  if(!file.exists(outputFile)) file.create(outputFile)

  #write(gseID,file =outputFile,append = TRUE)
  ksOut = cbind(colnames(EMTExp),finalScore)
  #write(colnames(EMTExp),file = outputFile,sep = '\t',ncolumns= length(colnames(EMTExp)),append = TRUE)
  write.table(ksOut,file = outputFile,sep = '\t',row.names=F,quote= F)
  return(finalScore)

}



getEMTScore = function(expVal,gseNum,gplID){

  score1 = EMT76GS(expVal,gseNum)
  score2 = KSScore(expVal,gseNum)
  
  
  return(list(score1,score2))

}


getEMTCor  = function(gse){
  expMat=geneExpGEOdata(gse)
  setwd("../")
  emtScore=getEMTScore(expMat[[1]],gse,expMat[[2]])
  return(emtScore)  
}
