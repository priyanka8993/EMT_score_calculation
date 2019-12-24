
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

emtPredExp = function(geneExpMat,gseID,fileLoc){

  #Predictor + Normalizer gene expression to run MLR model for EMT prediction 
  geneList =  scan(file = "F:/IISC/EMT/EMT_Score/EMT_Score_MLR/Code/R_Code/genes_for_EMT_score.txt",sep = '\n',what = "vector")
  idx = match(geneList,as.character(geneExpMat[,1]))
  geneExpSubMat = geneExpMat[na.omit(idx), ]
  notFound = setdiff(geneList,as.character(geneExpMat[,1]))
  print(paste("Gene not found in the dataset", notFound))
  sampleNum = ncol(geneExpSubMat)
  outFileName1 = paste(gseID,"_EMT_gene_explevels.txt",sep="")
  writeFileName = paste(fileLoc,gseID,outFileName1,sep = '/')
  write.table(geneExpSubMat[,-2],file =writeFileName ,sep = '\t',quote = F,row.names = F)
  write(writeFileName,file = "F:/matlab_input_temp.txt",append=T)
  return(geneExpSubMat)
}


mlrEMTPred = function(expSet,gse_id,gpl_id,in_dir,out_dir){
  emtPredExp(expSet,gse_id,in_dir)
  cat("EMT predictors expression file generated\n")
  cat("running MLR EMT predictions (matlab code)...\n")
  mlrOutfile  = paste(gse_id,"_emt_score_MLR.txt",sep  = "")
  writeOut = paste(out_dir,mlrOutfile,sep = '/')
  write(writeOut,file = "F:/matlab_input_temp.txt",append = T)
  if(gpl_id == "GPL570"){
  	print("MLR prediction for GPL570 data ")
  	run_matlab_script("F:/IISC/EMT/EMT_Score/EMT_Score_MLR/Code/Matlab_Code/EMT_GPL570_automated.m")	
  } else{
  	run_matlab_script("F:/IISC/EMT/EMT_Score/EMT_Score_MLR/Code/Matlab_Code/EMT_nonGPL570_automated.m")	
  }
  mlrOut = read.table(file  = writeOut,header = T,stringsAsFactors = F)
  file.remove("F:/matlab_input_temp.txt")
  return(mlrOut[,"ScoreEMT3_norm"])

}

all_scoreCor = function(scoreList,gse_data_id,out_path){
  
  count = 0
  corVal = NULL
  for(i in 1:2){
    for(j in (i+1):3){
      count = count + 1
      corEst = cor.test(scoreList[[i]],scoreList[[j]])
      corVal = c(corVal,c(corEst$estimate,corEst$p.value))
    }
  }
  
  names(corVal) = c("76GS-KS_Cor","76GS-KS_Pval","76GS-Jason_Cor","76GS-Jason_Pval","KS-Jason_Cor","KS-Jason_Pval" )

  outFileName =paste("./","76GS_KS_EMT_Cor_all.txt",sep = '/')
  if(!file.exists(outFileName)) file.create(outFileName)
  write(c(gse_data_id,corVal), file = outFileName,sep = '\t',ncolumns  = 7,append = T)

  return(corVal)
}


getEMTScore = function(expVal,gseNum,gplID){

  score1 = EMT76GS(expVal,gseNum)
  score2 = KSScore(expVal,gseNum)
  #score3 = mlrEMTPred(expVal,gseNum,gplID,dirIn,dirOut)
  #score_list = list(score1, score2, score3)
  #cor_val = all_scoreCor(score_list,gseNum)
  
  return(0)

}

remInputFiles=function(gsePath,inPath){
  curDirectory = paste(inPath,gsePath,sep  = "/")
  allFiles = list.files(curDirectory)
  sapply(allFiles,file.remove)
}

getEMTCor  = function(gse){
  expMat=geneExpGEOdata(gse)
  setwd("../")
  emtScore=getEMTScore(expMat[[1]],gse,expMat[[2]])
  #remInputFiles(gse,inputPath)
  return(emtScore)  
}