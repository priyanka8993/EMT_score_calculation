
library(readxl)
library(matlabr)

EMT76GS = function(finalGSEMat,gseID,outDir){

  cat("Calculating EMT Score by using 76 gene signatures ....\n")

  remIdx = which(apply(finalGSEMat,1,function(x) any(x == NaN | x == -Inf)) == TRUE)
  if(length(remIdx) > 0 ) finalGSEMat = finalGSEMat[-remIdx, ]

  finalGSEMat[is.na(finalGSEMat)] = 0

  sampleNum = ncol(finalGSEMat)
  genes = finalGSEMat[,1]
  exp = apply(finalGSEMat[ ,3:sampleNum],2,as.numeric)

  EMTSignature = data.frame(read_excel("../../Gene_signatures/76GS/EMT_signature_76GS.xlsx",col_names = TRUE))
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
  outFile = paste(outDir,paste(gseID, "_EMT_76GS.txt",sep=""),sep = '/')
  if(!file.exists(outFile)) file.create(outFile)
  emtWrite = cbind(colnames(exp),EMTScoreStd)
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
  
  EMTSignature = data.frame(read_excel("../../Gene_signatures/KS/EM_gene_signature_cellLine_KS.xlsx",col_names = FALSE))
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

  outputFile = paste(outDirectory,paste(gseID,"_EMT_KS.txt",sep = ''),sep = '/')
  if(!file.exists(outputFile)) file.create(outputFile)

  ksOut = cbind(colnames(EMTExp),finalScore)
  write.table(ksOut,file = outputFile,sep = '\t',row.names=F,quote= F)
  return(finalScore)

}


emtPredExp = function(geneExpMat,gseID,fileLoc){

  #Predictor + Normalizer gene expression to run MLR model for EMT prediction 
  geneList =  scan(file = "../../Gene_signatures/MLR/genes_for_EMT_score.txt",sep = '\n',what = "vector")
  idx = match(geneList,as.character(geneExpMat[,1]))
  geneExpSubMat = geneExpMat[na.omit(idx), ]
  notFound = setdiff(geneList,as.character(geneExpMat[,1]))
  print(paste("Gene not found in the dataset", notFound))
  sampleNum = ncol(geneExpSubMat)
  outFileName1 = paste(gseID,"_EMT_gene_explevels.txt",sep="")
  writeFileName = paste(fileLoc,outFileName1,sep = '/')
  write.table(geneExpSubMat[,-2],file =writeFileName ,sep = '\t',quote = F,row.names = F)
  #write(writeFileName,file = "../../MLR_Matlab_Code/matlab_input_temp.txt",append=T)
  return(list(geneExpSubMat,writeFileName))
}

getNCIindices = function(EMTexp,seriesID,outPath){
  nci60_data = read.table(file = "../../Gene_signatures/MLR/GPL570-55999.txt",sep = '\t',header = T,stringsAsFactors = F, fill = T, quote = "")
  nci60_indices = match(EMTexp[,1],nci60_data[,11])
  nci60_out = paste(seriesID,"_nci60_use_probe.txt",sep = "")
  nci60_indices_file = paste(outPath,nci60_out,sep = '/')
  write(nci60_indices,file = nci60_indices_file ,sep = '\n')
  #write(nci60_indices_file,file = "../../MLR_Matlab_Code/matlab_input_temp.txt",append=T)
  return(list(nci60_indices,nci60_indices_file))
}



mlrEMTPred = function(expSet,gse_id,gpl_id,in_dir,out_dir){
  
  #cat("EMT predictors expression file generated\n")
  emtPredExpOut = emtPredExp(expSet,gse_id,in_dir)
  emtExpMat = emtPredExpOut[[1]]
  emtWrite = emtPredExpOut[[2]]
  

  NCIindicesOut = getNCIindices(emtExpMat,gse_id,in_dir)
  nci60Indices = NCIindicesOut[[1]]
  nci60Write = NCIindicesOut[[2]]
  
  mlrOutfile  = paste(gse_id,"_emt_score_MLR.txt",sep  = "")
  writeOut = paste(out_dir,mlrOutfile,sep = '/')
  write(rbind(emtWrite, nci60Write,writeOut),file = "../../MLR_Matlab_Code/matlab_input_temp.txt", sep = '\n')

  cat("running MLR EMT predictions (matlab code)...\n")
  run_matlab_script("../../MLR_Matlab_Code/EMT_automated.m")	

  mlrOut = read.table(file  = writeOut,header = T,stringsAsFactors = F)
  file.remove("../../MLR_Matlab_Code/matlab_input_temp.txt")
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
  
  names(corVal) = c("76GS-KS_Cor","76GS-KS_Pval","76GS-MLR_Cor","76GS-MLR_Pval","KS-MLR_Cor","KS-MLR_Pval" )

  outFileName =paste("../","76GS_KS_EMT_Cor_all.txt",sep = '/')
  if(!file.exists(outFileName)) file.create(outFileName)
  write(c(gse_data_id,corVal), file = outFileName,sep = '\t',ncolumns  = 7,append = T)

  return(corVal)
}


getEMTScore = function(expVal,gseNum,gplID,dirIn){

  setwd("../../Output")
  dirOut = paste(getwd(),gseNum,sep = "/")
  dir.create(dirOut)
  setwd(dirOut)
  score1 = EMT76GS(expVal,gseNum,dirOut)
  score2 = KSScore(expVal,gseNum,dirOut)
  score3 = mlrEMTPred(expVal,gseNum,gplID,dirIn,dirOut)
  score_list = list(as.vector(score1), as.vector(score2[,1]), as.vector(score3))
  cor_val = all_scoreCor(score_list,gseNum,dirOut)
  
  return(score_list)
}

remInputFiles=function(curDirectory){
  allFiles = list.files(curDirectory)
  setwd(curDirectory)
  sapply(allFiles,file.remove)
}

getEMTCor  = function(gse){
  expMatInfo=geneExpGEOdata(gse)
  emtScore=getEMTScore(expMatInfo[[1]],gse,expMatInfo[[2]],expMatInfo[[3]])
  #remInputFiles(expMatInfo[[3]])
  return(emtScore)  
}