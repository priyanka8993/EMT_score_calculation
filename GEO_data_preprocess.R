library(readxl)
library(GEOquery)
library(R.utils)


downloadGEO = function(GSE){
  directoryStr = paste("./",GSE,sep = "")
  dir.create(directoryStr)
  getGEO(GEO = GSE, filename = NULL, destdir = directoryStr,
         GSElimits = NULL, GSEMatrix = TRUE, AnnotGPL = TRUE,
         parseCharacteristics = TRUE,getGPL = TRUE)
  setwd(directoryStr)
  seriesFile = list.files(pattern = "*_series_matrix.txt.gz")
  gunzip(seriesFile[1])
  annotFile = list.files(pattern = "*.annot.gz")
  if(length(annotFile) >= 1){
    gunzip(annotFile[1])
  }
  gzFiles = list.files(pattern = "*.gz")
  file.remove(gzFiles)
  return(directoryStr)
}

getGSE = function(fileDir) {
  name = list.files(pattern = "*_series_matrix.txt$")
  print(name)
  fileHeader = scan(file =name ,sep = '\n',what = "vector", nlines = 500)
  gsestart =  grep("!series_matrix_table_begin",fileHeader)
  lineNum = gsestart + 1
  print(paste("Skipping lines", lineNum,sep = '  '))
  gseMat  = read.table(file = name,sep = '\t',header =T,skip =lineNum,stringsAsFactors = F,fill = T)
  print("loading series matrix of dimension:")
  print(dim(gseMat))
  gseMat = gseMat[complete.cases(gseMat), ]
  return(gseMat)  
}


getAnnotMat = function(dir){
  platFilename = list.files(pattern = "*.annot$")
  if(length(platFilename) == 0){
    submitter_annot = TRUE
    platFilename = list.files(pattern = "*.soft$")
    platName = gsub(".soft","",platFilename)
  } else {
    submitter_annot = FALSE
    platName = gsub(".annot","",platFilename)
  }

  print(submitter_annot)

  fileBegin =  scan(file = platFilename,sep = '\n',what = "vector")
  annotStart = grep("^ID",fileBegin)-1
  print(annotStart)
  annotMat  = read.table(file = platFilename,sep = '\t',header =T,skip =annotStart,quote = "",stringsAsFactors = F,fill = T,check.names = F)
  print("loading probe annotation matrix of dimension:")
  print(dim(annotMat))

  return(list(annotMat,submitter_annot,platName))
}


getGeneCol = function(gplID){
  geneSymbolCol = read.table(file = "./gene_symbol_mapping.txt",header = F, stringsAsFactors = F,sep = '\t')
  gplIdx = which(geneSymbolCol[,1] == gplID)
  geneCol = geneSymbolCol[gplIdx,3] 
  return(geneCol)
}



readSampleNames = function(Dir){
  filename = list.files(pattern = "*_series_matrix.txt")
  fileHeader = scan(file = filename ,sep = '\n',what = "vector", nlines = 500)
  lineNumber =  grep("!Sample_title",fileHeader)
  print(lineNumber)
  sampleNames =  scan(file = filename,sep = '\t',what = "vector", nlines = 1,skip = lineNumber)
  sampleNames = sampleNames[-1]
  sampleNames=gsub("[[:punct:]]","_",gsub("\\s+","_",sampleNames))
  print(paste("Number of samples :", length(sampleNames), sep  = " "))
  return(sampleNames)
}

checkGeneSymbol = function(geneNames){
  slashGene = grep(".//.",geneNames)
  if(length(slashGene) >= 100 ){
    newGeneName = trimws(sapply(strsplit(geneNames,split = "//"),function(x) x[2]))
  } else{
    newGeneName = geneNames
  }
  return(newGeneName)
}


getGeneWiseExp = function(inputDir,gseMat,annotMat,otherAnnot,annotID,sampleName){
  
  commonGenes = intersect(gseMat[,1],annotMat[,1])
  cat(paste("# common probes:",length(commonGenes),"\n",sep = ' '))
  
  ## Match probe order
  probeIdx = match(commonGenes,gseMat[,1])
  gseMatFiltered = as.data.frame(gseMat[probeIdx,])
  probeAnnotIdx = match(commonGenes,annotMat[,1])
  annotMatFiltered = annotMat[probeAnnotIdx, ]
  
  gseMatNum = as.data.frame(apply(gseMatFiltered[ ,-1],2,as.numeric))  
  rownames(gseMatNum) = gseMatFiltered[,1]
  
    ## Substitute null with 0
  gseMatNum[is.null(gseMatNum)]  = 0

  ## Log2 transformation
  if(any(gseMatNum >= 100)) {
    cat("log2 transformation---\n") 
    gseMatNum = log2(gseMatNum + 1)
  }
  
  ## Unique genes 
  if(otherAnnot == TRUE){
    cat("Other platform annotation file\n")
    gene_name_col = getGeneCol(annotID)
    gene_symbol = checkGeneSymbol(annotMatFiltered[ ,gene_name_col])
    splitMat = split(gseMatNum,gene_symbol)
    cat(paste("# of unique genes per probe:", length(splitMat),"\n",sep = " "))
  } else {
    splitMat = split(gseMatNum,annotMatFiltered[,"Gene symbol"]) 
    cat(paste("# of unique genes per probe:", length(splitMat),"\n",sep = " "))
  }
  
  ## Mean expression of all the probes per gene
  cat("Getting mean expression of all probes per gene\n")
  expMat = matrix(0,nrow =length(splitMat),ncol = ncol(gseMatNum))
  probeNames=geneNames = list()
  
  for(i in 1:length(splitMat)){
    if(nrow(splitMat[[i]]) > 1 ){
      meanExp = apply(splitMat[[i]],2,mean)
      #maxIdx = which(rowSum == max(rowSum))
      expMat[i, ] = as.matrix(meanExp)
      probeNames[[i]] = paste(rownames(splitMat[[i]]),collapse= ",")
      geneNames[[i]] = names(splitMat)[i]
    } else {
      expMat[i, ] = as.matrix(splitMat[[i]][1,])
      probeNames[[i]] = paste(rownames(splitMat[[i]]),collapse= ",")
      geneNames[[i]] = names(splitMat)[i]
    }
  }

  ## Substitute NaNs with 0
  expMat[is.nan(expMat)]  = 0
  rowZero1 = which(apply(expMat,1,function(x) sum(x)) == 0)

  if(length(rowZero1) >= 1 ) expMat = expMat[-rowZero1, ]


  expMat[is.na(expMat)]  = 0

  expMatNew = cbind(toupper(unlist(geneNames)),unlist(probeNames),expMat)
  colnames(expMatNew) = c("Gene","Probe_ID",sampleName)
  return(expMatNew)
}

geneExpGEOdata = function(GSEID){
  path = downloadGEO(GSEID)
	#path = paste(inDir,GSEID,sep = '/')
  	seriesMat = getGSE(path)
	gplMatInfo = getAnnotMat(path)
  	sampleTitle = readSampleNames(path)
	geneExpMat = as.data.frame(getGeneWiseExp(path,seriesMat,gplMatInfo[[1]],gplMatInfo[[2]],gplMatInfo[[3]],sampleTitle))
	outFile = paste(GSEID,"_gene-exp.txt",sep = '')
	write.table(geneExpMat,outFile,sep = '\t',quote = F,row.names = F)
	return(list(geneExpMat,gplMatInfo[[3]]))
}

