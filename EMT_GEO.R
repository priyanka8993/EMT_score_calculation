source("./GEO_data_preprocess.R")
source("./GEO_Calculate_EMT_Score.R")

args = commandArgs(trailingOnly = TRUE)
print(args)

getEMTCor(args[1])