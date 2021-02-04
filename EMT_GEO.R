source("./GEO_data_functions_v4.R")
source("./GEO_EMT_Score_func_v2.R")

args = commandArgs(trailingOnly = TRUE)
print(args)

getEMTCor(args[1])