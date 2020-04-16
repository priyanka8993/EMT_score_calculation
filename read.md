# EMT Score calculation 

#This code downloads GEO data for given GSE ID, preprocesses it to generate gene wise expression and  calculates EMT score using metrics mentioned in *Byers et al.,2013* and *Tan et al., 2014*


	## setwd("./EMT_Score_Calculation/")
	## Rscript EMT_GEO.R <GSEID>

#This code generates one folder for each GSEID and generates following output files -
1. <GSEID>_gene-exp.txt - Gene wise expression matrix
2. <GSEID>_EMT_Score_76GS.txt - EMT scores calculated using *Byers et al.,2013* method (76 gene signature)
3. <GSEID>_EMT_KS_Score_GSE.txt - EMT scores calculated using *Tan et al., 2014* method (KS statistic)



Please consider citing *Chakraborty et al.,2020,Front.Bioeng.Biotech.*  article if you find this code useful in your research. 
