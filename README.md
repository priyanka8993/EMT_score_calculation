# EMT_score_calculation
#### This code downloads microarray GEO data for a given GSE ID, preprocesses it to generate gene wise expression and calculates EMT score using method mentioned in *Byers et al.,2013*, *Tan et al., 2014* and *George et al., 2017*

To run the code successfully you should have all the required R packages installed. 



##### Please run the code in R as follows :
 ###### setwd(“./EMT_score_calculation”) – depending upon where you save it in your desktop change the working directory
 ###### system(paste(“Rscript  ./EMT_GEO.R”, <"GSEID">)) 



If you find this code useful in your research, please cite Chakraborty P, George JT, Tripathi S, Levine H, and Jolly MK. Comparative Study of Transcriptomics-Based Scoring Metrics for the Epithelial-Hybrid-Mesenchymal Spectrum. Front. Bioeng. Biotechnol. 8:220. doi: 10.3389/fbioe.2020.00220 
