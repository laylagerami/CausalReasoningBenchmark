# Author: Layla Hosseini-Gerami
# This R code snipped runs a modified CausalR ScanR process, but instead of taking the top n nodes from each path length (1-5), will instead take nodes with a p.value <= 0.05 from each path length (consensus nodes). 
# The consensus nodes are then connected to the input nodes via concordant interactions to form subnetworks which can be concatenated to form a consensus subnetwork
# The ranked table at path length 5 is also output

# import packages
library(CausalR)
library(dplyr)

# define functions
rank_hypotheses <- function(ccg,expdata,delta) {
  hypothesis <- RankTheHypotheses(ccg, expdata,delta,correctPredictionsThreshold=1,doParallel=TRUE,numCores=32,writeFile=FALSE)
  hypothesis <- hypothesis[complete.cases(hypothesis), ]
  hypothesis <- data.frame(uniprot = row.names(hypothesis), hypothesis)
  return(hypothesis[!(hypothesis$'p.value'>0.05),])
}

rank_hypotheses_5 <- function(ccg,expdata,delta,results_dir) {
  hypothesis <- RankTheHypotheses(ccg, expdata,delta,correctPredictionsThreshold=1,doParallel=TRUE,numCores=32,writeFile=TRUE,outputDir=results_dir)
  hypothesis <- hypothesis[complete.cases(hypothesis), ]
  hypothesis <- data.frame(uniprot = row.names(hypothesis), hypothesis)
  return(hypothesis[!(hypothesis$'p.value'>0.05),])
}

# load network and input data 
ccg <- CreateCCG("path_to_network.sif")
expData <- ReadExperimentalData("path_to_experimental_data.txt",ccg)

# scan the path lengths 1-5
run1 <- rank_hypotheses(ccg,expData,1)
colnames(run1)[3:9] <- paste("PL1_", colnames(run1[,c(3:9)]), sep = "")
    
run2 <- rank_hypotheses(ccg,expData,2)
colnames(run2)[3:9] <- paste("PL2_", colnames(run2[,c(3:9)]), sep = "")
    
run3 <- rank_hypotheses(ccg,expData,3)
colnames(run3)[3:9] <- paste("PL3_", colnames(run3[,c(3:9)]), sep = "")
    
run4 <- rank_hypotheses(ccg,expData,4)
colnames(run4)[3:9] <- paste("PL4_", colnames(run4[,c(3:9)]), sep = "")
    
run5 <- rank_hypotheses_5(ccg,expData,5,results_dir="path_to_resuts_dir")
colnames(run5)[3:9] <- paste("PL5_", colnames(run5[,c(3:9)]), sep = "")
   
# put all of the results together
big_run <- setNames(data.frame(matrix(ncol = 30, nrow = 0)), c("PL1_Score","PL1_Correct","PL1_Incorrect","PL1_Ambiguous","PL1_p.value","PL1_Enrichment.p.value","PL2_Score","PL2_Correct","PL2_Incorrect","PL2_Ambiguous","PL2_p.value","PL2_Enrichment.p.value","PL3_Score","PL3_Correct","PL3_Incorrect","PL3_Ambiguous","PL3_p.value","PL3_Enrichment.p.value","PL4_Score","PL4_Correct","PL4_Incorrect","PL4_Ambiguous","PL4_p.value","PL4_Enrichment.p.value","PL5_Score","PL5_Correct","PL5_Incorrect","PL5_Ambiguous","PL5_p.value","PL5_Enrichment.p.value"))
all_sig_hyp_2 <- c(row.names(run1),row.names(run2),row.names(run3),row.names(run4),row.names(run5))
all_sig_hyp_2 <- unique(all_sig_hyp_2)
for (hyp in all_sig_hyp_2){
  big_run[hyp,"PL1_Score"]<-run1[hyp,"PL1_Score"]
  big_run[hyp,"PL1_Correct"]<-run1[hyp,"PL1_Correct"]
  big_run[hyp,"PL1_Incorrect"]<-run1[hyp,"PL1_Incorrect"]
  big_run[hyp,"PL1_Ambiguous"]<-run1[hyp,"PL1_Ambiguous"]
  big_run[hyp,"PL1_p.value"]<-run1[hyp,"PL1_p.value"]
  big_run[hyp,"PL1_Enrichment.p.value"]<-run1[hyp,"PL1_Enrichment.p.value"]
  big_run[hyp,"PL2_Score"]<-run2[hyp,"PL2_Score"]
  big_run[hyp,"PL2_Correct"]<-run2[hyp,"PL2_Correct"]
  big_run[hyp,"PL2_Incorrect"]<-run2[hyp,"PL2_Incorrect"]
  big_run[hyp,"PL2_Ambiguous"]<-run2[hyp,"PL2_Ambiguous"]
  big_run[hyp,"PL2_p.value"]<-run2[hyp,"PL2_p.value"]
  big_run[hyp,"PL2_Enrichment.p.value"]<-run2[hyp,"PL2_Enrichment.p.value"]
  big_run[hyp,"PL3_Score"]<-run3[hyp,"PL3_Score"]
  big_run[hyp,"PL3_Correct"]<-run3[hyp,"PL3_Correct"]
  big_run[hyp,"PL3_Incorrect"]<-run3[hyp,"PL3_Incorrect"]
  big_run[hyp,"PL3_Ambiguous"]<-run3[hyp,"PL3_Ambiguous"]
  big_run[hyp,"PL3_p.value"]<-run3[hyp,"PL3_p.value"]
  big_run[hyp,"PL3_Enrichment.p.value"]<-run3[hyp,"PL3_Enrichment.p.value"]
  big_run[hyp,"PL4_Score"]<-run4[hyp,"PL4_Score"]
  big_run[hyp,"PL4_Correct"]<-run4[hyp,"PL4_Correct"]
  big_run[hyp,"PL4_Incorrect"]<-run4[hyp,"PL4_Incorrect"]
  big_run[hyp,"PL4_Ambiguous"]<-run4[hyp,"PL4_Ambiguous"]
  big_run[hyp,"PL4_p.value"]<-run4[hyp,"PL4_p.value"]
  big_run[hyp,"PL4_Enrichment.p.value"]<-run4[hyp,"PL4_Enrichment.p.value"]
  big_run[hyp,"PL5_Score"]<-run5[hyp,"PL5_Score"]
  big_run[hyp,"PL5_Correct"]<-run5[hyp,"PL5_Correct"]
  big_run[hyp,"PL5_Incorrect"]<-run5[hyp,"PL5_Incorrect"]
  big_run[hyp,"PL5_Ambiguous"]<-run5[hyp,"PL5_Ambiguous"]
  big_run[hyp,"PL5_p.value"]<-run5[hyp,"PL5_p.value"]
  big_run[hyp,"PL5_Enrichment.p.value"]<-run5[hyp,"PL5_Enrichment.p.value"]
  
# get consensus nodes
k <- 0
if (!is.na(big_run[hyp,"PL1_Score"])){
  k <- k + 1
}
if (!is.na(big_run[hyp,"PL2_Score"])){
  k <- k + 1
}
if (!is.na(big_run[hyp,"PL3_Score"])){
  k <- k + 1
}
if (!is.na(big_run[hyp,"PL4_Score"])){
  k <- k + 1
}
if (!is.na(big_run[hyp,"PL5_Score"])){
  k <- k + 1
}
big_run[hyp,"Total"] <- k

top_nodes <- big_run[which(big_run$Total == max(big_run$Total)), ]
biggest_pl <- max(big_run$Total)
if(biggest_pl == 0){
  print("no consensus nodes were found")
}

# write the results
list_top_nodes <- rownames(top_nodes)
for(node in list_top_nodes){
  if (grepl('+',node,fixed=TRUE) == TRUE){
    protein <- unlist(strsplit(node,'+',fixed=TRUE))
    sign <- '+1'
  } else { 
    protein <- unlist(strsplit(node,'-',fixed=TRUE))
    sign <- '-1'
  }
  WriteExplainedNodesToSifFile(outputDir="path_to_results_dir",protein,sign,ccg,expData,delta=biggest_pl,correctlyExplainedOnly = TRUE)
}
  


