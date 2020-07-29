


# read in data
cnvs = read.csv("ASD_data.csv")
repeatmasker = read.table("repeat masker_t")
repeatmasker = unique(repeatmasker)
chromHmm = read.table("ChromHMM_t")
chromHmm = unique(chromHmm)
H3K27me3_file = read.table("H3K27me3_ENCFF559PMU.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE)
H3K9me3_file = read.table("H3K9me3_ENCFF260MTZ.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE)

# Define function 

density_function = function(cnv_dataset, annotation_file, column_name, column_number ){
  column = c(rep(0,158))
  cnv_dataset = cbind(cnv_dataset, column)
  colnames(cnv_dataset)[colnames(cnv_dataset) == "column"] <- column_name
  
  for (i in 1:nrow(cnv_dataset)){
    subgroup = 0
    subgroup = annotation_file[which((annotation_file$V1==cnv_dataset[i,1])&(annotation_file$V2 >= cnv_dataset[i,2] -50000)& (annotation_file$V3<= cnv_dataset[i,2]+50000) ),]
    
    cnv_dataset[i,column_number] = (sum(subgroup$V3)-sum(subgroup$V2))/100000
  }
  
  return(cnv_dataset)
}

# Checking surroundings
cnvs = density_function(cnvs, H3K27me3_file, "h3k27me3", 5)
cnvs = density_function(cnvs, H3K9me3_file, "h3k9me3", 6)

Weak_Txn = 0    
Weak_Enhancer   = 0
Txn_Elongation = 0
Repressed      = 0
Txn_Transition  = 0


cnvs = cbind(cnvs,Weak_Txn)
cnvs = cbind(cnvs,Weak_Enhancer)
cnvs = cbind(cnvs,Txn_Elongation)
cnvs = cbind(cnvs,Repressed)
cnvs = cbind(cnvs,Txn_Transition)

for (i in 1:nrow(cnvs)){
  counter_wt =0 
  counter_hc = 0
  counter_in = 0
  counter_wp = 0
  counter_ap = 0
  counter_we = 0
  counter_te = 0
  counter_pp= 0
  counter_r = 0
  counter_tt = 0
  counter_se = 0
  counter_rep = 0
  for (j in 1:nrow(chromHmm)){
    if(chromHmm[j,1]==cnvs[i,1] &  chromHmm[j,2] >= cnvs[i,2] -50000 & chromHmm[j,3]<= cnvs[i,2]+50000){
      if(chromHmm[j,4]=="11_Weak_Txn"){
        counter_wt = counter_wt+ ((chromHmm[j,3]-chromHmm[j,2])/100000)
        cnvs[i,7] = counter_wt
      }
      else if(chromHmm[j,4]=="7_Weak_Enhancer" | chromHmm[j,4]=="6_Weak_Enhancer"){
        counter_we = counter_we+ ((chromHmm[j,3]-chromHmm[j,2])/100000)
        cnvs[i,8] = counter_we
      }
      else if(chromHmm[j,4]=="10_Txn_Elongation"){
        counter_te = counter_te+ ((chromHmm[j,3]-chromHmm[j,2])/100000)
        cnvs[i,9] = counter_te
      }
      
      else if(chromHmm[j,4]=="12_Repressed"){
        counter_r = counter_r+ ((chromHmm[j,3]-chromHmm[j,2])/100000)
        cnvs[i,10] = counter_r
      }
      else if(chromHmm[j,4]=="9_Txn_Transition"){
        counter_tt = counter_tt+ ((chromHmm[j,3]-chromHmm[j,2])/100000)
        cnvs[i,11] = counter_tt
      }
      
      
    }
    
  }
}

cnvs = density_function(cnvs, repeatmasker, "repeat_masker", 12)

# save data to file

write.csv(cnv, "ASD_dataset.csv", row.names = FALSE)