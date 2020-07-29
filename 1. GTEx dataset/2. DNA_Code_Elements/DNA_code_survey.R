


# read in cnvs and annotation data
cnvs = read.csv("CNVs.csv")
CpG = read.table("CpG_1")
CpG = unique(CpG)
repeatmasker = read.table("repeatmasker_1")
repeatmasker = unique(repeatmasker)
simplerepeats = read.table("simple_repeats_1")
simplerepeats = unique(simplerepeats)

# Define function 
density_function = function(cnv_dataset, annotation_file, column_name, column_number ){
  column = c(rep(0,1035))
  cnv_dataset = cbind(cnv_dataset, column)
  colnames(cnv_dataset)[colnames(cnv_dataset) == "column"] <- column_name
  
  for (i in 1:nrow(cnv_dataset)){
    subgroup = 0
    subgroup = annotation_file[which((annotation_file$V1==cnv_dataset[i,2])&(annotation_file$V2 >= cnv_dataset[i,3] -50000)& (annotation_file$V3<= cnv_dataset[i,3]+50000) ),]
    
    cnv_dataset[i,column_number] = (sum(subgroup$V3)-sum(subgroup$V2))/100000
  }
  
  return(cnv_dataset)
}


# use function
cnvs = density_function(cnvs, simplerepeats, "simple_repeats", 7)
cnvs = density_function(cnvs, repeatmasker, "repeat_masker", 8)
cnvs = density_function(cnvs, CpG, "CpG", 9)


# save to csv file

write.csv(cnvs, "cnvs_DNA_code.csv", row.names = FALSE)

