# read in CNV dataset and annotation data
sv_total_coords = read.csv("CNVs.csv")
H3K27me3_file = read.table("H3K27me3_ENCFF559PMU.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE)
H3K36me3_file =  read.table("H3K36me3_ENCFF452OKA.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE)
H3K9ac_file = read.table("H3K9ac_ENCFF929GQP.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE)
H3K9me3_file = read.table("H3K9me3_ENCFF260MTZ.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE)
H3K9me2_file = read.table("H3K9me2_ENCFF491BTX.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE)
H3K4me3_file = read.table("H3K4me3_ENCFF192QQV.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE)
H3K27ac_file = read.table("H3K27ac_ENCFF507QTS.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE)
H3K14ac_file = read.table("H3K14ac_ENCFF879VSE.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE)


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



# Check histone modifications surrounding insertion sites

sv_total_coords = density_function(sv_total_coords, H3K27me3_file, "h3k27me3", 7)

sv_total_coords = density_function(sv_total_coords, H3K36me3_file, "h3k36me3", 8)

sv_total_coords = density_function(sv_total_coords, H3K9ac_file, "h3k9ac", 9)

sv_total_coords = density_function(sv_total_coords, H3K9me3_file, "h3k9me3", 10)

sv_total_coords = density_function(sv_total_coords, H3K9me2_file, "h3k9me2", 11)

sv_total_coords = density_function(sv_total_coords, H3K4me3_file, "h3k4me3", 12)

sv_total_coords = density_function(sv_total_coords, H3K27ac_file, "h3k27ac", 13)

sv_total_coords = density_function(sv_total_coords, H3K14ac_file, "h3k14ac", 14)


# save to file
write.csv(sv_total_coords, "cnvs_histone_mod.csv", row.names = FALSE)
