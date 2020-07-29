rep_classes = read.table("repeat_classes_test.txt")
rep_classes = unique(rep_classes)



SINE_MIR = rep_classes[which(rep_classes$V2 == "SINE/MIR"),1]
SINE_tRNA = rep_classes[which((rep_classes$V2 == "SINE/tRNA")| (rep_classes$V2 == "SINE?/tRNA")),1]
SINE_Alu = rep_classes[which(rep_classes$V2 == "SINE/Alu"),1]
LINE_Dong = rep_classes[which(rep_classes$V2 == "LINE/Dong-R4"),1]
LINE_L2 = rep_classes[which(rep_classes$V2 == "LINE/L2"),1]
LINE_CR1 = rep_classes[which(rep_classes$V2 == "LINE/CR1"),1]
LINE_RTE = rep_classes[which((rep_classes$V2 == "LINE/RTE-BovB") |(rep_classes$V2 == "LINE/RTE-X") ),1]
LINE_L1 = rep_classes[which((rep_classes$V2 == "LINE/L1")|(rep_classes$V2 == "LINE/L1-Tx1")),1]
LTR_ERVK = rep_classes[which(rep_classes$V2 == "LTR/ERVK"),1]
LTR_ERV1 = rep_classes[which((rep_classes$V2 == "LTR/ERV1") | (rep_classes$V2 == "LTR/ERV1?")),1]
LTR_ERVL = rep_classes[which((rep_classes$V2 == "LTR/ERVL") | (rep_classes$V2 == "LTR/ERVL?")|(rep_classes$V2 == "LTR/ERVL-MaLR")),1]
LTR_Gypsy = rep_classes[which((rep_classes$V2 == "LTR/Gypsy") | (rep_classes$V2 == "LTR/Gypsy?")),1]
RC_Helitron = rep_classes[which((rep_classes$V2 == "RC/Helitron") | (rep_classes$V2 == "RC/Helitron?")),1]
DNA_TcMar = rep_classes[which((rep_classes$V2 == "DNA/TcMar") | (rep_classes$V2 == "DNA/TcMar?")|(rep_classes$V2 == "DNA/TcMar-Mariner") | (rep_classes$V2 == "DNA/TcMar-Pogo") | (rep_classes$V2 == "DNA/TcMar-Tc1")|(rep_classes$V2 == "DNA/TcMar-Tc2")| (rep_classes$V2 == "DNA/TcMar-Tigger")),1]
DNA_PiggyBac = rep_classes[which((rep_classes$V2 == "DNA/PiggyBac") | (rep_classes$V2 == "DNA?/Piggybac?")),1]
DNA_hAT = rep_classes[which((rep_classes$V2 == "DNA/hAT") | (rep_classes$V2 == "DNA/hAT-Ac")|(rep_classes$V2 == "DNA/hAT-Blackjack") | (rep_classes$V2 == "DNA/hAT-Charlie") | (rep_classes$V2 == "DNA/hAT-Tag1")|(rep_classes$V2 == "DNA/hAT-Tip100")| (rep_classes$V2 == "DNA/hAT-Tip100?")|(rep_classes$V2 == "DNA?/hAT-Tip100?")|(rep_classes$V2 == "DNA/hAT?")),1]
SIMPLE_REPEAT = rep_classes[which(rep_classes$V2 == "Simple_repeat/Simple_repeat"),1]
Low_complexity = rep_classes[which(rep_classes$V2 == "Low_complexity/Low_complexity"),1]

repeats = read.table("repeat masker_t")
repeats$V4 = as.character(repeats$V4)
repeats = unique(repeats)


repeat_type_extract = function(repeat_list, repeat_dataset){
  ltr_int_rep = data.frame()
  
  for( i in 1:length(repeat_list)){
    
    ltr_int_rep = rbind(ltr_int_rep, repeat_dataset[repeat_dataset$V4 == repeat_list[i],])
  }
  return(ltr_int_rep)
}

SINE_MIR_reps = repeat_type_extract(SINE_MIR,repeats)
SINE_tRNA_reps = repeat_type_extract(SINE_tRNA,repeats)
SINE_Alu_reps = repeat_type_extract(SINE_Alu,repeats)
LINE_Dong_reps = repeat_type_extract(LINE_Dong,repeats)
LINE_L2_reps = repeat_type_extract(LINE_L2,repeats)
LINE_CR1_reps = repeat_type_extract(LINE_CR1,repeats)
LINE_RTE_reps = repeat_type_extract(LINE_RTE,repeats)
LINE_L1_reps = repeat_type_extract(LINE_L1,repeats)
LTR_ERVK_reps = repeat_type_extract(LTR_ERVK,repeats)
LTR_ER1_reps = repeat_type_extract(LTR_ERV1,repeats)
LTR_ERVL_reps = repeat_type_extract(LTR_ERVL,repeats)
LTR_Gypsy_reps = repeat_type_extract(LTR_Gypsy,repeats)
RC_Helitron_reps = repeat_type_extract(RC_Helitron,repeats)
DNA_TcMar_reps = repeat_type_extract(DNA_TcMar,repeats)
DNA_PiggyBac_reps = repeat_type_extract(DNA_PiggyBac,repeats)
DNA_hAT_reps = repeat_type_extract(DNA_hAT,repeats)
Simple_repeats_reps = repeat_type_extract(SIMPLE_REPEAT,repeats)
Low_complexity_reps = repeat_type_extract(Low_complexity,repeats)

density_function = function(cnv_dataset, annotation_file, column_name, column_number ){
  column = c(rep(0,157))
  cnv_dataset = cbind(cnv_dataset, column)
  colnames(cnv_dataset)[colnames(cnv_dataset) == "column"] <- column_name
  
  for (i in 1:nrow(cnv_dataset)){
    subgroup = 0
    subgroup = annotation_file[which((annotation_file$V1==cnv_dataset[i,1])&(annotation_file$V2 >= cnv_dataset[i,2] -50000)& (annotation_file$V3<= cnv_dataset[i,2]+50000) ),]
    
    cnv_dataset[i,column_number] = (sum(subgroup$V3)-sum(subgroup$V2))/100000
  }
  
  return(cnv_dataset)
}

cnvs = read.csv("CNVs_ASD.csv")



cnvs = density_function(cnvs, SINE_MIR_reps,"SINE/MIR",5)
cnvs = density_function(cnvs, SINE_tRNA_reps,"SINE/tRNA",6)
cnvs = density_function(cnvs, SINE_Alu_reps,"SINE/Alu",7)
cnvs = density_function(cnvs, LINE_Dong_reps,"SINE/Dong",8)
cnvs = density_function(cnvs, LINE_L2_reps,"LINE/L2",9)
cnvs = density_function(cnvs, LINE_CR1_reps,"LINE/CR1",10)
cnvs = density_function(cnvs, LINE_L1_reps,"LINE/L1",11)
cnvs = density_function(cnvs, LTR_ERVK_reps,"LTR/ERVK",12)
cnvs = density_function(cnvs, LTR_ER1_reps,"LTR/ER1",13)
cnvs = density_function(cnvs, LTR_ERVL_reps,"LTR/ERVL",14)
cnvs = density_function(cnvs, LTR_Gypsy_reps,"LTR/Gypsy",15)
cnvs = density_function(cnvs, RC_Helitron_reps,"RC/Helitron",16)
cnvs = density_function(cnvs, DNA_TcMar_reps,"DNA/TcMar",17)
cnvs = density_function(cnvs, DNA_PiggyBac_reps,"DNA/PiggyBac",18)
cnvs = density_function(cnvs, DNA_hAT_reps,"DNA/hAT",19)
cnvs = density_function(cnvs, Simple_repeats_reps,"Simple_repeats",20)
cnvs = density_function(cnvs, Low_complexity_reps,"Low_Complexity",21)



write.csv(cnvs,"cnvs_repeat_class_ASD.csv", row.names = FALSE)



