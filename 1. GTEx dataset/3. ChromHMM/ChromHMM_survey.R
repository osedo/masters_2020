
# read in data
cnvs = read.csv("CNVs.csv")
chromHmm = read.table("ChromHMM_1")
chromHmm = unique(chromHmm)


# prepare dataframe
Heterochrom = 0
Insulator = 0    
Weak_Txn = 0    
Weak_Promoter   = 0
Active_Promoter= 0
Weak_Enhancer   = 0
Txn_Elongation = 0
Poised_Promoter = 0
Repressed      = 0
Txn_Transition  = 0
Strong_Enhancer= 0
Repetitive= 0

cnvs = cbind(cnvs,Heterochrom)
cnvs = cbind(cnvs,Insulator)
cnvs = cbind(cnvs,Weak_Txn)
cnvs = cbind(cnvs,Weak_Promoter)
cnvs = cbind(cnvs,Active_Promoter)
cnvs = cbind(cnvs,Weak_Enhancer)
cnvs = cbind(cnvs,Txn_Elongation)
cnvs = cbind(cnvs,Poised_Promoter)
cnvs = cbind(cnvs,Repressed)
cnvs = cbind(cnvs,Txn_Transition)
cnvs = cbind(cnvs,Strong_Enhancer)
cnvs = cbind(cnvs,Repetitive)

# Survey around insertions
for (i in 1:3){
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
    if(chromHmm[j,1]==cnvs[i,2] &  chromHmm[j,2] >= cnvs[i,3] -50000 & chromHmm[j,3]<= cnvs[i,3]+50000){
      if(chromHmm[j,4]=="11_Weak_Txn"){
        counter_wt = counter_wt+ ((chromHmm[j,3]-chromHmm[j,2])/100000)
        cnvs[i,9] = counter_wt
      }
      else if(chromHmm[j,4]=="13_Heterochrom/lo"){
        counter_hc = counter_hc+ ((chromHmm[j,3]-chromHmm[j,2])/100000)
        cnvs[i,7] = counter_hc
      }
      else if(chromHmm[j,4]=="8_Insulator"){
        counter_in = counter_in+ ((chromHmm[j,3]-chromHmm[j,2])/100000)
        cnvs[i,8] = counter_in
      }
      else if(chromHmm[j,4]=="2_Weak_Promoter"){
        counter_wp = counter_wp+ ((chromHmm[j,3]-chromHmm[j,2])/100000)
        cnvs[i,10] = counter_wp
      }
      else if(chromHmm[j,4]=="1_Active_Promoter"){
        counter_ap = counter_ap+ ((chromHmm[j,3]-chromHmm[j,2])/100000)
        cnvs[i,11] = counter_ap
      }
      else if(chromHmm[j,4]=="7_Weak_Enhancer" | chromHmm[j,4]=="6_Weak_Enhancer"){
        counter_we = counter_we+ ((chromHmm[j,3]-chromHmm[j,2])/100000)
        cnvs[i,12] = counter_we
      }
      else if(chromHmm[j,4]=="10_Txn_Elongation"){
        counter_te = counter_te+ ((chromHmm[j,3]-chromHmm[j,2])/100000)
        cnvs[i,13] = counter_te
      }
      else if(chromHmm[j,4]=="3_Poised_Promoter"){
        counter_pp = counter_pp+ ((chromHmm[j,3]-chromHmm[j,2])/100000)
        cnvs[i,14] = counter_pp
      }
      else if(chromHmm[j,4]=="12_Repressed"){
        counter_r = counter_r+ ((chromHmm[j,3]-chromHmm[j,2])/100000)
        cnvs[i,15] = counter_r
      }
      else if(chromHmm[j,4]=="9_Txn_Transition"){
        counter_tt = counter_tt+ ((chromHmm[j,3]-chromHmm[j,2])/100000)
        cnvs[i,16] = counter_tt
      }
      else if(chromHmm[j,4]=="5_Strong_Enhancer" | chromHmm[j,4]=="4_Strong_Enhancer" ){
        counter_se = counter_se+ ((chromHmm[j,3]-chromHmm[j,2])/100000)
        cnvs[i,17] = counter_se
      }
      else if(chromHmm[j,4]=="14_Repetitive/CNV" | chromHmm[j,4]=="15_Repetitive/CNV" ){
        counter_rep = counter_rep+ ((chromHmm[j,3]-chromHmm[j,2])/100000)
        cnvs[i,18] = counter_rep
      }
    }
  }
}



# Save to file
write.csv(cnvs, "cnvs_ChromHMM.csv", row.names = FALSE)

