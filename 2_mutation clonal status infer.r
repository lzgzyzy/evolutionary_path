setwd("/pub5/xiaoyun/Jobs/J22/EvoClass1.0/Resources/PanCanAtlas")
my_file<-dir()[grep("mut.cff.RData",dir())]
setwd("/IJob/J34/rawResource/pan_cancer2")
cancer_all<-dir()
load("/IJob/J34/clonal_evolution/GBM/data/CCF_estimate/seg.table2.RData")
source("/IJob/J34/rawResource/pan_cancer2/km_cox_function2.r")


kk=4
cancer<-cancer_all[kk]
setwd(paste0("/IJob/J34/rawResource/pan_cancer2/",cancer))
temp_file<-my_file[grep(cancer,my_file)]
load(paste0("/pub5/xiaoyun/Jobs/J22/EvoClass1.0/Resources/PanCanAtlas/",temp_file))
patient<-intersect(com.data[,1],seg.table2[,"Sample"])
length(patient)
mut_clonal_data<-com.data##突变CCF数据
cna_clonal_data<-seg.table2[seg.table2[,"Sample"]%in%patient,]###拷贝数CCF数据



dim(mut_clonal_data)

mut_clonal_data<-mut_clonal_data[mut_clonal_data[,"obs.VAF"]>=0.05,]
flag<-rep(1,nrow(mut_clonal_data))
pos<-which(mut_clonal_data[,"absolute.ccf.0.95"]>=1)
#pos<-which(mut_clonal_data[,"prob.clonal"]>=0.5)
flag[pos]<-2
mut_clonal_data_new<-cbind.data.frame(mut_clonal_data,flag)
#save(mut_clonal_data_new,file="/IJob/J34/rawResource/pan_cancer2/GBM/version5/mut_clonal_data_new.RData")###所有突变的克隆性数据，用于异质性分析

patient_filter<-names(which(table(mut_clonal_data_new[,"patient"])<1000))##
mut_clonal_data_new<-mut_clonal_data_new[mut_clonal_data_new[,"patient"]%in%patient_filter,]
mut_clonal_data_new<-mut_clonal_data_new[mut_clonal_data_new[,"absolute.ccf"]>=0.1,]

#mut_clonal_data_new<-mut_clonal_data_new[mut_clonal_data_new[,"Variant_Classification"]!="Silent",]
mut_clonal_data_new<-mut_clonal_data_new[mut_clonal_data_new[,"Variant_Classification"]%in%c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Splice_Site","Translation_Start_Site"),]
dim(mut_clonal_data_new)


mutsig_gene<-read.table(paste0("/IJob/J34/rawResource/pan_cancer2/",cancer,"/data_mutsig.txt"),sep="\t",stringsAsFactors=F,head=T)
driver<-mutsig_gene[1:max(which(as.numeric(mutsig_gene[,"q"])<=0.25)),"gene"]

freq<-sapply(driver,function(x){length(unique(mut_clonal_data_new[mut_clonal_data_new[,"Hugo_Symbol"]%in%x,"patient"]))})/length(patient_filter)
CGC<-read.csv(file="/IJob/J34/rawResource/common/COSMIC.CGC.20180703.tsv",header = T,sep="\t",stringsAsFactors=F)[,"Gene.Symbol"]

driver_f<-intersect(names(which(freq>=0.05)),CGC)
mut_clonal_data_new2<-mut_clonal_data_new[mut_clonal_data_new[,"Hugo_Symbol"]%in%driver_f,]

length(unique(mut_clonal_data_new2[,"patient"]))# 
table(mut_clonal_data_new2[,"flag"])/nrow(mut_clonal_data_new2)# 
sort(table(mut_clonal_data_new2[,"Hugo_Symbol"]))
dim(mut_clonal_data_new2)