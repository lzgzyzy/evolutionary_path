###1. read GBM mutation and copy number data 
mut_data<-read.csv(file="/pub5/xiaoyun/Jobs/J22/OriginalData/PanCanAtlas/mc3.v0.2.8.PUBLIC.maf",header = T,sep="\t",stringsAsFactors=F)
mut_table<-mut_data[,c("Tumor_Sample_Barcode","Chromosome","Start_Position" ,"End_Position" ,"Reference_Allele" ,"Tumor_Seq_Allele2","t_alt_count" ,"t_ref_count" ,"Hugo_Symbol","Variant_Classification","HGVSp_Short")]
mut_table[,"Tumor_Sample_Barcode"]<-substr(mut_table[,"Tumor_Sample_Barcode"],1,15)
colnames(mut_table)<-c('Patient' ,'Chr','Start_position' ,'End_position' ,'Reference' ,'Alternate','Variant_freq' ,'Ref_freq' ,'Hugo_Symbol','Variant_Classification','HGVSp_Short')
write.table(mut_table,file="/IJob/J34/clonal_evolution/GBM/data/CCF_estimate/mut_table.txt",sep="\t",row.names=F,quote=F)

mutation.table<-read.csv(file="/IJob/J34/clonal_evolution/GBM/data/CCF_estimate/mut_table.txt",header=T,sep="\t",stringsAsFactors=F)
cna_data<-read.csv(file="/pub5/xiaoyun/Jobs/J22/OriginalData/PanCanAtlas/TCGA_mastercalls.abs_segtabs.fixed.txt",header = T,sep="\t",stringsAsFactors=F)
p2<-read.csv(file="/pub5/xiaoyun/Jobs/J22/OriginalData/PanCanAtlas/TCGA_mastercalls.abs_tables_JSedit.fixed.txt",header = T,sep="\t",stringsAsFactors=F)
ploidy_purity<-p2[,c("ploidy","purity")]
rownames(ploidy_purity)<-p2[,"array"]

######################################################

seg.table2<-cbind.data.frame(cna_data,ploidy_purity[cna_data[,"Sample"],])
setwd("/pub5/xiaoyun/Jobs/J22/EvoClass1.0/Resources/PanCanAtlas/TCGA_GBM")
patient<-dir()
patient_new<-intersect(intersect(mutation.table[,1],seg.table2[,1]),patient)

seg.table_gbm<-seg.table2[seg.table2[,"Sample"]%in%patient_new,]
save(seg.table_gbm,file="/IJob/J34/clonal_evolution/GBM/data/seg.table_gbm_clonal.RData")##
################################################

seg.table<-cna_data[,c("Sample","Chromosome","Start","End","Num_Probes","Modal_Total_CN","Modal_HSCN_1","Modal_HSCN_2")]
seg.table2<-cbind.data.frame(seg.table,ploidy_purity[seg.table[,"Sample"],])
colnames(seg.table2)<-c("SampleID","Chr","Start","End","nProbes","cn","nA","nB","Ploidy","Aberrant Cell Fraction")
rownames(seg.table2)<-NULL
mutation.table_gbm<-mutation.table[mutation.table[,"Patient"]%in%patient_new,]
save(mutation.table_gbm,file="/IJob/J34/clonal_evolution/GBM/data/mutation.table_gbm.RData")
seg.table_gbm<-seg.table2[seg.table2[,"SampleID"]%in%patient_new,]
save(seg.table_gbm,file="/IJob/J34/clonal_evolution/GBM/data/seg.table_gbm.RData")


###2. infer CCF 
TCGA.barcode<-patient
gender<-NA
ANALYSIS.DIR<-"/pub6/temp/zhangyong/clonal_evolution/zy_test/PRAD/"
source("/IJob/J34/clonal_evolution/GBM/data/CCF_estimate/GD.functions.r")
source("/IJob/J34/clonal_evolution/GBM/data/CCF_estimate/main.functions3.r")
source("/IJob/J34/clonal_evolution/GBM/data/CCF_estimate/ccf.wrapper.r")

patient_add<-setdiff(patient_new,mut_file2)
for(i in 10:length(patient_add)){
	print("#############")
	print(i)
	print(patient_add[i])
	clonality.estimation(mutation.table=mutation.table,seg.table=seg.table2,data.type="TCGA_GBM",TCGA.barcode=patient_add[i],sex.chr= FALSE,gender=gender,ANALYSIS.DIR=ANALYSIS.DIR,plotting=FALSE,min.var.prop=0 ,min.alt.reads=0,min.depth= 0)
}




###3. get CCF result
setwd("/pub6/temp/zhangyong/clonal_evolution/zy_test/PRAD/TCGA_GBM")
mut_file<-dir()##
flag<-sapply(mut_file,function(i){
  #i<-mut_file[1]
  setwd(paste0("/pub6/temp/zhangyong/clonal_evolution/zy_test/PRAD/TCGA_GBM/",i))
  length(dir())
})
mut_file2<-mut_file[flag==1]
mut_clonal_data<-NULL
for(i in mut_file2){
  #i<-mut_file2[1]
  temp<-read.csv(file=paste0("/pub6/temp/zhangyong/clonal_evolution/zy_test/PRAD/TCGA_GBM/",i,"/",i,".earlylate.tsv"),sep="\t",stringsAsFactors=F,head=T)
  if(i%in%patient_add){
		if(ncol(temp)==1){
		  print(i)
		  temp<-t(temp)
		  colnames(temp)<-colnames(mut_clonal_data)
		}
  }else{
		if(ncol(temp)==1){
		  print(i)
		  temp<-t(temp)
		  colnames(temp)<-colnames(mut_clonal_data)
		}
		prob.clonal<-apply(temp,1,function(x){sum(as.numeric(x[24:123][96:100]))})
		temp$"prob.clonal"<-prob.clonal
		temp<-temp[,1:23]
  }
  
  mut_clonal_data<-rbind.data.frame(mut_clonal_data,temp) 
}
save(mut_clonal_data,file="/IJob/J34/clonal_evolution/GBM/data/mut_clonal_data_gbm_ccf_p_0.95.RData")

save(mut_clonal_data,file="/IJob/J34/clonal_evolution/GBM/data/mut_clonal_data_gbm_ccf_p_0.85.RData")
save(mut_clonal_data,file="/IJob/J34/clonal_evolution/GBM/data/mut_clonal_data_gbm_ccf_dis.RData")
################################


