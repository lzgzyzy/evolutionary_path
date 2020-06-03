pos<-which(!is.na(cna_clonal_data[,"Length"]))
cna_clonal_data<-cna_clonal_data[pos,]
cn<-paste0(cna_clonal_data[,"Modal_HSCN_1"],"_",cna_clonal_data[,"Modal_HSCN_2"])
cna_clonal_data2<-cbind.data.frame(cna_clonal_data,cn)
cna_clonal_data2<-cna_clonal_data2[cna_clonal_data2["Chromosome"]!=23,]


# arm_loc1<-read.table(file="/IJob/J34/rawResource/common/ChromoSome_band_hg37.txt",sep="\t",stringsAsFactors=F,head=F)
# arm_loc2<-arm_loc1[arm_loc1[,1]%in%paste0("chr",1:22),]
# arm_loc3<-t(sapply(1:22,function(x){
	# pos1<-intersect(which(arm_loc2[,1]==paste0("chr",x)),grep("p",arm_loc2[,4]))
	# pos2<-intersect(which(arm_loc2[,1]==paste0("chr",x)),grep("q",arm_loc2[,4]))
	# c(max(arm_loc2[pos1,3]),min(arm_loc2[pos2,2]),max(arm_loc2[pos2,3]))
	
# }))
# colnames(arm_loc3)<-c("cen_start","cen_end","chr_end")
# rownames(arm_loc3)<-1:22

# p_arm<-cbind.data.frame(1:22,0,arm_loc3[,1])
# rownames(p_arm)<-paste0(rownames(p_arm),"p")
# q_arm<-cbind.data.frame(1:22,arm_loc3[,2],arm_loc3[,3])
# rownames(q_arm)<-paste0(rownames(q_arm),"q")
# colnames(p_arm)<-colnames(q_arm)<-c("Chromosome","Start","End")
# arm_info<-rbind.data.frame(p_arm,q_arm)

# cnv_info<-arm_info
# cnv_info[,2]<-as.integer(cnv_info[,2])
#save(cnv_info,file="/IJob/J34/rawResource/common/chr_arm_info.RData")

load("/IJob/J34/rawResource/common/chr_arm_info.RData")
#########2
cna_clonal_data3<-cna_clonal_data2
chazhi<-cna_clonal_data3[,"Modal_Total_CN"]-cna_clonal_data3[,"ploidy"]
index1<-which(chazhi<=-0.6)
index2<-which(chazhi>=0.6)
cna_flag<-c(rep("loss",length(index1)),rep("gain",length(index2)))
cna_clonal_data4<-cbind.data.frame(cna_flag,cna_clonal_data3[c(index1,index2),])
cna_clonal_data4<-cna_clonal_data4[cna_clonal_data4[,"Sample"]%in%patient_filter,]


source("/IJob/J34/software/rangeOverlap.R")
strand<-"*"
foi1<-cbind.data.frame(cna_clonal_data4[,c("Chromosome","Start","End")],strand)
gf1<-cbind.data.frame(cnv_info,strand,rownames(cnv_info))
colnames(gf1)[4]<-"region"
results1<-RegionOverlapping.gr(foi1,gf1)
dim(results1)

#save(results1,file="/IJob/J34/rawResource/pan_cancer2/GBM/version5/results1_arm_cna_overlap.RData")
cna_clonal_data_new<-cbind.data.frame(results1[,c("X","OLpercS")],cna_clonal_data4[results1$"Qindex",])
region<-paste0(cna_clonal_data_new[,"X"],"_",cna_clonal_data_new[,"cna_flag"])
cna_clonal_data_new<-cbind.data.frame(region,cna_clonal_data_new)

##
ccf_merge<-rep(0,nrow(cna_clonal_data_new))
ccf_merge_ci95_high<-rep(0,nrow(cna_clonal_data_new))
pos1<-which(cna_clonal_data_new[,"cna_flag"]=="loss")
ccf_merge[pos1]<-cna_clonal_data_new[pos1,"Cancer_cell_frac_a1"]
ccf_merge_ci95_high[pos1]<-cna_clonal_data_new[pos1,"Ccf_ci95_high_a1"]
pos2<-which(cna_clonal_data_new[,"cna_flag"]=="gain")
ccf_merge[pos2]<-cna_clonal_data_new[pos2,"Cancer_cell_frac_a2"]
ccf_merge_ci95_high[pos2]<-cna_clonal_data_new[pos2,"Ccf_ci95_high_a2"]
cna_clonal_data_new<-cbind.data.frame(cna_clonal_data_new,ccf_merge,ccf_merge_ci95_high,stringsAsFactors=F)

###
a1<-cna_clonal_data_new[pos1,"Ccf_ci95_high_a1"]
a2<-cna_clonal_data_new[pos1,"Ccf_ci95_high_a2"]
b1<-cna_clonal_data_new[pos2,"Ccf_ci95_high_a1"]
b2<-cna_clonal_data_new[pos2,"Ccf_ci95_high_a2"]
library(vioplot)
setwd("/IJob/J34/rawResource/pan_cancer2/GBM/version6")
pdf(file="allelic specific CCF2.pdf")
par(mfrow=c(2,2))
plot(cbind(a1,a2),col = densCols(cbind(a1,a2)),xlab="CCF of minor allele",ylab="CCF of major allele",main="copy number loss segments",pch=20)
plot(cbind(b1,b2),col = densCols(cbind(b1,b2)),xlab="CCF of minor allele",ylab="CCF of major allele",main="copy number gain segments",pch=20)
###
m<-length(which(a1>=a2))
n<-length(which(a1<a2))
m2<-length(which(b1>=b2))
n2<-length(which(b1<b2))
chisq.test(matrix(c(m,m2,n,n2),2,2),correct=F)$"p.value"

m<-length(which(a1<=a2))
n<-length(which(a1>a2))
m2<-length(which(b1<=b2))
n2<-length(which(b1>b2))
chisq.test(matrix(c(m,m2,n,n2),2,2),correct=F)$"p.value"
loss<-a2-a1
gain<-b2-b1
vioplot(loss,gain,ylab="CCF of major allele-CCF of minor allele",col=c("blue","red"),names=c("loss","gain"))
dev.off()


# ccf_merge<-apply(cna_clonal_data_new,1,function(x){max(as.numeric(x[c("Cancer_cell_frac_a1","Cancer_cell_frac_a2")]))})
# ccf_merge_ci95_high<-apply(cna_clonal_data_new,1,function(x){
# y<-as.numeric(x[c("Cancer_cell_frac_a1","Cancer_cell_frac_a2")])
# pos<-match(max(y),y) 
# x[c("Ccf_ci95_high_a1","Ccf_ci95_high_a2")][pos]
# })
# cna_clonal_data_new<-cbind.data.frame(cna_clonal_data_new,ccf_merge,ccf_merge_ci95_high,stringsAsFactors=F)




cna_clonal_data_new[which(is.na(cna_clonal_data_new[,"ccf_merge"])),"ccf_merge"]<-0


flag<-rep(1,nrow(cna_clonal_data_new))
#flag[which(cna_clonal_data_new[,"ccf_merge"]>=0.85)]<-2

flag[which(cna_clonal_data_new[,"ccf_merge_ci95_high"]>=1)]<-2
cna_clonal_data_new<-cbind.data.frame(cna_clonal_data_new,flag)
table(flag)
#save(cna_clonal_data_new,file="/IJob/J34/rawResource/pan_cancer2/GBM/version5/cna_clonal_data_new.RData")###

thr<-50
cna_clonal_data_new2<-NULL
event<-unique(as.character(cna_clonal_data_new[,"region"]))
for(i in event){
	#i <-event[2]
	for(j in patient_filter){
		#j <-"TCGA-RR-A6KA-01"  # "TCGA-14-0871-01"
			pos<-intersect(which(cna_clonal_data_new[,"region"]==i),which(cna_clonal_data_new[,"Sample"]==j))
			# print(j)
			# print(length(pos))
			if(length(pos)==0){
				temp<-NULL
			}else if(length(pos)==1){
				if(cna_clonal_data_new[pos,"OLpercS"]<thr) ###
					temp<-NULL
				else
					temp<-cna_clonal_data_new[pos,]
			}else{ 
				cna_clonal_data_new_tmp<-cna_clonal_data_new[pos,]
				a1<-sum(cna_clonal_data_new_tmp[cna_clonal_data_new_tmp[,"flag"]==2,"OLpercS"])
				a2<-sum(cna_clonal_data_new_tmp[cna_clonal_data_new_tmp[,"flag"]==1,"OLpercS"])
				if(a1>=thr){
					index1<-which(cna_clonal_data_new_tmp[,"flag"]==2)
					temp<-cna_clonal_data_new_tmp[index1[1],]
					temp[1,"OLpercS"]<-a1
					temp[1,"Start"]<-min(cna_clonal_data_new_tmp[index1,"Start"])
					temp[1,"End"]<-max(cna_clonal_data_new_tmp[index1,"End"])
					temp[1,"Length"]<-temp[1,"End"]-temp[1,"Start"]
					temp[1,"ccf_merge"]<-median(cna_clonal_data_new_tmp[index1,"ccf_merge"])
		
				} else if(a2>=thr){
					index2<-which(cna_clonal_data_new_tmp[,"flag"]==1)
					temp<-cna_clonal_data_new_tmp[index2[1],]
					temp[1,"OLpercS"]<-a2
					temp[1,"Start"]<-min(cna_clonal_data_new_tmp[index2,"Start"])
					temp[1,"End"]<-max(cna_clonal_data_new_tmp[index2,"End"])
					temp[1,"Length"]<-temp[1,"End"]-temp[1,"Start"]
					temp[1,"ccf_merge"]<-median(cna_clonal_data_new_tmp[index2,"ccf_merge"])
		
				} else{
					temp<-NULL
				}
					
			}
			cna_clonal_data_new2<-rbind.data.frame(cna_clonal_data_new2,temp)
	}
}
dim(cna_clonal_data_new2)

##
cna_clonal_data_new22<-cna_clonal_data_new2[which(cna_clonal_data_new2[,"ccf_merge"]>=0.1),]
dim(cna_clonal_data_new22)

###
cna_clonal_data_new22[,"region"]<-as.character(cna_clonal_data_new22[,"region"])
freq<-sapply(unique(cna_clonal_data_new22[,"region"]),function(x){length(unique(cna_clonal_data_new22[cna_clonal_data_new22[,"region"]%in%x,"Sample"]))})/length(patient_filter) 
cna_event<-names(which(freq>=0.1))
cna_clonal_data_new3<-cna_clonal_data_new22[cna_clonal_data_new22[,"region"]%in%cna_event,]


table(cna_clonal_data_new3[,"flag"])/nrow(cna_clonal_data_new3) 
dim(cna_clonal_data_new3)

table(cna_clonal_data_new3[,"region"])
tapply(cna_clonal_data_new3[,"flag"],cna_clonal_data_new3[,"region"],table)

# 
# cna_clonal_data_new3[cna_clonal_data_new3[,"region"]%in%c("10p_loss","10q_loss"),"region"]<-"10_loss"
# cna_clonal_data_new3[cna_clonal_data_new3[,"region"]%in%c("7p_gain","7q_gain"),"region"]<-"7_gain"
# cna_clonal_data_new3[cna_clonal_data_new3[,"region"]%in%c("19p_gain","19q_gain"),"region"]<-"19_gain"
# cna_clonal_data_new3[cna_clonal_data_new3[,"region"]%in%c("20p_gain","20q_gain"),"region"]<-"20_gain"


###
cna_sample<-unique(cna_clonal_data_new3[,"Sample"])
four_events<-c(list(c("7p_gain","7q_gain")),list(c("19p_gain","19q_gain")),list(c("20p_gain","20q_gain")),list(c("10p_loss","10q_loss")))
cna_clonal_data_new33<-NULL
for(i in cna_sample){
	#i<-cna_sample[1]
	temp_data<-cna_clonal_data_new3[cna_clonal_data_new3[,"Sample"]==i,]
	for(j in four_events){
	#j<-four_events[[1]]
	pos<-match(j,temp_data[,"region"])
	if(length(na.omit(pos))<2){
		xx<-NULL
	}else{
		xx<-temp_data[pos[1],]
		xx["region"]<-gsub("p_","_",j[1])
		xx["ccf_merge"]<-max(temp_data[pos,"ccf_merge"])
		xx["ccf_merge_ci95_high"]<-max(temp_data[pos,"ccf_merge_ci95_high"])
		xx["flag"]<-max(temp_data[pos,"flag"])
	}
	cna_clonal_data_new33<-rbind.data.frame(cna_clonal_data_new33,xx,stringsAsFactors=F)
	}
}
cna_clonal_data_new333<-cna_clonal_data_new3[-which(cna_clonal_data_new3[,"region"]%in%unlist(four_events)),]
cna_clonal_data_new3<-rbind.data.frame(cna_clonal_data_new333,cna_clonal_data_new33,stringsAsFactors=F)
table(cna_clonal_data_new3[,"region"])
tapply(cna_clonal_data_new3[,"flag"],cna_clonal_data_new3[,"region"],table)

patient_f<-unique(c(mut_clonal_data_new2[,"patient"],cna_clonal_data_new3[,"Sample"]))
length(patient_f) ##

####################################################################
cna_clonal_data_new3_temp<-cna_clonal_data_new3[,c("Sample","region","ccf_merge","flag")]
colnames(cna_clonal_data_new3_temp)<-c("patient","Hugo_Symbol","absolute.ccf","flag")
cna_clonal_data_new3_temp<-na.omit(cna_clonal_data_new3_temp)
alt_clonal_data_new2<-rbind.data.frame(mut_clonal_data_new2[,c("patient","Hugo_Symbol","absolute.ccf","flag")],cna_clonal_data_new3_temp)
table(alt_clonal_data_new2[,"Hugo_Symbol"])
table(alt_clonal_data_new2[,"flag"])/nrow(alt_clonal_data_new2) #
dim(alt_clonal_data_new2)
tapply(alt_clonal_data_new2[,"flag"],alt_clonal_data_new2[,"Hugo_Symbol"],table)

save(mut_clonal_data_new2,file=paste0("/IJob/J34/rawResource/pan_cancer2/",cancer,"/version6/mut_clonal_data_new2.RData"))
save(cna_clonal_data_new3,file=paste0("/IJob/J34/rawResource/pan_cancer2/",cancer,"/version6/cna_clonal_data_new3.RData"))
save(alt_clonal_data_new2,file=paste0("/IJob/J34/rawResource/pan_cancer2/",cancer,"/version6/alt_clonal_data_new2.RData"))
