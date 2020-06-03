specific_patient<-colnames(order_co_relation_matrix_filter)
share_event<-unique(c(order_co_relation_filter_pair[,1],order_co_relation_filter_pair[,2])[duplicated(c(order_co_relation_filter_pair[,1],order_co_relation_filter_pair[,2]))])
p_value_mu_list<-sapply(share_event,function(temp){
	#temp<-"TP53"
	order_co_relation_filter2<-order_co_relation_filter[grep(temp,order_co_relation_filter)]
	driver_com_temp<-t(combn(1:length(order_co_relation_filter2),2))
	rownames(driver_com_temp)<-apply(driver_com_temp,1,function(i){ paste0(order_co_relation_filter2[i[1]],":",order_co_relation_filter2[i[2]])})

	p_value_mu<-t(apply(driver_com_temp,1,function(x){
		# x<-driver_com_temp[1,]; c(order_co_relation_filter2[x[1]],order_co_relation_filter2[x[2]])
		
		p1<-intersect(colnames(order_co_relation_matrix_filter)[which(order_co_relation_matrix_filter[order_co_relation_filter2[x[1]],]==1)],specific_patient)
		non_p1<-setdiff(specific_patient,p1)
		p2<-intersect(colnames(order_co_relation_matrix_filter)[which(order_co_relation_matrix_filter[order_co_relation_filter2[x[2]],]==1)],specific_patient)
		non_p2<-setdiff(specific_patient,p2)
		# ###1)
		m<-length(intersect(p1,p2))
		n<-length(intersect(p2,non_p1))
		p<-length(intersect(p1,non_p2))
		q<-length(specific_patient)-m-n-p #0
		matrix(c(m,n,p,q),2)
		if(m>=2){##
			xx1<-round(fisher.test(matrix(c(m,n,p,q),2),alternative = "less")$p.value,3)
		}else if(m%in%c(1)){
			xx1<-0.001
		}else{
			xx1<-0
		}

		xx2<-xx1
		c(xx1,xx2)
	}))
	p_value_mu
})
lapply(p_value_mu_list,function(p_value_mu){
	c(names(p_value_mu[which(p_value_mu[,1]<0.05),1]))
})


#####
clinical_data2<-read.csv(file="/pub5/xiaoyun/Jobs/J22/OriginalData/PanCanAtlas/TCGA-CDR-SupplementalTableS1.txt",header=T,sep="\t",stringsAsFactors=F)
pos<-which(clinical_data2[,"type"]=="GBM")
clinical_data2<-clinical_data2[pos,]
library(splines)
library(survival)


# clinical_data2[as.numeric(clinical_data2[,"OS.time"])>365*2,"OS"]<-0
# clinical_data2[as.numeric(clinical_data2[,"OS.time"])>365*2,"OS.time"]<-365*2
# clinical_data2[as.numeric(clinical_data2[,"DSS.time"])>365*2,"DSS"]<-0
# clinical_data2[as.numeric(clinical_data2[,"DSS.time"])>365*2,"DSS.time"]<-365*2
# clinical_data2[as.numeric(clinical_data2[,"PFI.time"])>365,"PFI"]<-0
# clinical_data2[as.numeric(clinical_data2[,"PFI.time"])>365,"PFI.time"]<-365

xx<-c("OS.time","OS")
xx<-c("DSS.time","DSS")
xx<-c("PFI.time","PFI")
survival_type<-c(list(c("OS.time","OS")),list(c("DSS.time","DSS")),list(c("PFI.time","PFI")))


get_km_cox_result<-function(class_pos){
	sapply(1:3,function(i){
		xx<-survival_type[[i]]
		k1<-cbind(1,clinical_data2[class_pos[[1]],xx])
		colnames(k1)[1]<-"group"
		k2<-cbind(2,clinical_data2[class_pos[[2]],xx])
		colnames(k2)[1]<-"group"


		addicts_all<-rbind(k1,k2)
		pos<-which(nchar(rownames(addicts_all))==4)
		addicts_all<-addicts_all[pos,]
		colnames(addicts_all)<-c('clinic','survt','status')
		addicts_all[,"survt"]<-as.numeric(addicts_all[,"survt"])
		addicts_all[,"status"]<-as.numeric(addicts_all[,"status"])
		addicts_all<-na.omit(addicts_all)
		y_all<- Surv(addicts_all$survt,addicts_all$status)
		kmfit_all<- survfit(y_all~addicts_all$clinic)
		kmfit_all
		km_p_all<-round(1-pchisq(survdiff(Surv(survt,status)~clinic, data=addicts_all)$chisq,1),4)
		km_p_all

		a<-summary(kmfit_all)$"table"[,"median"]
		flag<-sign(a[1]-a[2])##衡量前者相对于后者是好的还是差的预后
		
		
		xx<-c(survival_type[[i]],"age_at_initial_pathologic_diagnosis","gender")
		k1<-cbind(1,clinical_data2[class_pos[[1]],xx])
		colnames(k1)[1]<-"group"
		k2<-cbind(2,clinical_data2[class_pos[[2]],xx])
		colnames(k2)[1]<-"group"


		addicts_all<-rbind(k1,k2)
		pos<-which(nchar(rownames(addicts_all))==4)
		addicts_all<-addicts_all[pos,]
		colnames(addicts_all)<-c('clinic','survt','status','age','gender')
		addicts_all[,"survt"]<-as.numeric(addicts_all[,"survt"])
		addicts_all[,"status"]<-as.numeric(addicts_all[,"status"])
		addicts_all[,"age"]<-as.numeric(addicts_all[,"age"])
		addicts_all<-na.omit(addicts_all)
		y_all<- Surv(addicts_all$survt,addicts_all$status)
		kmfit_all<- survfit(y_all~addicts_all$clinic)
		kmfit_all
		##进行cox回归，对年龄，肿瘤等级和组织类型进行校正
		Coxph_Result<-coxph(y_all~ clinic+age+gender,data=addicts_all)
		Multi_Result<-cbind(summary(Coxph_Result)$conf.int,summary(Coxph_Result)$coefficients)
		result<-c(km_p_all,flag,Multi_Result[1,c(9,1)])
		names(result)<-c("km_p","km_flag","cox_p","cox_hr")
		result
	})
}



###3. 对于互斥组生存曲线图绘制
setwd("/IJob/J34/rawResource/pan_cancer/GBM/focal/key_peak")
pdf(file="km_cox.pdf"))
par(mfrow=c(3,2))
for(k in 1:length(me)){
	print(k)
	XXX<-unlist(strsplit(me[k],":"))
	for(i in 1:3){
		xx<-survival_type[[i]]
		sample1<-names(which(order_co_relation_matrix[XXX[1],]==1))
		sample2<-names(which(order_co_relation_matrix[XXX[2],]==1))
		sample_class<-c(list(setdiff(sample1,sample2)),list(setdiff(sample2,sample1)))
		class_pos<-lapply(sample_class,function(x){
			sam<-clinical_data2[,"bcr_patient_barcode"]
			pos<-match(substr(x,9,12),substr(sam,9,12))
		})
		sapply(class_pos,length)
		get_km_cox_plot(class_pos)
		}
}
dev.off()

