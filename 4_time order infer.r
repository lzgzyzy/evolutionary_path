library(lme4)
library(Matrix)
library(BradleyTerry2)

cna_clonal_data_new2_temp<-cna_clonal_data_new2[,c("Sample","region_more","ccf_merge","flag")]
colnames(cna_clonal_data_new2_temp)<-c("patient","Hugo_Symbol","absolute.ccf","flag")
cna_clonal_data_new2_temp<-na.omit(cna_clonal_data_new2_temp)
alt_clonal_data_new2<-rbind.data.frame(mut_clonal_data_new2[,c("patient","Hugo_Symbol","absolute.ccf","flag")],cna_clonal_data_new2_temp)

setwd("/IJob/J34/clonal_evolution/GBM/results/arm0.6(3)")
save(alt_clonal_data_new2,file="alt_clonal_data_new2.RData")


#alt_clonal_data_new2<-alt_clonal_data_new2[alt_clonal_data_new2[,"flag"]==2,]


driver_event<-unique(alt_clonal_data_new2[,"Hugo_Symbol"])
driver_com<-t(combn(driver_event,2))
order_info<-NULL
driver_com_all<-NULL
for(i in 1:nrow(driver_com)){
    #i=1
	gene1<-driver_com[i,1]
	gene2<-driver_com[i,2]
	pos1<-which(gene1==alt_clonal_data_new2[,"Hugo_Symbol"])
	pos2<-which(gene2==alt_clonal_data_new2[,"Hugo_Symbol"])
	
	temp_patient<-intersect(alt_clonal_data_new2[pos1,"patient"],alt_clonal_data_new2[pos2,"patient"])
	if(length(temp_patient)==0){
		a<-c(NA,NA)
	} else{
		ccf1<-sapply(temp_patient,function(i){xx<-which(i==alt_clonal_data_new2[pos1,"patient"]);max(alt_clonal_data_new2[pos1,][xx,"absolute.ccf"])})
		ccf2<-sapply(temp_patient,function(i){xx<-which(i==alt_clonal_data_new2[pos2,"patient"]);max(alt_clonal_data_new2[pos2,][xx,"absolute.ccf"])})
		
		a<-c(length(which(ccf1<ccf2)),length(which(ccf1>ccf2)))
	}
	driver_com_all<-rbind(driver_com_all,cbind(gene1,gene2))
	order_info<-rbind(order_info,a)
}
rownames(order_info)<-NULL
order_info_new<-na.omit(cbind.data.frame(driver_com_all,order_info))
colnames(order_info_new)<-c("gene1","gene2","ccf1","ccf2")
order_info_new[,1]<-as.character(order_info_new[,1])
order_info_new[,2]<-as.character(order_info_new[,2])


pos<-which(apply(order_info_new[,3:4],1,function(i){max(i)})>=3)###
order_info_new<-order_info_new[pos,]
driver_event<-unique(c(order_info_new[,1],order_info_new[,2]))
dim(order_info_new)
length(driver_event)
order_info_new[,1]<-factor(order_info_new[,1],levels=driver_event)###必须转化成factor格式
order_info_new[,2]<-factor(order_info_new[,2],levels=driver_event)
head(order_info_new)


baseballModel <- BTm(cbind(ccf1, ccf2), gene1, gene2, refcat="TP53",data=order_info_new)
index<-which(BTabilities(baseballModel)[,2]<2)###
relative_time<-sort(BTabilities(baseballModel)[index,1])
error<-BTabilities(baseballModel)[index,2][order(BTabilities(baseballModel)[index,1])]
BTabilities(baseballModel)
cbind(relative_time,error)

setwd("/IJob/J34/clonal_evolution/GBM/results/arm0.6(3)")

event_number<-tapply(alt_clonal_data_new2[,1],alt_clonal_data_new2[,2],function(x){length(unique(x))})
event_number2<-event_number[match(names(relative_time),names(event_number))]
pdf(file="Bradley-Terry timeLine.pdf")
plot(relative_time,length(relative_time):1,pch=21,cex=0.6,xlim=c(-1.3,4),ylab="event rank")
std1<-relative_time-error/2
std2<-relative_time+error/2
for(i in length(relative_time):1){
lines(c(std1[i],std2[i]),c(length(relative_time)-i+1,length(relative_time)-i+1))
} ###
text(relative_time+0.8,length(relative_time):1,paste0(names(relative_time)," (n=",event_number2,")"),cex=0.6)
dev.off()
#

