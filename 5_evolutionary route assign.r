
patient_f<-unique(alt_clonal_data_new2[,"patient"])
clonal_event<-sapply(patient_f,function(x){
    #x<-patient_f[1]
	pos<-which(alt_clonal_data_new2[,"patient"]==x)
	temp<-alt_clonal_data_new2[pos,]
	unique(temp[temp[,"flag"]==2,"Hugo_Symbol"])

})

subclonal_event<-sapply(patient_f,function(x){
    #x<-patient_f[1]
	pos<-which(alt_clonal_data_new2[,"patient"]==x)
	temp<-alt_clonal_data_new2[pos,]
	a<-unique(temp[temp[,"flag"]==1,"Hugo_Symbol"])
	index<-which(x==patient_f)
    #intersect(c(b,a),clonal_event[[index]])
	setdiff(a,clonal_event[[index]])
})

#### 
cand_pair<-sapply(1:length(patient_f),function(i){
	a<-clonal_event[[i]]
	b<-subclonal_event[[i]]
	if(length(a)*length(b)==0){
		temp<-NULL
	}else{
		temp<-unlist(lapply(a,function(z){
			paste0(z,">",b)
		}))
	}

	temp
})


cand_pair_all<-names(which(table(unlist(cand_pair))>=5))
length(cand_pair_all)
p_value<-sapply(cand_pair_all,function(x){
	#x<-cand_pair_all[77]
	pos<-grep(x,cand_pair)
	x_rev<-paste0(rev(unlist(strsplit(x,">"))),collapse=">")
	pos2<-grep(x_rev,cand_pair)
	if(length(pos)>length(pos2))
	p<-binom.test(c(length(pos),length(pos2)),p=0.5,alternative="greater")$"p.value"
	else
	p<-1
	
})

order_co_relation<-intersect(names(which(p.adjust(p_value,method="BH")<0.05)),names(which(p_value<0.05)))
#order_co_relation<-names(which(p_value<0.05))
length(order_co_relation)
order_co_relation_pair<-t(sapply(order_co_relation,function(x){unlist(strsplit(x,">"))}))

#############################################################
##################################################
library(pheatmap)
all_event<-unique(alt_clonal_data_new2[,2])
patient_f<-unique(alt_clonal_data_new2[,1])

ccf_m<-rep(0,length(all_event))
alt_ccf_matrix<-sapply(patient_f,function(x){
	#x=patient_f[1]
	pos<-which(alt_clonal_data_new2[,1]==x)
	temp_event<-unique(alt_clonal_data_new2[pos,2])
	temp_ccf<-sapply(temp_event,function(y){
		max(alt_clonal_data_new2[pos[which(alt_clonal_data_new2[pos,2]==y)],"absolute.ccf"])
	})
	index<-match(temp_event,all_event)
	ccf_m[index]<-temp_ccf
	ccf_m
})
rownames(alt_ccf_matrix)<-all_event

flag<-rep(0,length(all_event))
alt_matrix<-sapply(patient_f,function(x){
	pos<-which(alt_clonal_data_new2[,1]==x)
	index<-match(alt_clonal_data_new2[pos,2],all_event)
	flag[index]<-1
	flag
})
rownames(alt_matrix)<-all_event
sort(rowSums(alt_matrix)/ncol(alt_matrix))


order_co_relation_matrix_raw<-sapply(patient_f,function(x){
	#x=patient_f[1]
	pos<-which(alt_clonal_data_new2[,1]==x)
	apply(order_co_relation_pair,1,function(y){
		index1<-which(alt_clonal_data_new2[pos,2]==y[1])
		index2<-which(alt_clonal_data_new2[pos,2]==y[2])
		if(length(index1)*length(index2)>0)
		return(1)
		else
		return(0)
	})
})
table(colSums(order_co_relation_matrix_raw))
table(rowSums(order_co_relation_matrix_raw))
sort(rowSums(order_co_relation_matrix_raw))
length(which(colSums(order_co_relation_matrix_raw)!=0))/ncol(order_co_relation_matrix_raw)


flag<-rep(0,length(all_event))
order_co_relation_matrix<-sapply(patient_f,function(x){
	#x=patient_f[1]
	pos<-which(alt_clonal_data_new2[,1]==x)
	apply(order_co_relation_pair,1,function(y){
		#y=order_co_relation_pair[5,]
		
		a<-alt_clonal_data_new2[pos[which(alt_clonal_data_new2[pos,2]==y[1])],]
		b<-alt_clonal_data_new2[pos[which(alt_clonal_data_new2[pos,2]==y[2])],]
		if((nrow(a)*nrow(b))==0){
			return(0)
		}else{
		###
			flag<-((max(a[,3])+max(b[,3]))>1)&&((max(a[,3])-(max(b[,3])))>=0)
			if(flag){
				return(1)
			}else{
				return(0)
		    }
		}
	})
})
table(colSums(order_co_relation_matrix))
table(rowSums(order_co_relation_matrix))
length(which(colSums(order_co_relation_matrix)!=0))/ncol(order_co_relation_matrix)
sort(rowSums(order_co_relation_matrix))


###
order_co_relation_filter<-order_co_relation[which(rowSums(order_co_relation_matrix)>=19)]###占总样本数的5%以上
pos<-which(colSums(order_co_relation_matrix[order_co_relation_filter,])!=0)
order_co_relation_matrix_filter<-order_co_relation_matrix[order_co_relation_filter,pos]
order_co_relation_filter_pair<-t(sapply(order_co_relation_filter,function(x){unlist(strsplit(x,">"))}))

length(which(colSums(order_co_relation_matrix_filter)!=0))/ncol(order_co_relation_matrix)
sort(rowSums(order_co_relation_matrix_filter))


pdf(file="gene_path_matrix.pdf")
pheatmap(alt_matrix)
pheatmap(alt_ccf_matrix)
pheatmap(order_co_relation_matrix_filter)
dev.off()

save(order_co_relation,file="order_co_relation.RData")
save(alt_matrix,file="alt_matrix.RData")
save(order_co_relation_matrix,file="order_co_relation_matrix.RData")
save(order_co_relation_matrix_filter,file="order_co_relation_matrix_filter.RData")
