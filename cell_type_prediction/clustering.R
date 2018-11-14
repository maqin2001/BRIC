setwd("C:/Users/Juan.XIe/Desktop/test1/") # you may change to your own folder
library(igraph)
library(mclust)
library(MCL)
library(clues)
library(anocva)

RAW <-read.table("Yan_RPKM",header=T,sep="\t")   # ground truth cell type
CellNum <-dim(RAW)[2]-1  # the number of cells 


Graph <-read.csv('Yan_RPKM_graph.csv',header=T,sep=",")
names(Graph) <-c('Node1','Node2','weight')
G <-graph.data.frame(Graph,directed = FALSE)  # convert file into graph
A <- as_adjacency_matrix(G,type="both",attr="weight",names=TRUE,sparse=FALSE)  # convert graph into adjacency matrix
V_name <-rownames(A)   # the vertix
Covered <-length(V_name)  # the #of covered cells
	
MCL_cs <-function(A){
	CLUST <-list()
	for (i in 1:100){
		CLUST[[i]] <-mcl(A,addLoops = FALSE,inflation =i,max.iter=200)
	}
	KK <- as.data.frame(do.call(rbind,lapply(CLUST,'[[',1)))  # extract the number of clusters
	CAN_I <-c(which(as.numeric(as.character(KK$V1))>=2)) 	# results that has more than 5 clusters
	tt <-as.numeric(as.character(KK$V1))
	tt <-sort(table(tt),decreasing=T)[1]
	Final_K <-as.numeric(names(tt))
		
	if (length(CAN_I)!=0){
		MATRIX <-rep(0,Covered)%o%rep(0,Covered)
		for (k in 1:length(CAN_I)){	
			MCL_label <-CLUST[[CAN_I[k]]]$Cluster  # record the label
			ClusterNum <-unique(MCL_label)   # record the number of clusters
			TEMP <-rep(0,Covered)%o%rep(0,Covered)
			temp <-rep(0,Covered) %o% rep(0,length(ClusterNum))
			for (n in 1:length(ClusterNum)){
				index <-which(MCL_label==ClusterNum[n])
				temp[index,n] <-1
				TEMP <-TEMP+temp[,n]%o%temp[,n] 
			}
		MATRIX <-MATRIX+TEMP
		}
		MATRIX <-MATRIX/length(CAN_I)
		rownames(MATRIX) <-colnames(MATRIX) <-rownames(A)
		hc <-hclust(dist(MATRIX))
		memb <-cutree(hc,k=Final_K)
		if (length(rownames(A)) ==CellNum){
			label <-memb
		}else{
			LEFT <-setdiff(names(RAW)[-1],V_name)
			LEFT_Cluster <-rep(Final_K+1,length(LEFT))
			df_cell_label <-data.frame(cell=c(names(memb),LEFT),cluster=c(memb,LEFT_Cluster),K=rep(Final_K+1,CellNum))				
			label <-df_cell_label$cluster
		}	
	}
return(label)
}

SC_cs <-function(A,K){   
	sc <-spectralClustering(A,k=K)
	names(sc) <-rownames(A)
	if (length(rownames(A)) ==CellNum){
		label <-sc
	}else{
		LEFT <-setdiff(names(RAW)[-1],V_name)
		LEFT_Cluster <-rep(K+1,length(LEFT))
		df_cell_label <-data.frame(cell=c(names(sc),LEFT),cluster=c(sc,LEFT_Cluster),K=rep(K+1,CellNum))				
		label <-df_cell_label$cluster
	}
return(label)
}

CLUSTERING <-function(A,K,method){
	if (method=='MCL_cs'){
		label <-MCL_cs(A)
		return(label)
	}else if (method =='SC_cs'){
		label <-SC_cs(A,K)
		return(label)
	}
}

label <-CLUSTERING(A,7,'SC_cs') # get the predicted cell cluster labels

## if you have some reference labels, you can calculate the ARI, RI, FM and JI
target <-read.table('Yan_cell.csv',header=T,sep=',')
head(target)
sorted <-label[match(target$Cell_type,names(label))]  # sort predicted cell cluster labels

ARI <-adjustedRandIndex(sorted,target$Cluster)  # calculate ARI
RI <-adjustedRand(sorted,target$Cluster,randMethod='Rand')
FM <-adjustedRand(sorted,target$Cluster,randMethod='FM')
JI <-adjustedRand(sorted,target$Cluster,randMethod='Jaccard')



