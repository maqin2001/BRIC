### step2: construct weighted graph  ##

setwd('C:/Users/Juan.XIe/Desktop/test1/BC/')
files <-list.files(path=getwd(),pattern='*blocks',full.names=T, recursive=FALSE)

fOUT <-'C:/Users/Juan.XIe/Desktop/test1/BC/Graph/'  # define the path for the output weighted graph

for (i in 1:length(files)){
	F <-readLines(files[i])  # read .blocks file
	TEMP <-grep('Conds',F,value=T) ## extract condition lines in each BC
	BC <-sapply(strsplit(TEMP,':',2),'[',2) # only keep cell names
	temp1 <-gsub('#','_',BC)  # replace # with _
	temp2 <-gsub('-','_',temp1)  
	temp3 <-gsub('c','C',temp2) 
	BC <-temp3

	CONDS <-as.character()   # store the conditions 
	label_C <-as.numeric()   # store the occurence of one condistions

	for (j in 1:length(BC)){
		BCcond <-unlist(strsplit(BC[j], split = " "))
		BCcond <-BCcond[BCcond!=""]  # exclude the blank string
		CONDS <-c(BCcond,CONDS)
		label_C <-c(label_C,rep(j,length(BCcond)))
	}

	df_C <-data.frame(conds=CONDS,label=label_C)
	uniq_C <-df_C$conds[!duplicated(df_C$conds)]   # unique conditions
	Node <-t(combn(uniq_C,2))
	
	Wt <-rep(-1,dim(Node)[1])
	for (k in 1:dim(Node)[1]){
		member1 <-df_C[which(df_C$conds %in% Node[k,1]),]   # identify which BC the k th Node appear
		member2 <-df_C[which(df_C$conds %in% Node[k,2]),]
		Wt[k] <-length(intersect(member1[,2],member2[,2])) # the weight between two node
	}
	GRAPH <-data.frame(Node[,1],Node[,2],Wt)
	if (dim(GRAPH)[1]!=0)	{
	write.csv(subset(GRAPH,Wt!=0),paste(fOUT,basename(files[i]),".csv",sep=""))
	}
}



