

#sample = "HSU6348"
#path = "micros from 6348/RPYTHON/Serie Files/"

sample= "HSU1836"
path = "micros from 1836/RPYTHON/series/"

#setwd(path6348)


fileList <-list.files(path,pattern = "FINAL_",full.names = TRUE)

table=NULL
table=as.data.frame(matrix(data=NA,nrow=length(fileList),ncol=3))

colNamesTable=c("Microsat","allele1","allele2")
rowNamesTable= sapply(strsplit(basename(fileList),"_"), `[`)[4,]
colnames(table)=colNamesTable
rownames(table)=rowNamesTable

alleles=c("","")
for(f in 1:length(fileList)){
  data=read.table(paste("",fileList[f], sep=""), h=T)
  table[f,1] = levels(data[1,1])
  table[f,2] = levels(data[1,2])
}
table[is.na(table)]=""
tableName = paste(sample,"_FINAL_ALLELES.csv",sep="")

#write.csv(table,tableName,row.names = 1)
write.table(table,tableName,quote=FALSE,sep="\t")


