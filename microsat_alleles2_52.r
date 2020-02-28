

#path = "micros from 1836/RPYTHON/series/"
path = "micros from 6348/RPYTHON/Serie Files/"


filename="MICROSAT.PCR_glsa52_serie.tab"
#filename="MICROSAT.PCR_glsa65_serie.tab"
#filename = file.choose()
cutoff=10
motifLength = 2

### read csv data file
data=read.table(paste(path,filename, sep=""), h=T)

######### filtering
data2=data[data$count>=cutoff,]
data2=data2[order(data2$series),]
data2=data2[duplicated(data2$series) | duplicated(data2$series, fromLast = TRUE),]

#Check for existence of observations, if none quit
if (dim(data2)[1] <= 1) {
  stop("Not enough observations. Try running again with smaller cutoff. \n If the ")
}



tab = (table(data2$series))
size = dim(table(data2$series))

dataList <- list()
alleles = c()
allelesCount=c()
for (i in 1:size){
  dataList[[i]] = data2[data2$series==names(tab[i]),]
  dataList[[i]] = subset(dataList[[i]],select = c("count","series","seq_length"))
}



for (d in dataList){
  d=d[order(d$count),]
  curSeries=d$series[1]
  highCnt = max(d$count)
  for(row in 1:(nrow(d))){
    curSeqLength = d[row,]$seq_length
    curCount = d[row,]$count
    if(curCount >= 0.25*highCnt){
      alleles = append(alleles, sprintf("%i_%i",curSeqLength,curSeries))
      allelesCount = append(allelesCount, curCount)
    }
  }
}

dataColNames=colnames(data)

###### prepare table2 to transfer genotype data id + 6 geno + 6 counts
### create and initialize table 2
table=NULL
table=as.data.frame(matrix(data=NA,nrow=1,ncol=12))

rownames(table)=dataColNames[3]
colNamesTable=c("allele1","allele2","allele3","allele4","allele5","allele6","count_allele1","count_allele2","count_allele3","count_allele4","count_allele5","count_allele6")
colnames(table)=colNamesTable

for(i in 1:length(alleles)){
  table[1,i]=alleles[i]
  table[1,i+6]=allelesCount[i]
}


fileoutname2=paste("table_",filename,sep="")
### write table2 to csv file
write.table(table,fileoutname2,quote=FALSE,sep="\t")











	
