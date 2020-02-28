
output_table <- function(alleles, allelesCount,sampleName, tableName, path){
  table=NULL
  table=as.data.frame(matrix(data=NA,nrow=1,ncol=2*length(alleles)))
  
  colNamesTable=rep("",2*length(alleles))
  for (i in 1:length(alleles)){
    colNamesTable[i] = paste("allele",i, sep="")
    colNamesTable[i+length(alleles)] = paste("alleleCount",i, sep="")
  }
  colnames(table)=colNamesTable
  rownames(table)=sampleName
  
  for(j in 1:length(alleles)){
    table[1,j]=alleles[j]
    table[1,j+length(alleles)]=allelesCount[j]
  }
  
  fileoutname2=paste(path,tableName,filename,sep="")
  write.table(table,fileoutname2,quote=FALSE,sep="\t")
}

homozygote_check <- function(alleles, allelesCount, finalAlleles, finalCounts){
  if(length(allelesCount) == 0){
    finalAlleles[2]=finalAlleles[1]
    finalCounts[2]=finalCounts[1]
    table2Name="FINAL_table_"
    output_table(alleles=finalAlleles, allelesCount = finalCounts,
                 sampleName = sample,tableName = table2Name, path=path)
    stop("homozygote")
  }
}




#path = "micros from 6348/RPYTHON/Serie Files/"
path = "micros from 1836/RPYTHON/series/"

#setwd(path6348)

#### collect and check arguments provided
#filename="MICROSAT.PCR_gs2_serie.tab"
#filename="MICROSAT.PCR_gs4_serie.tab"
#filename="MICROSAT.PCR_gs8_serie.tab"
#filename="MICROSAT.PCR_gs10_serie.tab"
#filename="MICROSAT.PCR_gs13_serie.tab"
#filename="MICROSAT.PCR_gs16_serie.tab"
#filename="MICROSAT.PCR_glsa22_serie.tab"
filename="MICROSAT.PCR_glsa65_serie.tab"

cutoff=1
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
    if(curCount >= 0.4*highCnt && ((curSeqLength - motifLength) %in% d$seq_length)){
      alleles = append(alleles, sprintf("%i_%i",curSeqLength,curSeries))
      allelesCount = append(allelesCount, curCount)
    }
  }
}


###### prepare table2 to transfer genotype data id + 6 geno + 6 counts
### create and initialize table 2

sample=names(data)[3]
table1Name = "table_"
output_table(alleles=alleles, allelesCount = allelesCount,
             sampleName = sample,tableName = table1Name, path = path)

#############################     Intermediary Output      ###################################

#############################     FINAl OUTPUT     ###################################


finalAlleles=c("","")
finalCounts = c(0,0)
indexOfMax = which.max(allelesCount)

finalAlleles[1] = alleles[indexOfMax]
seriesOfMax = strsplit(finalAlleles[1],"_")[[1]][2]
lengthOfMax = strsplit(finalAlleles[1],"_")[[1]][1]
finalCounts[1] = allelesCount[indexOfMax]
allelesCount = allelesCount[-indexOfMax]
alleles = alleles[-indexOfMax]

homozygote_check(alleles, allelesCount, finalAlleles, finalCounts)


#REMOVE NON DOMINANT SERIES
for (a in length(alleles)){
  curSeries = strsplit(alleles[a],"_")[[1]][2]
  if(curSeries != seriesOfMax){
    allelesCount[a] = NA
    alleles[a] = NA
  }
}

alleles<-alleles[!is.na(alleles)]
allelesCount<-allelesCount[!is.na(allelesCount)]

homozygote_check(alleles, allelesCount, finalAlleles, finalCounts)


indexOfSecond = which.max(allelesCount)
finalAlleles[2] = alleles[indexOfSecond]
finalCounts[2] = allelesCount[indexOfSecond]

table2Name="FINAL_table_"
output_table(alleles=finalAlleles, allelesCount = finalCounts,
             sampleName = sample,tableName = table2Name, path=path)



	
