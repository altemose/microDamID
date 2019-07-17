#R Commands for Normalization of Single-Cell DamID Data from Altemose et al. 2019
#Nicolas Altemose, UC Berkeley, 2019

#read in bulk plus SC coverage across all 250kb bins in the genome (see "0_Alignment*.README.txt" for generation of this file")
data0 = read.table("coverage_logFC_250kb_ALLscCov_goldflags.sorted.bed")

#create a vector of cell identifiers
cellnames1 = vector()
cellnames2 = vector()
cellnames3 = vector()
for(i in 1:18){
	cell = paste0("NA",sprintf("%03d", i))
	cellnames1[i]=cell
	cellnames2[i]=paste0(cell,"E") #E for expected values
	cellnames3[i]=paste0(cell,"OE") #OE for observed-expected values
}
keep = c(3,4,5,6,8,9,10,11,13,14,15) #vector of indices, after removing Dam-only cells, anomalous NA007, and cells with >5% plasmid DNA

#column names; e.g. SumBulkDL1 is the number of q1 bulk Dam-LMNB1 reads mapping to that bin from replicate 1; SumBulkD1 is the number of q1 bulk Dam-only reads mapping to that bin from replicate 1
names(data0) = c("chr","start","stop","logFC","SumBulkDL1","SumBulkDL2","SumBulkD1","SumBulkD2",cellnames1,"cLADgold","ciLADgold")

#remove outlier bins by bulk DamID coverage
data=subset(data0,(SumBulkDL1+SumBulkDL2)>0 & ((SumBulkDL1+SumBulkDL2)<=quantile((data0$SumBulkDL1+data0$SumBulkDL2),0.99) | (SumBulkD1+SumBulkD2)<=quantile((data0$SumBulkD1+data0$SumBulkD2),0.99)))
#only keep bins where one of the 11 cells has >0 coverage (eliminates centromeres and other unmappable regions)
data = data[rowSums(data[,keep+8])>0,]

#sum total coverage from both bulk Dam-only replicates and compute log
total =  sum(data$SumBulkD1+data$SumBulkD2)
logtot = log(total)

#for each bin in each cell, multiply the total number of unique fragments per cell [sum(data[,8+i])] by the fraction of bulk Dam-only reads found in that bin [exp(log((data$SumBulkD1+data$SumBulkD2))-logtot)]
#nm: done in log space to avoid overflow
 for(i in 1:18){
 	data[,28+i] = sum(data[,8+i])*exp(log((data$SumBulkD1+data$SumBulkD2))-logtot)
 }
#now compute difference between observed and expected values in each bin
 for(i in 1:18){
  	data[,46+i] = data[,8+i]-data[,28+i]
 }
#update column names
  names(data) = c("chr","start","stop","logFC","SumBulkDL1","SumBulkDL2","SumBulkD1","SumBulkD2",cellnames1,"cLADgold","ciLADgold",cellnames2,cellnames3)

#save results
save(data,file="ObservedExpectedTable.r")
