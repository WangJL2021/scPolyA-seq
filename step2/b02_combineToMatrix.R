# get exp matrix from htseq-count output: Right207

setwd("/home/wangjl/data/apa/190517R/script/")

#1 read cell id: 207
cellid=readLines("/data/jinwf/wangjl/apa/190515/01NO_empty.Right.cellID")
length(cellid) #207 id

#2 get exp matrix
dt=NULL
cell1=cellid[1];#c12_ROW02
tmp=read.table(paste0("../htseq/",cell1,".bam.count"),header = F); #head(tmp)
dt=data.frame(
  enid=tmp$V1
  #paste0(cell1)=tmp$V2
)
dt[,cell1]=tmp$V2
row.names(dt)=dt$enid
head(dt) #only one column
#
for(i in 2:length(cellid)){
  if(i%%10==0){
    print(i)
  }
  cell1=cellid[i];#c12_ROW02
  tmp=read.table(paste0("../htseq/",cell1,".bam.count"),header = F); #head(tmp)
  dt[,cell1]=tmp$V2
}
#rm 1st col
dt=dt[,2:ncol(dt)]
#
dt[1:4,1:5] #########
dim(dt) #[1] 57792   207

#
#rm last 5 rows: htseq log info
dt[(nrow(dt)-4):nrow(dt),1:5] 
# save
write.csv( dt[(nrow(dt)-4):nrow(dt),], file="../last5Row_HTseq.207R.csv" )


#3 get gene names intersected with HGNC
HGNC=readLines("/data/jinwf/wangjl/ref/hg19/HGNC.uniq.37397.genenames")
head(HGNC);length(HGNC) #37397
length(row.names(dt))   #57792 
geneList=intersect( row.names(dt), HGNC );
length(geneList) #37143


#4 the final df 
dim(dt) #[1] 57792   207
dt.2=dt[geneList,]
dim(dt.2) #37143   207
#
dt.2[1:5,1:5]; # head
dt.2[ (nrow(dt.2)-4):nrow(dt.2) ,1:5] #tail

# save
write.csv( dt.2, file="../expression_HTseq.207R.csv" )


# read file, and save as rds 
dd=read.csv("../expression_HTseq.207R.csv",row.names = 1);
dim(dd); dd[1:4,1:4] #37143   207

BC_HeLa.count=dd
save(BC_HeLa.count,file="../BC_HeLa.count.rds")
