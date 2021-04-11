# get top cell id, from gene names.

setwd("~/data/apa/20200701Fig/f1/track/")
getwd()

#1. load data
# (1)
apaM=read.csv("/data/jinwf/wangjl/apa/190705PAS/bed/apa_Matrix_afterFilterByCounts_cell_20222.csv", 
              row.names = 1)
dim(apaM) #20222   225
apaM[1:3,1:4]
#                c01ROW07 c01ROW12 c01ROW24 c01ROW31
#chr10:320403:-         0        0       30      159
#chr10:855485:-         0        0        0        0
#chr10:1085971:-        2       68      101       33

# (2)
apaSite=read.table("/data/jinwf/wangjl/apa/190705PAS/bed/pasPostions_Location_transcriptName-noChrM_noInnerPrime_PY-filterCountCell-motif-Header.bed")
dim(apaSite) # 20222    12
apaSite[1:4,]
#            PASid   chr     pos strand count region transcript                down20   gene
#1  chr10:320403:- chr10  320403      - 10945   UTR3  DIP2C-201 TGAAAATAGCAGTTTCTTAAT  DIP2C
#2  chr10:855485:- chr10  855485      -  3583     PA LARP4B-201 TCAATTTTTCTTTTGATTTTT LARP4B

# (3)
cellInfo=read.table("/data/jinwf/wangjl/apa/20200701Fig/f3/cell_cycle/cellInfo.V6.txt")
dim(cellInfo) #222 8
cellInfo[1:4,]
#              cid cellType0 cluster cellType cellType2 countsPerCell geneNumber cellCycle
#c12ROW03 c12ROW03   unknown       0     BC_0      BC_0       1820666       6659       MG1
#c12ROW04 c12ROW04   unknown       0     BC_0      BC_0       2816130       7930       MG1
table(cellInfo$cellType)
#BC_0        BC_1 HeLa_normal   HeLa_sync 
#92          73          30          27

cid.1=rownames( cellInfo[ which(cellInfo$cellType=="BC_0"), ] )
cid.2=rownames( cellInfo[ which(cellInfo$cellType=="BC_1"), ] )
cid.3=rownames( cellInfo[ which(cellInfo$cellType=="HeLa_normal"), ] )
cid.4=rownames( cellInfo[ which(cellInfo$cellType=="HeLa_sync"), ] )


#2. for a given gene, select cells whose reads > 5
getCids=function(gene="CCND1"){
  #1. get pA site of the gene
  pa.sites=apaSite[which(apaSite$gene==gene), ]
  #2. order by pos, depend on strand
  if(pa.sites$strand[1]=="+"){
    pa.sites=pa.sites[order(pa.sites$pos),]
  }else{
    pa.sites=pa.sites[order(-pa.sites$pos),]
  }
  pa.sites=as.character(pa.sites$PASid)
  #3. get sub from apaMatrix
  apaM.sub=apaM[pa.sites,]
  # 4. filter out cells whose total counts <5
  df=apaM.sub[,apply(apaM.sub, 2,sum)>=5]
  cid=apply(df,2, sum)
  cid=cid[order(-cid)]
  print(dim(df))
  df[, names(cid) ]
}
#
showCounts=function(gene){
  rs=getCids(gene);cids=colnames(rs[,1:20]);
  print("BC:")
  print(rs[, intersect( colnames(rs[,1:20]), c(cid.1, cid.2) ) ]) #BC
  print("HeLa:")
  print(rs[, intersect( colnames(rs[,1:20]), c(cid.3, cid.4) ) ]) #HeLa
  
  return(cids)
}
showCounts("KRAS")

rs=getCids("KRAS");cids=colnames(rs[,1:20]);cids
rs[, intersect( colnames(rs[,1:20]), c(cid.1, cid.2) ) ] #BC
rs[, intersect( colnames(rs[,1:20]), c(cid.3, cid.4) ) ] #HeLa

#
showCounts("IFT20")
showCounts("EMP1")
showCounts("CCDC43") # good no
showCounts("NPM1P39")

########
# get sub bams
cids=showCounts("KRAS"); cids
cids
#[1] "c16ROW10" "c15ROW30" "c10ROW34" "c12ROW16" "c14ROW23" "c19ROW34" "c12ROW10" "c15ROW17"
#[9] "c16ROW38" "c14ROW10" "c13ROW32" "c16ROW32" "c3ROW35"  "c3ROW02"  "c16ROW20" "c12ROW35"
#[17] "c7ROW12"  "c12ROW12" "c15ROW34" "c12ROW37"

# vim KRAS.bed
# chr12   25356111        25380691        KRAS    0       -

# test
#system("~/bin/bedtools intersect -a KRAS.bed -b /data/jinwf/wangjl/apa/190705PAS/hg19/CutA_c16_ROW10_Aligned.sortedByCoord.out.bam  >c16_ROW10_KRAS.bam")
system("~/bin/bedtools intersect -s -a KRAS.bed -b /data/jinwf/wangjl/apa/190705PAS/hg19/CutA_c15_ROW30_Aligned.sortedByCoord.out.bam  >c15_ROW30_KRAS-s.bam")
system("~/bin/bedtools intersect -s -a /data/jinwf/wangjl/apa/190705PAS/hg19/CutA_c15_ROW30_Aligned.sortedByCoord.out.bam -b KRAS.bed  >c15_ROW30_KRAS-s2.bam")

# batch
dirRoot="/data/jinwf/wangjl/apa/190705PAS/hg19/"
for(cid in cids){
  cid=sub("ROW","_ROW",cid)
  print(cid)
  cmd=paste0("~/bin/bedtools intersect -s -a ",dirRoot,"CutA_",cid,"_Aligned.sortedByCoord.out.bam -b KRAS.bed >",cid,"_KRAS.bam" )
  system( cmd )
}

# (2) transfer to Yi
# scp *bam wangjl@y.biomooc.com:/data/wangjl/igv/KRAS/
# $ ls | while read id; do echo $id; samtools index $id;  done;
# to IGV