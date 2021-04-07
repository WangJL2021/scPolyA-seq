#4_2_DEapa—notDEG_HeLa_syncVSnormal.R

setwd('/data/jinwf/wangjl/apa/191111Figure/f4/DEA_not_DEG/')
getwd()

library('Cairo')

#################
#1. read data:
deg=read.csv('/data/jinwf/wangjl/apa/191111Figure/f5/DEG_RNA_all/02-DESeq2_ALL_sync_VS_normal_HeLa(RNA).csv',row.names = 1);
dim(deg) #[1] 8095    6
head(deg)
#        baseMean log2FoldChange     lfcSE      stat       pvalue         padj
#ACTA2   399.35413      10.564362 0.6307483  16.74893 5.765039e-63 4.666799e-59
#S100A7  166.37318      -9.150554 0.6970454 -13.12763 2.287044e-39 9.256810e-36

deApa=read.csv('/data/jinwf/wangjl/apa/191111Figure/f4/DEG_APA/DESeq2_ALL_sync_VS_normal_HeLa(APA).csv', row.names = 1)
dim(deApa) #8766    6
head(deApa)
#                 baseMean log2FoldChange     lfcSE     stat       pvalue         padj
#chr6:153986769:+ 651.8086       9.398069 0.5130505 18.31802 5.943827e-75 5.210359e-71
#chr2:143857505:- 518.0942       8.766839 0.5217277 16.80348 2.301483e-63 1.008740e-59

keyword='DEapa_notDEG_syncVSnormal_HeLa'


#####
# define sig
deg.sig=deg
deg.sig$threshold='ns'
deg.sig[which( deg$log2FoldChange>log2(1.5) & deg$padj<0.05 ),]$threshold='up'
deg.sig[which( deg$log2FoldChange< -log2(1.5) & deg$padj<0.05 ),]$threshold='down'
table(deg.sig$threshold)
#down n.s.   up 
#170 7635  290

#
deApa.sig=deApa
deApa.sig$threshold='ns'
deApa.sig[which( deApa.sig$log2FoldChange>log2(1.5) & deApa.sig$padj<0.05 ),]$threshold='up'
deApa.sig[which( deApa.sig$log2FoldChange< -log2(1.5) & deApa.sig$padj<0.05 ),]$threshold='down'
table(deApa.sig$threshold)
#down   ns   up 
#153 8345  268


#####
# annotate apa with gene symbols
#source('/home/wangjl/pylib/con_mysql.R')
#
apaInfo=read.table("/data/jinwf/wangjl/apa/190705PAS/bed/pasPostions_Location_transcriptName-noChrM_noInnerPrime_PY-filterCountCell-motif-Header.bed",
                   row.names = 1, header = T)
head(apaInfo)
table(apaInfo$gene=="")




# add gene symbol and region
index2=match(  rownames(deApa.sig), apaInfo$PASid )
table(is.na(index2)) #all valid
length(index2) #8766
#
deApa.sig$gene=apaInfo[index2, 'gene']
deApa.sig$region=apaInfo[index2, 'region']

table(deApa.sig$region)

# add mRNA geneLog2FC, shreshold, padj
head(deApa.sig)
index3=match(deApa.sig$gene, row.names(deg.sig) )
head(index3)
table(is.na(index3))
#FALSE  TRUE 
#7068  1698 
deApa.sig$geneLog2FC=deg.sig[index3, ]$log2FoldChange
deApa.sig$geneThreshold=deg.sig[index3, ]$threshold
deApa.sig$genePadj=deg.sig[index3, ]$padj
table( is.na(deApa.sig$geneLog2FC) )
#FALSE  TRUE 
#7068  1698

dim(deApa.sig) #8766   12

#save csv
write.csv(deApa.sig, paste0('01_annotate_deApa_with_symbol_geneFC',keyword,'.csv') )
#


###################
#2. up, down gene
deg.up= deg.sig[which(deg.sig$threshold=='up'),]; dim(deg.up) #290
deg.down=deg.sig[which(deg.sig$threshold=='down'),]; dim(deg.down) #170

# up,down apa
deApa.up= deApa.sig[which(deApa.sig$threshold=='up'),]; dim(deApa.up) #268
deApa.down= deApa.sig[which(deApa.sig$threshold=='down'),]; dim(deApa.down) #153

table(deApa.up$region)
table(deApa.down$region)

#remove na,
deApa.sig2=deApa.sig[!is.na(deApa.sig$geneThreshold),]
head(deApa.sig2)
deApa.sig2$threshold=factor(deApa.sig2$threshold, levels = c('up','ns','down') )
deApa.sig2$geneThreshold=factor(deApa.sig2$geneThreshold, levels = c('up','ns','down') )
dim(deApa.sig2) #7161   11
#
g=ggplot(deApa.sig2, aes( log2FoldChange,  geneLog2FC, color=geneThreshold))+
  geom_point( alpha=0.5,size=0.4)+
  geom_vline(aes(xintercept=log2(1.5)), colour="#aaaaaa", linetype="dashed")+
  geom_vline(aes(xintercept=-log2(1.5)), colour="#aaaaaa", linetype="dashed")+
  #
  geom_hline(aes(yintercept=log2(1.5)), colour="#aaaaaa", linetype="dashed")+
  geom_hline(aes(yintercept=-log2(1.5)), colour="#aaaaaa", linetype="dashed")+
  #
  labs(title="Gene up, but polyA site may down", x= "log2FC of polyA sites", y="log2FC of gene")+
  scale_color_manual('Gene',values=c('red','grey','blue'))+
  theme_test(); g#
#
library(ggrepel)
#
dd_text=deApa.sig2[which(deApa.sig2$threshold!='ns' | deApa.sig2$geneThreshold!='ns'),]; dim(dd_text) #635
head(dd_text)

# reverse direction: 2nd, geneUp while apa down
dd_text.2=dd_text[which( dd_text$log2FoldChange< (-log2(1.5)) & dd_text$geneLog2FC>log2(1.5) ), ]; dim(dd_text.2) #17  11
#
dd_text.2s=dd_text.2[which(dd_text.2$padj<0.05 & dd_text.2$genePadj<0.05),]; dim(dd_text.2s)
dd_text.2=dd_text.2[order(dd_text.2$log2FoldChange),]
head(dd_text.2s)
g2=g+geom_text_repel(data=dd_text.2s, aes(log2FoldChange,  geneLog2FC, label=gene),
                     alpha=0.4,colour="black",size=3)+labs(subtitle = "couter trend: g+*, a-*(n=1)");g2#
#
# reverse direction: 4th, geneDown while apa Up
dd_text.4=dd_text[which( dd_text$log2FoldChange>log2(1.5) & dd_text$geneLog2FC< (-log2(1.5)) ), ]; dim(dd_text.4) #16 11
dd_text.4s=dd_text.4[which(dd_text.4$padj<0.05 & dd_text.4$genePadj<0.05),]; dim(dd_text.4s) #1
head(dd_text.4s)
g4=g+geom_text_repel(data=dd_text.4s, aes(log2FoldChange,  geneLog2FC, label=gene),
                     alpha=0.4,colour="black",size=3)+labs(subtitle = "couter trend: g-*, a+*(n=1)");g4#
#
CairoPDF(file= paste0("02-volcano_counterSig_" , keyword,'.pdf'), width=7, height=3)
grid.arrange(
  g2,g4,
  nrow=1
)
dev.off()




# sig 1st
dd_text.s1=dd_text[which( dd_text$log2FoldChange>log2(1.5) & dd_text$geneLog2FC> log2(1.5) ), ]; dim(dd_text.s1) #265 12
dd_text.s1.1=dd_text.s1[which(  (dd_text.s1$threshold=='ns' & dd_text.s1$geneThreshold=='up')  ),]; dim(dd_text.s1.1) #130
dd_text.s1.2=dd_text.s1[which( (dd_text.s1$threshold=='up' & dd_text.s1$geneThreshold=='ns') ),]; dim(dd_text.s1.2) #59
head(dd_text.s1.2)
#g+* while apa+ not
g1=g+geom_text_repel(data=dd_text.s1.1[115:130,], aes(log2FoldChange,  geneLog2FC, label=gene),
                     alpha=0.4,colour="black",size=3)+labs(subtitle = "couter significant: g+*, a+(n=130)");g1#
#g+ not while apa+*
g1.2=g+geom_text_repel(data=dd_text.s1.2[1:15,], aes(log2FoldChange,  geneLog2FC, label=gene),
                     alpha=0.4,colour="black",size=3)+labs(subtitle = "couter significant: g+, a+*(n=59)");g1.2#
# sig, 3rd
dd_text.s3=dd_text[which( dd_text$log2FoldChange< (-log2(1.5)) & dd_text$geneLog2FC< (-log2(1.5)) ), ]; dim(dd_text.s3) #163 12
dd_text.s3.1=dd_text.s3[which(  (dd_text.s3$threshold=='ns' & dd_text.s3$geneThreshold=='down')  ),]; dim(dd_text.s3.1) #77
dd_text.s3.2=dd_text.s3[which( (dd_text.s3$threshold=='down' & dd_text.s3$geneThreshold=='ns') ),]; dim(dd_text.s3.2) #27
head(dd_text.s3.2)
#g-* while apa- not
g3=g+geom_text_repel(data=dd_text.s3.1[62:77,], aes(log2FoldChange,  geneLog2FC, label=gene),
                     alpha=0.4,colour="black",size=3)+labs(subtitle = "couter significant: g-*, a-(n=77)");g3 #
#g+ not while apa+*
g3.2=g+geom_text_repel(data=dd_text.s3.2[1:15,], aes(log2FoldChange,  geneLog2FC, label=gene),
                       alpha=0.4,colour="black",size=3)+labs(subtitle = "couter significant: g-, a-*(n=27)");g3.2#
#

CairoPDF(file= paste0("02-volcano_sigOrNot_" , keyword,'.pdf'), width=7, height=6)
grid.arrange(
  g1, g1.2,
  g3,g3.2,
  nrow=2
)
dev.off()
#




###################
#3. check

# APA matrix
APA.counts=read.table('/data/jinwf/wangjl/apa/190705PAS/bed/freq2/all225_matrix_APACountsV4.txt',
                      header = T,row.names = 1)
dim(APA.counts) #267961   225
colnames(APA.counts)= sub("_","", colnames(APA.counts) )
APA.counts[1:8,1:10]

#
APA.counts.norm=apply(APA.counts, 2, function(x){
  log2( 1e6*x/sum(x) +1 )
} )
dim(APA.counts.norm) #
APA.counts.norm[1:10,1:10]


#
# read matrix
rna.counts=read.csv('/home/wangjl/data/apa/190530Mix/BC_HeLa.225cells.count.V3.csv',header = T,row.names = 1)
dim(rna.counts) #18679   225
rna.counts[1:8,1:10]

rna.counts.norm=apply(rna.counts, 2, function(x){
  log2( 1e6*x/sum(x) +1 )
} )
dim(rna.counts.norm) #
rna.counts.norm[1:10,1:10]

# cell info
cellInfo=read.table("/data/jinwf/wangjl/apa/20200701Fig/f3/cell_cycle/cellInfo.V6.txt", 
                    row.names = 1, stringsAsFactors = F)
head(cellInfo)

#
# 2 groups: sync and normal
cid.norm=cellInfo[which(cellInfo$cellType=="HeLa_normal"),"cid"] #29
cid.sync=cellInfo[which(cellInfo$cellType %in% c("HeLa_sync","HeLa_syncMix") ),"cid"] #27
length(cid.norm) #29
length(cid.sync) #27



#
checkGeneCounts=function(gDF,gN){
  print(gDF[gN, ])
  gene1=gDF[gN,]$gene
  apa1=row.names(gDF[gN, ])
  
  #print(deg[gene1,])
  #
  #gene
  df.g=data.frame(
    counts=c(
      as.numeric(rna.counts.norm[gene1,cid.norm]),
      as.numeric(rna.counts.norm[gene1,cid.sync])
    ),
    type=c(rep('normal', length(cid.norm)),   rep('sync', length(cid.sync)) )
  )
  
  #library(ggplot2) 
  g=ggplot(df.g, aes(type,log10(counts+1),color=type))+
    theme_bw()+
    geom_boxplot()+geom_jitter(size=0.5, alpha=1)+
    scale_color_manual(values=c('grey','#93BBFD'))+
    labs(title=gene1, x='', y="log10(RNAcounts+1)"); #g
  #apa
  df.a=data.frame(
    counts=c(
      as.numeric(APA.counts.norm[apa1,cid.norm]),
      as.numeric(APA.counts.norm[apa1,cid.sync])
    ),
    type=c(rep('normal', length(cid.norm)),   rep('sync', length(cid.sync)) )
  )
  #library(ggplot2) 
  g2=ggplot(df.a, aes(type,log10(counts+1),color=type))+
    theme_bw()+
    geom_boxplot()+geom_jitter(size=0.5, alpha=1)+
    scale_color_manual(values=c('grey','#93BBFD'))+
    labs(title=paste0( apa1  ,'(',gene1,')' ), x='', y="log10(RNAcounts+1)"); #g2
  grid.arrange(g2,g, nrow=1)
}

# 2 and 4
checkGeneCounts(dd_text.2, 1)
checkGeneCounts(dd_text.4s, 1)

t.test(APA.counts.norm['chr2:232577972:+',cid.norm], APA.counts.norm['chr2:232577972:+',cid.sync])
#

# 1 and 3
dim(dd_text.s1.1)
dd_text.s1.2

checkGeneCounts(dd_text.s1.2,1)
checkGeneCounts(dd_text.s3.2,1)
