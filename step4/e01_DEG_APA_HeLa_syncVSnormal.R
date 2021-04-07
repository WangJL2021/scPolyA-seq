# HeLa sync vs normal cells. (APA data)
# DE APA and GO plot
# https://blog.csdn.net/hfcao_bnu/article/details/52848363
setwd('/home/wangjl/data/apa/191111Figure/f4/DEG_APA/')
getwd()

# load
cellInfo=read.table("/data/jinwf/wangjl/apa/20200701Fig/f3/cell_cycle/cellInfo.V6.txt", 
                  row.names = 1, stringsAsFactors = F)
head(cellInfo)
#              cid cellType0 cluster cellType cellType2 countsPerCell geneNumber cellCycle
# c12ROW03 c12ROW03   unknown       0     BC_0      BC_0       1820666       6659       MG1
# c12ROW04 c12ROW04   unknown       0     BC_0      BC_0       2816130       7930       MG1
# c12ROW05 c12ROW05   unknown       0     BC_0      BC_0       2109777       6846         M

# APA matrix
APA.counts=read.csv('/data/jinwf/wangjl/apa/190705PAS/bed/apa_Matrix_afterFilterByCounts_cell_20222.csv',
                      header = T,row.names = 1)
dim(APA.counts) #20222   225
#colnames(APA.counts)= sub("_","", colnames(APA.counts) )
APA.counts[1:8,1:10]
#

#get hela cell
table(cellInfo$cellType)
#BC_0        BC_1 HeLa_normal   HeLa_sync 
#92          73          30          27
HeLa.id=cellInfo[which(cellInfo$cellType %in% c('HeLa_normal', 'HeLa_sync') ),'cid']
HeLa.id
br.data=APA.counts[, HeLa.id]
dim(br.data) #267961    57
br.data[1:4,1:4]
#
br.data2=br.data




#### NO Use
# get gene info
query<-dbSendQuery(con, 'select * from feature_gene;');
geneInfo=fetch(query, n=-1);
dim(geneInfo) #[1] 60627    15
head(geneInfo)
#  id              gene_id   chr source     start       end strand      gene_type gene_status gene_name level tag
#1  1 ENSG00000000003.14_2  chrX HAVANA  99882106  99894988      - protein_coding                TSPAN6     2    
#2  2  ENSG00000000005.6_3  chrX HAVANA  99839933  99854882      + protein_coding                  TNMD     2   
#
#  stop_codon nextGenePos   utr_end
#1   99885795    99665271  99881106
#2   99854714    99899192  99855882



# load apa info
apaInfo=read.table("/data/jinwf/wangjl/apa/190705PAS/bed/pasPostions_Location_transcriptName-noChrM_noInnerPrime_PY-filterCountCell-motif-Header.bed",
                   row.names = 1, header = T)
head(apaInfo)
table(apaInfo$gene=="")
#FALSE 
#20222


dim(br.data2) #20222    57
br.data2[1000:1005,1:10]
# all 0 rows
table(apply(br.data2, 1, sum)==0)
#FALSE  TRUE 
#20179    43 
# filter out all 0 rows
br.data2=br.data2[apply(br.data2, 1, sum)!=0,]
dim(br.data2) #20179    57

#save
write.table(br.data2, file="01-apa.HeLa.matrix_150626_56.txt")






# 2. divide to 2 group: sync vs normal.
# discard all 0 rows in either group; keep genes >=5 cell expressed.
cid.norm=cellInfo[which(cellInfo$cellType=="HeLa_normal"),"cid"] #30
apa.normal=br.data2[, cid.norm]
dim(apa.normal) #[1] 20179    30

cid.sync=cellInfo[which(cellInfo$cellType %in% c("HeLa_sync") ),"cid"] #27
apa.sync=br.data2[, cid.sync]
dim(apa.sync) #[1] 20179    27

apa.HeLa=br.data2[, c(cid.norm, cid.sync)]
dim(apa.HeLa) #20179    27




##########################
# 1 get DE gene list, with DESeq2
##########################
library(DESeq2)

# matrix
countData <- apa.HeLa
dim(countData) #20179    27
countData[1:4,1:4]
# keep genes >=5 cell expression
df=data.frame(
  gene=row.names(countData),
  num.N=apply(countData[, cid.norm]>0,1,sum),
  num.S=apply(countData[, cid.sync]>0,1,sum),
  row.names = 1
)
head(df)
df2=df[which(df$num.N>5 & df$num.S>5),]
dim(df2) #8457    2
#
countData=countData[rownames(df2),]
dim(countData) #[1] 8457    57

condition <- factor(c( rep('normal', length(cid.norm)), rep('sync', length(cid.sync))  ), 
                    levels=c('normal','sync'))
condition
#
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition )
# 
nrow(dds)
dds2 <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds2)

keyword='sync_VS_normal_HeLa_APA'
#
#
dds3 <- DESeq(dds2) #1min
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
## -- replacing outliers and refitting for 5164 genes
## -- DESeq argument 'minReplicatesForReplace' = 7 
## -- original counts are preserved in counts(dds)
## estimating dispersions
## fitting model and testin

#get result
res <- results(dds3)
head(res)
#log2 fold change (MLE): condition sync vs normal 
#Wald test p-value: condition sync vs normal 
#DataFrame with 6 rows and 6 columns
#baseMean     log2FoldChange             lfcSE              stat             pvalue              padj
#<numeric>          <numeric>         <numeric>         <numeric>          <numeric>         <numeric>
#chr10:320403:-  24.7996983833884   1.45583037737819  0.96463214508465  1.50920782061481    0.1312456822689 0.536363473552064
#chr10:1085971:- 55.8440517623289  -1.54994565965727 0.687224704120947 -2.25536953250227 0.0241101434880722 0.210927662491458
#chr10:1705831:- 22.7976953084813   2.69560799769561  1.41433607889631  1.90591758063554 0.0566608997265875 0.353337759654335


# replace p-value NA as 1
res$padj[is.na(res$padj)] = 1

# order
res <- res[order(res$padj),]
dim(res) #8457 6

# set cutoff
resSig <- subset(res, abs(log2FoldChange)>log2(1.5) & padj < 0.05)
dim(resSig) #[1] 428  6
head(resSig)

#save to file
resSig<-data.frame(resSig)
dim(resSig) #[1] 428   6
head(resSig)
write.csv(resSig, file=paste0("DESeq2_DE_APA_",keyword,".csv") )


#save all file for GSEA
res2=data.frame(res)
dim(res2) #8766 6
head(res2)
write.csv(res2, file=paste0("DESeq2_ALL_",keyword,".csv") )




##########################
# 2 volcano plot
##########################
library('ggplot2')
dif=data.frame(res)
dif$threshold= factor( abs(dif$log2FoldChange) > log2(1.5) & dif$padj < 0.05, #1.5, 0.05
                       levels=c(TRUE,FALSE) )
str(dif)
head(dif)
tb=table(dif$threshold);tb
#TRUE FALSE 
#428  8029

#
dif$threshold2="ns";
dif[which(dif$log2FoldChange > log2(1.5) & dif$padj < 0.05),]$threshold2="up";
dif[which(dif$log2FoldChange < (-log2(1.5)) & dif$padj < 0.05),]$threshold2="down";
dif$threshold2=factor(dif$threshold2, levels=c('down','ns','up'))
tb2=table(dif$threshold2);tb2
#down   ns   up 
#155 8029  273

# save up and down gene list
geneUp=row.names(dif[which(dif$threshold2=='up'),]);length(geneUp) #268
head(geneUp)
writeLines(geneUp, paste0('DESeq2_', keyword,'_polyASite_UP.txt') )
#
geneDown=row.names(dif[which(dif$threshold2=='down'),]);length(geneDown) #153
head(geneDown)
writeLines(geneDown, paste0('DESeq2_', keyword,'_polyASite_DOWN.txt') )

##############
g = ggplot(data=dif, aes(x=log2FoldChange, y=-log10(padj), color=threshold2)) +
  geom_point(alpha=0.4, size=0.1) +
  theme_bw() +
  theme(legend.box = "horizontal",
        legend.position="bottom") +
  #scale_color_manual('Significant',labels=c("TRUE","FALSE"), values=c("red", "grey") )+ 
  scale_color_manual('Significant',labels=c(paste0("down(",tb2[[1]],')'),'ns',
                                            paste0("up(",tb2[[3]],')' )),
                     values=c("blue", "grey",'red') )+
  xlab("log2(FoldChange)") + ylab("-log10(p.adj)") +
  labs(title= paste0("DE: ",keyword) ); g
# add text to a few genes
dd_text = dif[ ((abs(dif$log2FoldChange) > 4) & (dif$padj < 1e-45) ) | 
                 abs(dif$log2FoldChange) > 7,]; dim(dd_text)
head(dd_text)
#add text
library(ggrepel)
g2=g + geom_text_repel(data=dd_text, aes(x=log2FoldChange, y=-log10(padj), label=row.names(dd_text)), 
                       size=2.5, colour="#333333",alpha=0.6); g2
#
CairoPDF(file=paste0('volcano_plot_',keyword,'.pdf'), width=3.6,height=4)
print(g2)
dev.off()
#

##########################
# 3 check the counts
##########################
showCounts=function(gene1){
  print(gene1)
  c0=as.numeric(br.data2[gene1,cid.norm]);
  c1=as.numeric(br.data2[gene1,cid.sync]);
  deltaCV=sd(c1)/mean(c1)-sd(c0)/mean(c0)
  df=data.frame(
    counts=c(c0,c1),
    type=c(rep('normal', length(c0)),   rep('sync', length(c1)) )
  )
  #library(ggplot2)
  g=ggplot(df, aes(type,log10(counts+1),color=type))+
    theme_bw()+
    geom_boxplot()+geom_jitter(size=0.5, alpha=1)+
    scale_color_manual(values=c('grey','#93BBFD'))+
    labs(title=gene1,x=paste0("deltaCV:",round(deltaCV,2) ) , y="log10(RNAcounts+1)")
  g
}
dd_text=dd_text[order(-dd_text$log2FoldChange),]
head(dd_text)
dim(dd_text) #27
geneUp2=rownames(dd_text[which(dd_text$log2FoldChange>0),] );length(geneUp2) #26
geneDown2=rownames(dd_text[which(dd_text$log2FoldChange<0),] );length(geneDown2) #1
#
CairoPDF(file=paste0("02-Check_counts_", keyword,".pdf"),width=7,height=2.5)
grid.arrange(
  showCounts(geneUp2[1]),
  showCounts(geneUp2[2]),
  #showCounts(geneUp2[13]),
  
  showCounts(geneDown2[1]),
  #showCounts(geneDown2[7]),
  #showCounts(geneDown2[8]),
  nrow=1
)
dev.off()
#











##########################
# 4 GO and plot
##########################
library('GOstats')
library('org.Hs.eg.db')

head(resSig)
dim(resSig) #421 6

head(apaInfo)
index2=match( row.names(resSig), apaInfo$PASid)
head(index2)
resSig$gene=apaInfo[index2, ]$gene
head(resSig)

#
DEapa=rownames(resSig)
head(DEapa)# "chr6:153986769:+" "chr2:143857505:-" 

#
DEgenes2=resSig$gene
length(DEgenes2) #421
DEgenes=unique(DEgenes2)
length(DEgenes) #350
#
head(DEgenes) #"MIR7641-2" "KYNU"      "SPATA17"   "CORIN"     "BCHE"      "CXCR4"



####
#gene symbole to id 
# https://www.biostars.org/p/255657/
library(org.Hs.eg.db)
library(biomaRt)
hs <- org.Hs.eg.db
#my.symbols <- c("MESP1")
my.symbols=DEgenes
head(my.symbols) #"MIR7641-2" "KYNU"      "SPATA17"   "CORIN"     "BCHE"      "CXCR4" 
rs=select(hs, 
          keys = my.symbols,
          columns = c("ENTREZID", "SYMBOL"),
          keytype = "SYMBOL")
head(rs)
#        SYMBOL  ENTREZID
#1         MESP1     55897
#
rs2=rs[!is.na(rs$ENTREZID),]
dim(rs2) #339 2
head(rs2$ENTREZID) #"8942"   "128153" "10699"



#GO
ego <- enrichGO(OrgDb="org.Hs.eg.db", 
                gene = rs2$ENTREZID, 
                ont = "ALL", 
                pvalueCutoff = 0.05, 
                readable= TRUE)
head(ego)
dim(ego) #307  10
table(ego$ONTOLOGY)
#BP  CC  MF 
#223  67  17

#
write.csv(summary(ego),paste0('GO-enrich_',keyword,'.csv'),row.names =F)

# 可视化
dotplot(ego,showCategory=20,title=paste0("Enrichment GO Top20",'\n', keyword) )
barplot(ego, showCategory=20,title="EnrichmentGO")  #柱状图


#
ego_BP <- enrichGO(OrgDb="org.Hs.eg.db", 
                   gene = rs2$ENTREZID, 
                   ont = "BP", 
                   pvalueCutoff = 0.05, 
                   readable= TRUE) #GO富集分析

#very slow
#BiocManager::install('topGO') #
library(topGO)
CairoPDF(file=paste0('GO_',keyword,'_plotGOgraph.pdf'), width=10,height=10)
plotGOgraph(ego_BP)
dev.off()
#
#文字清晰，无遮挡
CairoPDF(file=paste0('GO_',keyword,'_goplot.pdf'), width=12,height=8)
#
goplot(ego_BP)
dev.off()

#
CairoPDF(file=paste0('GO_',keyword,'_emapplot.pdf'), width=8,height=6)
emapplot(ego, showCategory = 10)
dev.off()

#
CairoPDF(file=paste0('GO_',keyword,'_cnetplot.pdf'), width=10,height=10)
cnetplot(ego, showCategory = 5)
dev.off()





#####################
head(resSig)
dim(resSig) #428 6
#order
rs1=resSig
rs1=rs1[order(rs1$log2FoldChange),]
head(rs1)
dim(rs1) #428   7

########## up APA, down gene
gUp= rs1[which(rs1$log2FoldChange>0),]$gene; length(gUp) #273
gDn= rs1[which(rs1$log2FoldChange<0),]$gene; length(gDn) #155
# no na
table(is.na(gUp))
table(is.na(gDn))
# de dup
gUp2=unique(gUp);length(gUp2) #225
gDn2=unique(gDn); length(gDn2) #130

# 6 overlap
gCommon=intersect(gUp2,gDn2); length(gCommon) #6
gCommon=gCommon[which(gCommon!=".")]
gCommon #"PTMA"   "DNAJA1"    "ACTG1"  "RPS7"   "RPL28"
# 


# save
writeLines(gUp2, paste0('DESeq2_', keyword,'_gene_UP.txt') )
writeLines(gDn2, paste0('DESeq2_', keyword,'_gene_DOWN.txt') )
writeLines(gCommon, paste0('DESeq2_', keyword,'_gene_UP_DOWN-5.txt') )



#############
## browse up and down gene expression
dim(resSig) #428 7
head(resSig)
#                 baseMean log2FoldChange     lfcSE     stat       pvalue         padj      gene
#chr6:153986769:+ 647.5427       9.455781 0.5019714 18.83729 3.735904e-79 3.159454e-75 MTCO2P31
#chr2:143857505:- 507.8949       8.851964 0.5191605 17.05053 3.463539e-65 1.464557e-61  MTCO3P5

scanByGene=function(i=1){
  tmp=resSig[ grep(gCommon[i], resSig$gene),];
  print(dim(tmp)) #39 7
  tmp
}

scanByGene(1) #39(32+, 7-)
scanByGene(2) #3(1+, 2-)
scanByGene(3) #2
scanByGene(4) #2
scanByGene(5) #2
