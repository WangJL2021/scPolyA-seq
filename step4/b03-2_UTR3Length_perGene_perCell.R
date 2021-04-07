#f08_UTR3Length_perGene_perCell.R

setwd('/data/jinwf/wangjl/apa/190705PAS/UTR3_len/')

keyword="sync_vs_normal_HeLa"

library(ggplot2)
library(Cairo)

######################
#step1 get apa and gene list
#data from mysql
source('/home/wangjl/pylib/con_mysql.R')
query<-dbSendQuery(con, 'select pasID, an_gene, UTR3_len_an from feature_apa where UTR3_len_an is not null;');
data=fetch(query, n=-1);
row.names(data)=data$pasID
dim(data) #35230     3
head(data)
#                        pasID an_gene UTR3_len_an
#chr10:320134:- chr10:320134:-   DIP2C        3131
#chr10:320403:- chr10:320403:-   DIP2C        2862

#step2 get uniq gene list
geneList=unique(data$an_gene)
geneList[1:10]
length(geneList) #10633

#step3 get apa count matrix
apaCounts=read.table("/data/jinwf/wangjl/apa/190705PAS/bed/freq2/all225_matrix_APACountsV4.txt",header = T)
colnames(apaCounts)=sub('_','',colnames(apaCounts))
apaCounts[1:4,1:4]
dim(apaCounts) #[1] 267961    225
#
#apaID both in DB anno and matrix
apaID=intersect(row.names(apaCounts),data$pasID)
length(apaID) #35230
#
apaCounts2=apaCounts[data$pasID,]
dim(apaCounts2)  #35230   225










#step4, calculate 3UTR length weight by counts, per gene per cell
meanLengthDF=NULL
for(i in seq(1: length(geneList)) ){ # 耗时22:32 - 22:37
  if(i%%500==0){
    print(i)
  }
  #i=2
  # get gene symbol
  gene=geneList[i];#gene
  #apaID belong to this gene
  apaList.gene=data[which(data$an_gene==gene),]$pasID
  #apaID UTR3 length
  length.gene=data[which(data$an_gene==gene),]$UTR3_len_an
  #
  #mean length, weighted by counts at this PolyA site
  df=apply(apaCounts2[apaList.gene,],2,function(x){
    if(sum(x)==0) return(0);
    sum( x*length.gene)/sum(x)
  } )
  df=as.data.frame( t(df) )
  rownames(df)=gene
  #df
  meanLengthDF=rbind(meanLengthDF,df)
}
dim(meanLengthDF) #[1] 10633   225
meanLengthDF[1:4,1:8]
#        c01ROW07 c01ROW12 c01ROW24 c01ROW31 c01ROW35 c04ROW06 c04ROW14 c05ROW02
#DIP2C         0        0     2862     2862 2862.000        0        0     2862
#LARP4B     2105        0        0        0 2374.961        0        0        0
#ADARB2        0        0        0        0    0.000        0        0        0
#PITRM1      300        0        0        0  300.000        0        0      300
View(meanLengthDF[,c(1,8,20,30,50,100,200)] )


#
#write to file
getwd()
write.table(meanLengthDF,'01-UTR3_length_perGenePerCell.txt')
meanLengthDF0=meanLengthDF






################################
# DE length use LIMMA

#db conn
source('/home/wangjl/pylib/con_mysql.R')

#sync vs normal
query<-dbSendQuery(con, 'select * from cell_c1;');
cellInfo=fetch(query, n=-1);
dim(cellInfo) #[1] 225   9
head(cellInfo)
#  id      cid readCounts geneNumber     cellType cellTypeRough clusterRNA batch cellCycle
#1  1 c01ROW07    2543685       5815    HeLa_sync     HeLa_sync          2     2         M
#2  2 c01ROW12    2703371       6007    HeLa_sync     HeLa_sync          2     2       G2M

HeLa_sync.cid=cellInfo[which(cellInfo$cellTypeRough=='HeLa_sync'),]$cid;length(HeLa_sync.cid) #27
HeLa_normal.cid=cellInfo[which(cellInfo$cellTypeRough=='HeLa_normal'),]$cid;length(HeLa_normal.cid) #29


## at leat 5 cell expression, then calc length
tmp.1=data.frame(
  row.names=row.names(meanLengthDF),
  n.N=apply(meanLengthDF[, HeLa_normal.cid]>0,1,sum),
  n.S=apply(meanLengthDF[, HeLa_sync.cid]>0,1,sum)
)
head(tmp.1)
dim(tmp.1) #10633     2
tmp.2=tmp.1[which(tmp.1$n.N>5 & tmp.1$n.S>5),]
dim(tmp.2) #5040    2

#meanLengthDF0=meanLengthDF
meanLengthDF=meanLengthDF[rownames(tmp.2),]
dim(meanLengthDF) #5040  225

#compare the length
compareL=NULL;
var0Gene=NULL;
for(i in 1:nrow(meanLengthDF) ){
  if(i%%500==0){print(i)}
  gene=row.names(meanLengthDF)[i];gene
  # test
  #if( var(as.numeric( meanLengthDF[gene,HeLa_sync.cid]) )==0  & var( as.numeric(meanLengthDF[gene,HeLa_normal.cid]) )==0 ){
  #  var0Gene=c(var0Gene,gene)
  #  next;
  #}
  #(1) t test, the distribution not like normal
  #t.output=t.test(  as.numeric( meanLengthDF[gene,HeLa_sync.cid]), as.numeric( meanLengthDF[gene,HeLa_normal.cid]) )
  #p.value=t.output$p.value
  #
  #(2)
  test.output=wilcox.test(  as.numeric( meanLengthDF[gene,HeLa_sync.cid]), as.numeric( meanLengthDF[gene,HeLa_normal.cid]) )
  p.value=test.output$p.value
  #
  if(is.na( mean( as.numeric(meanLengthDF[gene,HeLa_normal.cid])) ) ){
    FC=NA
  }else{
    FC= mean( as.numeric(meanLengthDF[gene,HeLa_sync.cid]) ) / mean( as.numeric(meanLengthDF[gene,HeLa_normal.cid]) )
  }
  df=data.frame(
    id=gene,
    p.value=p.value,
    FC=FC
  )
  #df
  compareL=rbind(compareL,df)
}
dim(compareL) #[1] 5040    3
#
head(var0Gene) #NULL
length(var0Gene)
#gene=var0Gene[2];print(gene); wilcox.test(  as.numeric( meanLengthDF[gene,HeLa_sync.cid]), as.numeric( meanLengthDF[gene,HeLa_normal.cid]) )

##replace na with 1
table(is.na(compareL$p.value))
#FALSE  TRUE 
#5029    11 
compareL[is.na(compareL$p.value), ]
compareL[is.na(compareL$p.value), ]$p.value=1;

#order
compareL=compareL[order(compareL$p.value),]
head( compareL )





#adj p value
############
n=nrow(compareL);n #5040
qvalues=c()
for(i in seq(1,n)){
  if(i%%1000==0){
    print(i)
  }
  q=compareL$p.value[i]*n/i
  qvalues=c(qvalues, q)
}
length(qvalues) #
qvalues[1:10] 
# 1.237026e-05 2.178873e-05 9.390084e-05 2.327849e-04 1.435291e-03
#
############ R fn: test for adj.p value function
qvalues2=p.adjust(p=compareL$p.value, method = "fdr")
qvalues2[1:10]
#1.234326e-05 2.174118e-05 9.369590e-05 2.322768e-04 
## result: same as before.

compareL$q.value=qvalues2
compareL$logFC=log2(compareL$FC)

#add sig
compareL$threshold = as.factor( abs(compareL$logFC) > log2(1.5) & compareL$q.value < 0.05 )
compareL$threshold=factor(compareL$threshold, levels=c(T, F) )
table(compareL$threshold)
#TRUE FALSE 
#78  4962 
head(compareL)
str(compareL)

#
compareL$threshold2 = "ns"
compareL[which( compareL$logFC > log2(1.5) & compareL$q.value < 0.05 ),]$threshold2 ='up';
compareL[which( compareL$logFC < (-log2(1.5)) & compareL$q.value < 0.05 ),]$threshold2 ='down';
compareL$threshold2=factor(compareL$threshold2, levels=c('down','ns','up'))
tb2=table(compareL$threshold2);tb2
#down   ns   up 
#59 4962   19 

#save
head(compareL)
dim(compareL) #5040 7
getwd()
write.table(compareL,'02-UTR3_length_DE_HeLa_sync_vs_normal.txt')

#
gene.lengthen=as.character( compareL[which(compareL$threshold2=='up'),]$id ); length(gene.lengthen) #19
writeLines(gene.lengthen, paste0('03-UTR3_lengthen_genes_19.txt') )
#
gene.shorten=as.character( compareL[which(compareL$threshold2=='down'),]$id ); length(gene.shorten) #59
writeLines(gene.shorten, paste0('03-UTR3_shorten_genes_59.txt') )
#






################
#3. plot 3UTR length volcano plot
dif=compareL
g = ggplot(data=dif, aes(x=log2(FC), y=-log10(q.value), color=threshold2)) +
  geom_point(alpha=0.8, size=0.4) +
  theme_bw() +
  theme(legend.box = "horizontal", 
        legend.margin=margin(t = -0.5, unit='line'), 
        legend.spacing.x = unit(5, 'pt'),
        legend.key.size = unit(5, "pt"), 
        legend.position="bottom") +
  scale_color_manual("3'UTR",labels=c(paste0("shorten(",tb2[[1]],')'),'ns',
                                            paste0("lengthen(",tb2[[3]],')' )),
                     values=c("blue", "grey",'red') )+
  labs(x="log2(FoldChange)",y="-log10(q-value)",title= paste0("3'UTR: ",keyword) ); g
#
# add text to a few genes
dd_text = dif[ which( ((abs( log2(dif$FC) ) > log2(2) ) & (dif$q.value < 0.01) ) | 
                 (abs( log2(dif$FC) ) > 3 & dif$q.value<0.05 ) ),]; dim(dd_text); #17
dd_text$id=as.character(dd_text$id)
head(dd_text)
#View(dd_text)
library(ggrepel)
g2=g + geom_text_repel(data=dd_text, aes(x=log2(FC), y=-log10(q.value), label=dd_text$id),
              size=2.5, colour="black",alpha=0.4)+xlim(-3.5,2.8);g2
#
CairoPDF(file=paste0('04-volcano_plot_3UTR_',keyword,'.pdf'), width=3.6,height=4)
print(g2)
dev.off()






############
#check1
#
checkByGene=function(gene){
  #gene='MMP7'
  print( meanLengthDF[gene,HeLa_normal.cid])
  print('after sync:')
  print(meanLengthDF[gene,HeLa_sync.cid])
}
checkByGene(dd_text[4,'id'])

#
head(compareL)
dim(dd_text) #17 7
compareL[which(compareL$id==dd_text[5,'id']),]
compareL[which(compareL$id==dd_text[17,'id']),]
#
head(compareL)
table(compareL$threshold)
# TRUE FALSE 
# 78  4962
table(compareL$threshold2)
#down   na   up 
#59 4962   19


compareL[which(compareL$id==dd_text[2,'id']),] #
#   id      p.value       FC        q.value     logFC     threshold threshold2
#  295 LTBP3 8.646322e-09 2.363135 2.178873e-05 1.240702      TRUE         up


compareL[which(compareL$id==dd_text[10,'id']),]
#         id      p.value       FC     q.value    logFC threshold threshold2
#GRB10 GRB10 1.671573e-05 2.871443 0.005616487 1.521776      TRUE         up

compareL[which(compareL$id==dd_text[4,'id']),]
#         id      p.value       FC      q.value    logFC threshold threshold2
#SPARC SPARC 1.847499e-07 2.118285 0.0002327849 1.082897      TRUE         up

compareL[which(compareL$id==dd_text[5,'id']),]
compareL[which(compareL$id==dd_text[17,'id']),]



# down
compareL[which(compareL$id==dd_text[3,'id']),]
#           id      p.value        FC      q.value    logFC threshold threshold2
#CLDN10 CLDN10 5.589336e-08 0.1656591 9.390084e-05 -2.59371      TRUE       down

checkByGene( 'SPARC' )

gene='RAC1'
compareL[which(compareL$id==gene),]
checkByGene(gene)

#
gene='CCND1'
checkByGene(gene)
compareL[which(compareL$id==gene),]

#
gene='CDK2'
checkByGene(gene)
compareL[which(compareL$id==gene),]

#

#
t1=as.numeric( ( meanLengthDF[gene,HeLa_normal.cid]) )
t2=as.numeric( (meanLengthDF[gene,HeLa_sync.cid]) )
t.test(t1,t2)
wilcox.test(t1,t2)


#不行，
#left: sync only one cell have length: 'TSTD1','LCN2','EFHD1','LRRC1'
#rigth: normal only 1-2 cells have length: 'FBN2','MAP2','CCND2'
gene='MMP7'
meanLengthDF[gene,HeLa_normal.cid]
meanLengthDF[gene,HeLa_sync.cid]
#


#########################################
# check 2
# IGV


##############
# 5. DEG and UTR length
dim(compareL) #5040    7
row.names(compareL)=compareL$id
head(compareL)
#        id      p.value        FC      q.value      logFC threshold threshold2
#3151    FN1 2.454416e-09 4.8333333 1.237026e-05  2.2730185      TRUE         up
#295   LTBP3 8.646322e-09 2.3631349 2.178873e-05  1.2407020      TRUE         up
#

# read deg
deg=read.csv("/home/wangjl/data/apa/191111Figure/f5/DEG_RNA/DESeq2_ALL_sync_VS_normal_HeLa.csv",header = T,row.names = 1)
dim(deg) #7487 6
head(deg)
#         baseMean log2FoldChange     lfcSE       stat       pvalue         padj
#ACTA2   395.05595      10.538784 0.6231554  16.911969 3.672328e-64 2.642607e-60
#S100A7  165.70262      -9.186365 0.6914571 -13.285516 2.809139e-40 1.010728e-36

# common genes
gene2=intersect( row.names(deg), compareL$id  );length(gene2) #4806
#
mtL=compareL[gene2,]
mtD=deg[gene2,]
head(mtL)
head(mtD)
#
df2=data.frame(
  row.names = gene2,
  #
  lengthening=mtL$FC,
  L.qvalue=mtL$q.value,
  L.thresh=mtL$threshold2,
  #
  geneLogFC=mtD$log2FoldChange,
  D.padj=mtD$padj
)
df2$g.thresh="ns";
df2[which( df2$geneLogFC>1.5 & df2$D.padj<0.05 ), ]$g.thresh="up";
df2[which( df2$geneLogFC<(-1.5) & df2$D.padj<0.05 ), ]$g.thresh="down";
df2$g.thresh=factor(df2$g.thresh, levels=c('up', 'ns', 'down'))
table(df2$g.thresh) #
#up   ns down 
#74 4695   37 

# new loose standard, define lengthen
df2$L.thresh="ns";
df2[which(  log2(df2$lengthening) > log2(1.5)   & df2$L.qvalue<0.1 ), ]$L.thresh="lengthen";
df2[which(  log2(df2$lengthening) < (-log2(1.5)) & df2$L.qvalue<0.1 ), ]$L.thresh="shorten";
#
df2$L.thresh=factor(df2$L.thresh, levels=c('lengthen','ns','shorten' ))
table(df2$L.thresh) #
#lenthen      ns shorten 
# 34     4630      142 

# SEC61G   0.9812221 1.073529e-02       ns -2.031922 2.757264e-05     down     1 SEC61G

table(df2$g.thresh, df2$L.thresh)
#    lengthen   ns shorten
#up          9   65       0
#ns         25 4538     132
#down        0   27      10

# long-up
row.names(df2[which(df2$g.thresh=='up' & df2$L.thresh=="lengthen"),]) #9
#[1] "ENO3"   "FBN1"   "TUBA1A" "NUAK1"  "FN1"    "DDR2"   "GRB10"  "FOXO3"  "COPG1" 
# short-down
row.names(df2[which(df2$g.thresh=='down' & df2$L.thresh=="shorten"),]) #10
#[1] "PDCL3"   "KRT19"   "S100A2"  "EGFR"    "UBE2J2"  "CLDN7"   "HSPBP1"  "AP1M2"   "S100A14" "SRSF1"  
#
writeLines(row.names(df2[which(df2$g.thresh=='up' & df2$L.thresh=="lengthen"),]), '05-up-lengthen.txt')
writeLines(row.names(df2[which(df2$g.thresh=='down' & df2$L.thresh=="shorten"),]), '05-down-shorten.txt')



# order, place sig in front
df2$order=0;
df2[which(df2$L.thresh!='ns' | df2$g.thresh!='ns'), ]$order=1;

df2=df2[order(df2$order),]
head(df2)


# up-lengthen, down-shorten plot:
g=ggplot(df2, aes(log2(lengthening), geneLogFC, color=g.thresh  ) )+geom_point(alpha=0.8, size=0.1)+
  theme_bw()+
  theme(#legend.box = "vertical", 
        legend.margin=margin(t = -0.5, unit='line'),
        legend.spacing.x = unit(5, 'pt'), 
        legend.key.size = unit(5, "pt"), 
        legend.position="bottom") +
  guides(col = guide_legend(ncol = 1, 
                            # inset=-0.5,
                            byrow=T))+ 
  scale_color_manual("Differential \nexpression",
                     #labels=c(paste0("shorten(",tb2[[1]],')'),'ns',paste0("lengthen(",tb2[[3]],')' )),
                     labels=c('Up regulated', 'No change','Down regulated'),
                     values=c("red", "grey",'blue') )+
  geom_hline(aes(yintercept=log2(1.5)), colour="#666666", linetype="dashed")+
  geom_hline(aes(yintercept=-log2(1.5)), colour="#666666", linetype="dashed")+
  #
  geom_vline(aes(xintercept=log2(1.5)), colour="#666666", linetype="dashed")+
  geom_vline(aes(xintercept=-log2(1.5)), colour="#666666", linetype="dashed")+
  scale_x_discrete(limits=df2$order, labels = NULL)+ 
  labs(x="log2(3'UTR length fold change)",y="log2(Gene expression fold change)",title= paste0("3'UTR length and gene FC") ); 
g
#
# select some genes: FC large and p small
dd_text2 = df2[which(df2$g.thresh!='ns' & abs(df2$lengthening) > log2(1.5) & df2$L.qvalue < 0.1 ),]; dim(dd_text2) #12 7
dd_text2$gene=row.names(dd_text2)
head(dd_text2)
#
g2=g+geom_text_repel(data=dd_text2, aes(x=log2(lengthening), y=geneLogFC, label= gene), alpha=0.5,
                  colour="black",size=3);g2

CairoPDF(file=paste0('05-plot_DE_3UTRLength_',keyword,'.pdf'), width=3.6,height=4)
print(g2)
dev.off()


