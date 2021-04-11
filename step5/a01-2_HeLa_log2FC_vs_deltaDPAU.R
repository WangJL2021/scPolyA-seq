# log2FC vs delta DPAU(using cell cycle info): G1S+S

setwd("/data/jinwf/wangjl/apa/20200701Fig/f4/log2FC_deltaDPAU/")
getwd()


keyword="_HeLa_G1S_S_"


library(ggplot2)

pP=function(...){   print(paste(...)) }
pP0=function(...){   print(paste0(...)) }
pd=function(...){
  print(dim(...))
}
pL=function(x, say=''){print( paste( say, 'length=', length(x))  )}
###################
# load data
###################
#(1) cell info
cellInfo=read.table("/data/jinwf/wangjl/apa/20200701Fig/f3/cell_cycle/cellInfo.V6.txt")
dim(cellInfo) #222   8
head(cellInfo)
table(cellInfo$cellType)
#BC_0        BC_1 HeLa_normal   HeLa_sync 
#92          73          30          27
table(cellInfo$cellType, cellInfo$cellCycle)
#            G1S G2M  M MG1  S
#BC_0         29  20 13  20 10
#BC_1         18   9 12  19 15
#HeLa_normal  10   3  3  11  3
#HeLa_sync     3   6  7   0 11



cid.sync=rownames(cellInfo[which( cellInfo$cellType=='HeLa_sync' & 
                                  cellInfo$cellCycle %in% c("G1S","S") ),]); length(cid.sync) #14
cid.normal=rownames(cellInfo[which( cellInfo$cellType=='HeLa_normal' &
                                      cellInfo$cellCycle %in% c("G1S","S")),]); length(cid.normal) #13

#(2) add RNA read counts for each genes
rnaM=read.csv('/data/jinwf/wangjl/apa/20200701Fig/f2/BC_HeLa/BC_HeLa.222cells.count.V4.csv', row.names = 1)
dim(rnaM) #18662 222
rnaM[1:10,1:10]
#
rnaM.logcpm=as.data.frame(apply(rnaM, 2, function(x){
  log2(x/sum(x)*1e6+1)
}))
dim(rnaM.logcpm)
rnaM.logcpm[1:10,1:5]

# get sub
rnaM.logcpm=rnaM.logcpm[, c(cid.sync, cid.normal )]
dim(rnaM.logcpm) #18662    27

# rm all 0 rows
table(apply(rnaM.logcpm, 1, sum)>0)
rnaM.logcpm=rnaM.logcpm[which(apply(rnaM.logcpm, 1, sum)>0), ]
dim(rnaM.logcpm)

#(3) gDPAU
gDPAU=read.table('/home/wangjl/data/apa/20200701Fig/f4/DPAU/Matrix_01_generalDPAU.txt', row.names = 1)
gDPAU=gDPAU*100
dim(gDPAU) # 4103  222
gDPAU[1:5,1:5]
max(gDPAU, na.rm = T) #100
min(gDPAU, na.rm = T) #0

# get sub 
gDPAU=gDPAU[, c(cid.sync, cid.normal )]
dim(gDPAU) #4103   27

# at leat 2 cells in sync or normal
df=data.frame(
  gene=rownames(gDPAU),
  syncNum=apply(gDPAU, 1, function(x){
    sum( is.na(x[1:14]) )
  } ),
  normNum=apply(gDPAU, 1, function(x){
    sum( is.na(x[15:27]) )
  } )
)
head(df)
#
df2=df[which(df$syncNum>1 & df$normNum>1),]
dim(df2) #2977    3
head(df2)




###################
# process1: filter out low expression genes
###################
gene.highRNA=(function(){
  thresh=1
  
  dt1=data.frame(
    gene=row.names(rnaM.logcpm),
    cv2=apply(rnaM.logcpm, 1, function(x){
      (sd(x)/mean(x) )**2
    }),
    mean=apply(rnaM.logcpm, 1, mean)
  )
  g1=ggplot(dt1, aes(mean, cv2))+geom_point(size=0.1)+
    geom_vline(xintercept = thresh, color="red", linetype=2)+
    theme_bw()+labs(title="mRNA noise~mean plot")
  pdf("01_mRNA_noise_mean_plot_G1S_S_HeLa.pdf", width=3, height=3)
  print(g1)
  dev.off()
  print('end==')
  
  dt1[which(dt1$mean>=thresh),]
})()
dim(gene.highRNA)
head(gene.highRNA)
#       gene       cv2     mean
#AAAS   AAAS 2.2037756 1.673852
#AACS   AACS 3.5674072 1.267416
write.table(gene.highRNA, "01_highRNA_gene_G1S_S_HeLa.df.txt")
gene.highRNA=read.table('01_highRNA_gene_G1S_S_HeLa.df.txt', row.names = 1)




###################
# process2: get DEG_mRNA of HeLa(G1S+S), use high exp genes
###################

library(DESeq2)

# input raw counts
countData <- cbind( rnaM[,cid.sync], rnaM[,cid.normal] )
countData=countData[ rownames(gene.highRNA),]
dim(countData) #7914 27

condition <- factor(c( rep('sync', length(cid.sync)),
                       rep('normal', length(cid.normal)) ) )
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition )
# filter
nrow(dds) #7914
dds2 <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)
#one step method
dds3 <- DESeq(dds2) #time consuming: 1min
#
# get result
res <- results(dds3)
head(res)
# log2 fold change (MLE): condition sync vs normal 
# Wald test p-value: condition sync vs normal 

# assign padj NA as 1
res$padj[is.na(res$padj)] = 1

hist(res$padj, n=100)
# order by adjp
res <- res[order(res$pvalue),]
#set cutoff
resSig <- subset(res, abs(log2FoldChange)>log2(1.5) & padj < 0.05)
dim(resSig) #[1] 614  6
head(resSig)

# save to file
resSig<-data.frame(resSig)
dim(resSig)
write.csv(resSig, file=paste0("02_DESeq2_DEG_",keyword,".csv") )

#save all, for GSEA
res2=data.frame(res)
head(res2)
dim(res2)
write.csv(res2, file=paste0("02_DESeq2_ALL_",keyword,".csv") )
#



#
DEG_mRNA=res2
dim(DEG_mRNA) #7914     6
head(DEG_mRNA)
#          baseMean log2FoldChange     lfcSE      stat       pvalue         padj
#EFL1    43.44250      -24.65834 1.7711439 -13.92227 4.639529e-44 3.528825e-40

# vocano
dif=DEG_mRNA
#
dif$p.adj=dif$padj
dif$log2FC=dif$log2FoldChange
dif$thresh="n.s."
dif[which(dif$p.adj<0.05 & dif$log2FC>log2(2)),]$thresh="up"
dif[which(dif$p.adj<0.05 & dif$log2FC< ( -log2(2) ) ),]$thresh="down"
tb=table(dif$thresh);tb
#down  n.s.    up 
#267  7319    328 

write(rownames(dif[which(dif$thresh=="up"), ]), '02_RNA_up.gene.txt')
write(rownames(dif[which(dif$thresh=="down"), ]), '02_RNA_down.gene.txt')


head(dif)

g=ggplot(dif, aes(log2FC, -log10(p.adj), color=factor(thresh, 
                                                    levels=c('up','n.s.','down')) ))+
  geom_point(size=0.2)+theme_bw()+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=1)))+
  scale_color_manual('', values=c('red', 'grey', 'blue'))+
  labs(title=paste0('mRNA DEG, sync vs normal\n',keyword ) )
g
#add gene symbols
dd_text = dif[ ((abs(dif$log2FoldChange) > 4) & (dif$padj < 0.05e-30) ) | 
                 (abs(dif$log2FoldChange) > 20 & dif$padj < 0.05e-20),]; dim(dd_text)

head(dd_text)
library(ggrepel)
g2=g + geom_text_repel(data=dd_text, aes(x=log2FoldChange, y=-log10(padj), label=row.names(dd_text)), 
                       size=2.5, colour="black",alpha=0.6); g2
# save
pdf('02_DEG_mRNA_volcano.pdf', width=3.5, height=3)
print(g)
print(g2)
dev.off()
#





###################
# process3: get DEG_APA of HeLa G1S+S, use High Exp genes
###################

dim(gDPAU)

genes=intersect(rownames(gDPAU), rownames(gene.highRNA) )
length(genes) #3801
gDPAU=gDPAU[genes,]

#
# /home/wangjl/data/apa/20200701Fig/f4/script/f4_1_get_DPAU_gDPAU.R
####
# HeLa: sync vs normal
# normalize gRatio to (0,1)
# method 1(log, (x-min)/(max-min) ),
# method 2(original raw ratio)
# method 3: x/(x+1)
compareA_B=function(dt, cid.A, cid.B){
  dt=dt[, c(cid.A, cid.B)]
  #filter out 0
  keep=apply(dt, 1, function(x){ sum(x, na.rm = T)})>0
  dt=dt[keep,]
  pd(dt)
  
  #
  dt=as.data.frame(dt)
  
  genes= rownames(dt);    
  df=NULL;
  j=0
  for(i in 1:length(genes)){
    gene=genes[i]
    arrA=as.numeric(dt[gene,cid.A])
    arrB=as.numeric(dt[gene,cid.B])
    # filter out 0 ratio
    arrA=arrA[!is.na(arrA)]
    arrB=arrB[!is.na(arrB)]
    #
    if(length(arrA)==0 | length(arrB)==0){
      pP('>>>Discard(either lenght is 0):',i, ' ', gene)
      next;
    }
    p=1
    if(length(arrA)>=2 & length(arrB)>=2){
      if(sd(arrA)!=0 & sd(arrB)!=0){
        p=wilcox.test(arrA, arrB)$p.value;
        j=j+1
      }
    }
    df=rbind(df, data.frame(
      #gene=gene,
      nA=length(arrA), nB=length(arrB),
      meanA=mean(arrA),meanB=mean(arrB),
      sdA=sd(arrA),sdB=sd(arrB),
      p=p,
      row.names=gene
    ))
  }
  pP("=========>>>j=",j)
  #rownames(df)=df$gene
  df$cvA=df$sdA / df$meanA
  df$cvB=df$sdB / df$meanB
  #
  df$adj.p=p.adjust(df$p, method ='fdr')
  df=df[order(df$adj.p),]
  df$fc=df$meanA / df$meanB
  df$delta=df$meanA - df$meanB
  df
}

# compare: sync vs normal, HeLa G1S+S
gDPAU.sync_normal=compareA_B(gDPAU, cid.sync, cid.normal)
dim(gDPAU.sync_normal)
head(gDPAU.sync_normal)

table(gDPAU.sync_normal$adj.p<0.05)
# only 1: EIF3K

dif2=gDPAU.sync_normal
head(dif2)
dif2$thresh="ns"
dif2[which(dif2$delta>30), ]$thresh="lengthen"
dif2[which(dif2$delta<(-30) ), ]$thresh="shorten"
#

delta_gDPAU_df=(function(){
  df=gDPAU.sync_normal
  df$delta=df$delta
  #plot(log2(df$cvA), log2(df$cvB) )
  #plot(df$meanA, df$meanB)
  hist(df$delta, n=100)
  
  rownames(df)
  m1=mean(df$delta)
  s1=sd(df$delta)
  pP(m1, s1)
  t1=m1+2*s1
  t2=m1-2*s1
  pP(t1, t2)
  abline(v=t1, col='red', lty=2)
  abline(v=t2, col='red', lty=2)
  tmp=rbind(
    df[which(df$delta>t1),],
    df[which(df$delta<t2),]
  )
  tmp
})()

dim(delta_gDPAU_df) # 249 12
head(delta_gDPAU_df)

writeLines(rownames(delta_gDPAU_df), "03_APA_changed.gene.txt")

hist(gDPAU.sync_normal$delta, n=100)
write.table(gDPAU.sync_normal, paste0("03_APA_changed_",keyword,'.txt') )










###################
# process4: ana x rna
###################

genes=intersect(rownames(dif), rownames(gDPAU.sync_normal))
length(genes) #3611

newDt=data.frame(
  gene=genes,
  rnaFC=dif[genes, 'log2FoldChange'],
  rnaP=DEG_mRNA[genes, 'pvalue'],
  rnaPadj=DEG_mRNA[genes, 'padj'],
  apaFC=log2( gDPAU.sync_normal[genes, 'fc'] ),
  apaDelta=gDPAU.sync_normal[genes, 'delta'],
  sigRNA=dif[genes, 'thresh']
)
row.names(newDt)=newDt$gene
dim(newDt)
head(newDt)

hist(newDt$apaDelta, n=100)

table(newDt$apaDelta>30)
table(newDt$apaDelta<(-30))


write.table(newDt, '04_foldChangePlot_APA_RNA.df.txt')



# Adjust p again, by the current row numbers.
dif_=newDt
head(dif_)
dif_$rnaPadj=p.adjust(dif_$rnaP, method = "fdr")
#
newDt=dif_
newDt$sigRNA="n.s."
newDt[which(newDt$rnaPadj<0.05 & newDt$rnaFC> 1), ]$sigRNA="up"
newDt[which(newDt$rnaPadj<0.05 & newDt$rnaFC< (-1) ), ]$sigRNA="down"
table(newDt$sigRNA)
#down n.s.   up 
#73 3475   63
writeLines(rownames(newDt[ which(newDt$sigRNA!="n.s."),]),  "geneList/new_RNA_sig.gene.txt")
#
pdf('04_foldChangePlot_APA_RNA.pdf', width=4,height=3.2)
ggplot(dif_, aes(apaDelta, rnaFC, color=factor(sigRNA, levels=c("up",'n.s.', 'down')) ) )+
  geom_point(alpha=0.8,size=0.3)+
  scale_color_manual("", values = c('red','grey','blue'),
                     labels=c('up', 'n.s.','down') )+
  theme_bw()+#ylim(-5,5)+
  #xlim(-5,5)+
  geom_hline(yintercept = c(-1,1), linetype=2, color='#555555')+
  geom_vline(xintercept = c(-30,30), linetype=2, color='#555555')+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=1)))+
  labs(x="Distal Usage Change of APA level",
       y="log2(Fold Change of mRNA level)",
       title="Sync vs normal HeLa")
dev.off()




# RNA as up, down
# define sig by RNA: |log2FC|>log2(2), FDR<0.05;
df_text=newDt[which(newDt$sigRNA!="n.s." & abs(newDt$apaDelta)>30 ),]
dim(df_text) #31 8
head(df_text)

library(ggrepel)
g=ggplot(newDt, aes(apaDelta, rnaFC, color=factor(sigRNA, levels = c("up",'n.s.','down')) ) )+
  geom_point(alpha=0.5,size=1)+
  scale_color_manual("", values = c('red','grey','blue'),
                     labels=c('Up', 'n.s.','Down') )+
  theme_bw()+#ylim(-7, 7)+
  #xlim(-5,5)+
  geom_hline(yintercept = c(-1,1), linetype=2, color='#555555')+
  geom_vline(xintercept = c(-20,20), linetype=2, color='#555555')+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=3)))+
  labs(x="Defference of distal polyA usage\nsync-normal", #expression( paste(Delta, "DPAU", sep="" ) ),
       y="log2(Fold Change of expression)",
       title="Sync vs normal HeLa"); g
#
pdf('04_foldChangePlot_APA_RNA-2.pdf', width=4,height=3.5)
g+geom_text_repel(data=df_text, aes(x=apaDelta, y=rnaFC, label=gene), alpha=1,
                  colour="black",size=2.5)
dev.off()





# plot again, only threshold no P cutoff.
newDt2=newDt
head(newDt2)
newDt2$RNAthresh="ns"
newDt2[which(newDt2$rnaFC>log2(2) & abs(newDt2$apaDelta)>30 ),"RNAthresh"]='up'
newDt2[which(newDt2$rnaFC< (-log2(2)) & abs(newDt2$apaDelta)>30 ),"RNAthresh"]='down'
table(newDt2$RNAthresh)
#
#
#df_text=newDt2[which(newDt2$RNAthresh!="ns"),] # by threshold
#df_text=newDt2[which(newDt2$sigAPA!="ns" & newDt2$sigAPA!="n.s."),] #by P and threshold
df_text=newDt2[which( newDt2$RNAthresh!="ns"),] # by RNA threshold
dim(df_text)
# df_text=df_text[order(-abs(df_text$apaDelta) ),] #order
df_text=df_text[order(-abs(df_text$rnaFC) ),] #order

#
g=ggplot(newDt2, aes(apaDelta, rnaFC, color=factor(RNAthresh, levels = c("up",'ns','down')) ) )+
  geom_point(alpha=0.5,size=1)+
  scale_color_manual("", values = c('red','grey','blue'),
                     labels=c('Up', 'n.s.','Down') )+
  theme_bw()+ylim(-15, 10)+
  #xlim(-5,5)+
  geom_hline(yintercept = c(-1,1), linetype=2, color='#555555')+
  geom_vline(xintercept = c(-30,30), linetype=2, color='#555555')+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=3)))+
  labs(x="Difference of distal polyA usage\nsync-normal", #expression( paste(Delta, "DPAU", sep="" ) ),
       y="log2(Fold Change of expression)",
       title="Sync vs normal HeLa"); g
#
g2=g+geom_text_repel(data=df_text[1:15,], aes(x=apaDelta, y=rnaFC, label=gene), alpha=1,
                     colour="black",size=2.5);g2
pdf('04_foldChangePlot_APA_RNA-3.pdf', width=4,height=3.5)
g2
dev.off()
#
length( rownames( newDt2[which(newDt2$rnaFC>1 & newDt2$apaDelta>30),]) ) #75
length( rownames( newDt2[which(newDt2$rnaFC>1 & newDt2$apaDelta<(-30) ),]) ) #96
length( rownames( newDt2[which(newDt2$rnaFC <(-1) & newDt2$apaDelta>30),]) ) #85
length( rownames( newDt2[which(newDt2$rnaFC<(-1) & newDt2$apaDelta<(-30) ),]) ) #69
#


# left: shorten
barDf=data.frame(
  x=c('up','down'),
  y=c(96,69)
)
barDf$label=paste0(barDf$x, "\n", barDf$y)

p1=ggplot(barDf, aes(x, y, fill=factor(x, levels = c('up','down'))))+
  geom_bar(stat = "identity")+
  scale_fill_manual("",values=c('red','blue'), labels=NULL)+
  labs(title="", x="", y="")+
  theme_classic()+
  geom_text( mapping=aes(x, y=y+7, label=label, angle=0 ), size=2 )+
  scale_y_continuous(expand = c(0,0), limits=c(0,120), breaks = seq(0,100,20))+
  theme(legend.position="none",
        
        axis.text = element_text(size = rel(0.5)),
        
        # no bg
        panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA),
        
        # no x axis
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank());p1


# right: lengthen
barDf2=data.frame(
  x=c('up','down'),
  y=c(75,85)
)
barDf2$label=paste0(barDf2$x, "\n", barDf2$y)

p2=ggplot(barDf2, aes(x, y, fill=factor(x, levels = c('up','down'))))+
  geom_bar(stat = "identity")+
  scale_fill_manual("",values=c('red','blue'), labels=NULL)+
  labs(title="", x="", y="")+
  theme_classic()+
  geom_text( mapping=aes(x, y=y+10, label=label, angle=0 ), size=2 )+
  scale_y_continuous(expand = c(0,0), limits=c(0,100), breaks = seq(0,100,20))+
  theme(legend.position="none",
        
        axis.text = element_text(size = rel(0.5)),
        
        # no bg
        panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA),
        
        # no x axis
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank());p2

#
print(g2)
library(grid)
vp <- viewport(width = 0.2, height = 0.4, x = 0.08,y = 0.5,just=c("left","top"))
# width\height表示插入图形的大小，x\y表示插入图形相对于图片底层just的位置
print(p1,vp=vp)
#
vp2 <- viewport(width = 0.2, height = 0.4, x = 0.6,y = 0.5,just=c("left","top"))
print(p2,vp=vp2)







## venn plot of RNA and APA: up down
library("VennDiagram")
 
A <- rownames( dif[which(dif$thresh=="up"),] ); length(A)
B <- rownames( dif[which(dif$thresh=="down"),] ); length(B)
C <- rownames( dif2[which(dif2$thresh=="lengthen"),] ); length(C)
D <- rownames( dif2[which(dif2$thresh=="shorten"),] ); length(D)
#
writeLines(A, "geneList/RNA_up.gene.txt")
writeLines(B, "geneList/RNA_down.gene.txt")
writeLines(C, "geneList/APA_lengthen.gene.txt")
writeLines(D, "geneList/APA_shorten.gene.txt")

#
g1=intersect(rownames(dif), rownames(dif2))
length(g1)
df3=data.frame(
  gene=g1, 
  tRNA=dif[g1, "thresh"],
  tAPA=dif2[g1, 'thresh']
)
table(df3$tRNA, df3$tAPA)
#




#
dat=list(
  A = A,
  B = B,
  C = C,
  D = D
)
names(dat) <- c('up', 'down', 'lengthen','shorten')

venn.plot <- venn.diagram(
  #数据列表
  x = dat,
  filename =NULL, #Venn_4set_pretty.pdf",    #保存路径
  col = "transparent",      #指定图形的圆周边缘颜色  transparent 透明           
  fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),  #填充颜色
  alpha = 0.50,                                      #透明度
  label.col = c("orange", "white", "darkorchid4", "white",
                "white", "white", "white", "white", "darkblue", "white",
                "white", "white", "white", "darkgreen", "white"),
  cex = 1.5,    #每个区域label名称的大小
  fontfamily = "serif",  #字体
  fontface = "bold",     #字体格式
  cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),  #分类颜色 
  cat.cex = 1.5,      #每个分类名称大小
  cat.pos = 0,        #
  cat.dist = 0.07,    #
  cat.fontfamily = "serif",     #分类字体
  rotation.degree = 0,        #旋转角度
  margin = 0.2               #在网格单元中给出图周围空白量的编号
);
plot(NULL, ann=F)
grid.draw(venn.plot)

