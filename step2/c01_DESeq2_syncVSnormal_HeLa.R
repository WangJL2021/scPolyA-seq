# aim: DEG of sync vs normal HeLa


setwd("/data/jinwf/wangjl/apa/20200701Fig/f2/HeLa_DEG/")
keyword='sync_VS_normal_HeLa'

##########################
# 1 get DE gene list, with DESeq2
##########################
library(DESeq2)
# read matrix
rnaM=read.csv('/data/jinwf/wangjl/apa/20200701Fig/f2/BC_HeLa/BC_HeLa.222cells.count.V4.csv',
                    header = T,row.names = 1)
dim(rnaM) #18662   222
rnaM[1:8,1:10]

####
# load cell info
cellInfo=read.table('/data/jinwf/wangjl/apa/20200701Fig/f3/cell_cycle/cellInfo.V6.txt',header = T,row.names = 1)
dim(cellInfo) #222   7
head(cellInfo)

# order
table(colnames(rnaM) == row.names(cellInfo)) #222 T
#cellInfo[colnames(rna.counts),]
#
table(cellInfo$cellType)
# BC_0        BC_1 HeLa_normal   HeLa_sync 
#  92          73          30          27 


##########
## DEG:
cid.normal=row.names(cellInfo[cellInfo$cellType=='HeLa_normal',]);length(cid.normal) #30
cid.sync=row.names(cellInfo[cellInfo$cellType=='HeLa_sync',]);length(cid.sync) #27
#
countData <- cbind( rnaM[,cid.normal], rnaM[,cid.sync] )
condition <- factor(c( rep('normal', length(cid.normal)),     rep('sync', length(cid.sync)) ) )
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition )
# filter
nrow(dds) #18662
dds2 <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds2)#[1] 15991




#one step method
dds3 <- DESeq(dds2) #time consuming: 1min? 14:26 to 14:27



# get result
res <- results(dds3)
head(res)
#log2 fold change (MLE): condition sync vs normal 
#Wald test p-value: condition sync vs normal 
#DataFrame with 6 rows and 6 columns
#              baseMean     log2FoldChange             lfcSE               stat             pvalue              padj
#              <numeric>          <numeric>         <numeric>          <numeric>          <numeric>         <numeric>
#A1BG-AS1  0.134031507022109 -0.662348822587702  2.93875759303932 -0.225383959587728  0.821680589045084                NA
#A2M      0.0444748982089599   0.30714792231751  2.93949565545738  0.104490007238918  0.916780497938929                NA

# assign padj NA as 1
res$padj[is.na(res$padj)] = 1

# order by adjp
res <- res[order(res$pvalue),]
#set cutoff
resSig <- subset(res, abs(log2FoldChange)>log2(1.5) & padj < 0.05)
dim(resSig) #[1] 1109  6
head(resSig)
#                 baseMean    log2FoldChange             lfcSE              stat               pvalue                 padj
#               <numeric>         <numeric>         <numeric>         <numeric>            <numeric>            <numeric>
#MESP1    18.9832252819049 -24.7374075092372  1.33593206804302 -18.5169651219425 1.50695970625321e-76 1.62239281975221e-72
#SLC39A11 15.4296355120073 -24.6518145747351  1.43711436806004 -17.1536901464651  5.8985028554193e-66 3.17516408707221e-62
#SIPA1    12.9644828764506 -24.4212767346591   1.4946030319017 -16.3396408366615  5.1555699794424e-60 1.85016221328923e-56
#ACTA2    396.592633144536  10.6230937104073 0.655447129149554  16.2073998618177 4.47087762889729e-59  1.2033367138177e-55

# save to file
resSig<-data.frame(resSig)
dim(resSig) #[1] 1109   6
write.csv(resSig, file=paste0("DESeq2_DEG_",keyword,".csv") )

#save all, for GSEA
res2=data.frame(res)
head(res2)
write.csv(res2, file=paste0("DESeq2_ALL_",keyword,".csv") )





##########################
# 2 volcano plot
##########################
library('ggplot2')
dif=data.frame(res)
dif$threshold="ns"
dif[which( res$log2FoldChange > log2(2) & res$padj < 0.05 ),]$threshold= "up"
dif[which( res$log2FoldChange < (-log2(2)) & res$padj < 0.05 ),]$threshold= "down"
tb=table(dif$threshold);tb
#down    ns    up 
#483 14972   536
# dif$threshold= factor( abs(res$log2FoldChange) > log2(1.5) & res$padj < 0.05, levels=c(TRUE,FALSE) )
dif$threshold=factor(dif$threshold, levels=c('down','ns','up'))
str(dif)
head(dif)
#          baseMean log2FoldChange     lfcSE      stat       pvalue         padj threshold
#MESP1     18.98323      -24.73741 1.3359321 -18.51697 1.506960e-76 1.622393e-72      down
#SLC39A11  15.42964      -24.65181 1.4371144 -17.15369 5.898503e-66 3.175164e-62      down



g = ggplot(data=dif, aes(x=log2FoldChange, y=-log10(padj), color=threshold)) +
  geom_point(alpha=0.5, size=0.1) +
  #opts(legend.position = "none") + 
  theme_bw() +
  theme(legend.box = "horizontal",
        legend.margin=margin(t = -0.5, unit='line'), #图例整体上边距,缩减n行
        legend.spacing.x = unit(4, 'pt'), #图例之间的x距离
        legend.key.width=unit(0.5,"line"), #图例方块的宽度度
        plot.title = element_text(size=10),
        legend.position="bottom") +
  scale_shape(guide = guide_legend(title.position = "top")) +
  #scale_colour_discrete(guide = guide_legend(title.position = "top", nrow = 1))+
  scale_color_manual('Sig.',labels=c(
    paste0('down(',tb[1],')'),
    "ns",
    paste0('up(',tb[3],')')
  ), 
                     values=c("blue", "grey",'red') )+
  #xlim(c(-10, 10)) + ylim(c(0, 15)) +
  xlab("log2(FoldChange)") + ylab("-log10(p-value)") +
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=1))) +#放大图例的点
  labs(title= paste0("DEG: ",keyword) ); g
#add gene symbols
dd_text = dif[ ((abs(dif$log2FoldChange) > 4) & (dif$padj < 0.05e-30) ) | 
                 (abs(dif$log2FoldChange) > 20 & dif$padj < 0.05e-20),]; dim(dd_text)
head(dd_text)
library(ggrepel)
g2=g + geom_text_repel(data=dd_text, aes(x=log2FoldChange, y=-log10(padj), label=row.names(dd_text)), 
                       size=2.5, colour="black",alpha=0.6); g2
#save to file
#CairoPDF(file=paste0('DEG_volcano_plot_',keyword,'.pdf'), width=3.5,height=4)
CairoPDF(file=paste0('DEG_volcano_plot_',keyword,'-small.pdf'), width=3.2,height=3.5)
print(g2)
dev.off()
#

# save to fiel
dim(dif)
head(dif) #15991     7
write.csv(dif[which(dif$threshold!="ns"),], file=paste0("DESeq2_DEG_",keyword,".with_threshold.csv") )


length(rownames(dif[which(dif$threshold=='up'),])) #536
length(rownames(dif[which(dif$threshold=='down'),])) #483
writeLines(rownames(dif[which(dif$threshold=='up'),]), "02_DEG.up.gene.txt")
writeLines(rownames(dif[which(dif$threshold=='down'),]), "02_DEG.down.gene.txt")

#end






