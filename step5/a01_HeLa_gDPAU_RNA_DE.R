#aim: HeLa, RNA_DEG, APA_DEG

setwd("/data/jinwf/wangjl/apa/20200701Fig/f4/HeLaDE_apa_RNA")
getwd()


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

cid.sync=rownames(cellInfo[which( cellInfo$cellType=='HeLa_sync'),]); length(cid.sync) #27
cid.normal=rownames(cellInfo[which( cellInfo$cellType=='HeLa_normal'),]); length(cid.normal) #30

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

#(3) gDPAU
gDPAU=read.table('/home/wangjl/data/apa/20200701Fig/f4/DPAU/Matrix_01_generalDPAU.txt', row.names = 1)
gDPAU=gDPAU*100
dim(gDPAU) # 4103  222
gDPAU[1:5,1:5]
max(gDPAU, na.rm = T) #100
min(gDPAU, na.rm = T) #0


###################
# process1: filter out low expression genes
###################
gene.highRNA=(function(){
  dt1=data.frame(
    gene=row.names(rnaM.logcpm),
    cv2=apply(rnaM.logcpm, 1, function(x){
      (sd(x)/mean(x) )**2
    }),
    mean=apply(rnaM.logcpm, 1, mean)
  )
  g1=ggplot(dt1, aes(mean, cv2))+geom_point(size=0.1)+
    geom_vline(xintercept = 0.5,color="red", linetype=2)+
    theme_bw()+labs(title="mRNA noise~mean plot")
  pdf("01_mRNA_noise_mean_plot.pdf", width=3, height=3)
  print(g1)
  dev.off()
  print('end==')
  
  dt1[which(dt1$mean>=0.5),]
})()
dim(gene.highRNA) # 7819    3
head(gene.highRNA)
#       gene       cv2     mean
#AAAS   AAAS 2.2037756 1.673852
#AACS   AACS 3.5674072 1.267416
write.table(gene.highRNA, "01_highRNA_gene.df.txt")
gene.highRNA=read.table('01_highRNA_gene.df.txt', row.names = 1)



###################
# process2: get DEG_mRNA of HeLa, use all genes
###################
DEG_mRNA=read.csv("/home/wangjl/data/apa/20200701Fig/f2/HeLa_DEG/DESeq2_ALL_sync_VS_normal_HeLa.csv",
                  row.names = 1)

dim(DEG_mRNA) #15991     6
head(DEG_mRNA)
#          baseMean log2FoldChange     lfcSE      stat       pvalue         padj
#MESP1     18.98323      -24.73741 1.3359321 -18.51697 1.506960e-76 1.622393e-72

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
#483 14972   536
head(dif)
pdf('02_DEG_mRNA_volcano.pdf', width=3.5, height=2.8)
ggplot(dif, aes(log2FC, -log10(p.adj), color=factor(thresh, 
                                                    levels=c('up','n.s.','down')) ))+
  geom_point(size=0.2)+theme_bw()+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=1)))+
  scale_color_manual('', values=c('red', 'grey', 'blue'))+
  labs(title='mRNA DEG, sync vs normal')
dev.off()


# check: how many high expression gene?
rs1=intersect(rownames(gene.highRNA), rownames(DEG_mRNA))
length(rs1) #7819
dim(DEG_mRNA) #15991     6
dim(gene.highRNA) #7819    3




###################
# process3: get DEG_APA of HeLa, use all genes
###################
gDPAU=read.table('/data/jinwf/wangjl/apa/20200701Fig/f4/cycles/sync_normal/02_volcano_gDPAU_chisqP_sync_normal.txt',
                 row.names = 1)
dim(gDPAU) #3940   16
gDPAU[1:4,]

# volcano
dif2=gDPAU
dif2$p.adj=dif2$padj
dif2$log2FC=log2(dif2$fc)
dif2$thresh=dif2$sig #"n.s."
dif2$thresh=as.character(dif2$thresh)
#dif[which(dif$p.adj<0.05 & dif$log2FC>log2(2)),]$thresh="up"
#dif2[which(dif2$thresh=='ns'),'thresh'] ="n.s."
tb=table(dif2$thresh);tb
#down  n.s.    up 
#203 3526  211
head(dif)

pdf('03_DEG_APA_volcano.pdf', width=3.5, height=2.8)
ggplot(dif2, aes(log2FC, -log10(p.adj), color=factor(thresh, 
                                                     levels=c('up','ns','down')) ))+
  geom_point(size=0.2)+theme_bw()+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=1)))+
  scale_color_manual('', values=c('red', 'grey', 'blue'))+
  labs(title='APA DEG, sync vs normal')
dev.off()


# check: how many high expression gene?
rs2=intersect(rownames(gene.highRNA), rownames(gDPAU))
length(rs2) #3766
dim(DEG_mRNA) #15991     6
dim(gDPAU) #3940   16


rs3=intersect(rs1, rs2); length(rs3) #3766
length(rs1)  #10038



###################
# process4: ana x rna
###################

newDt=data.frame(
  gene=rs3,
  rnaFC=DEG_mRNA[rs3, 'log2FoldChange'],
  rnaPadj=DEG_mRNA[rs3, 'padj'],
  apaFC=log2( gDPAU[rs3, 'fc'] ),
  apaDelta=gDPAU[rs3, 'delta'],
  sigRNA=dif[rs3, 'thresh'],
  sigAPA=dif2[rs3, 'thresh']
)
row.names(newDt)=newDt$gene
dim(newDt) #3766    5
head(newDt)
#             gene      rnaFC       apaFC  apaDelta sigRNA sigAPA
#AARS         AARS  1.9638733 -0.24566915 -6.595807     up     ns
#AASDHPPT AASDHPPT -0.6211772  0.68939308 25.616201   n.s.     ns
#AATF         AATF -0.3379638 -0.15774448 -9.986034   n.s.     ns

table(newDt$apaDelta>30, newDt$sigAPA)
table(newDt$apaDelta<(-30), newDt$sigAPA)


table(newDt$sigRNA, newDt$sigAPA)
#     down   ns   up
#down   7   100    6
#n.s.  186 3249  192
#up      6   95    6
# output the gene symbol
newDt[which(newDt$sigRNA=="up" & newDt$sigAPA=="up"),'gene'] #1. GRB10   MAP2    MRPS22  NPM1P27 PI4KA   SSBP2 
newDt[which(newDt$sigRNA=="up" & newDt$sigAPA=="down"),'gene'] #2 DESI2    EMP1     HIST1H4C HIST3H2A INTS6    NUAK1
newDt[which(newDt$sigRNA=="down" & newDt$sigAPA=="down"),'gene'] #3 LINC00842 MAPK8     RDH10     SLC39A11  TSTD3     UNC50     WDR55
newDt[which(newDt$sigRNA=="down" & newDt$sigAPA=="up"),'gene'] #4 ARHGDIB  DYNC1LI1 KIAA1147 KRCC1    LEMD2    MRPL49 
#
writeLines( as.character(newDt[which(newDt$sigRNA!="n.s." & newDt$sigAPA!="ns"),]$gene), "geneList/sig_RNA_or_APA.gene.txt" )




table( newDt$rnaFC>log2(1.5), newDt$apaDelta>20) #  93
table( newDt$rnaFC>log2(1.5), newDt$apaDelta<=(-20) ) #126
#
table( newDt$rnaFC< (-log2(1.5)), newDt$apaDelta>20) #160
table( newDt$rnaFC< (-log2(1.5)), newDt$apaDelta<=(-20) ) #128
#

# define sig as log2(1.5)
newDt$sigRNA2="ns"
newDt[which(newDt$rnaPadj<0.05 & newDt$rnaFC>log2(1.5)),"sigRNA2"]="up"
newDt[which(newDt$rnaPadj<0.05 & newDt$rnaFC<(-log2(1.5)) ),"sigRNA2"]="down"
#
table( newDt$sigRNA2 ) 
#down   ns   up 
#148 3566  133 
table(newDt$sigRNA2, newDt$sigAPA)


write.table(newDt, '04_foldChangePlot_APA_RNA.df.txt')

#
pdf('04_foldChangePlot_APA_RNA.pdf', width=4,height=3.2)
ggplot(newDt, aes(apaDelta, rnaFC, color=factor(sigRNA, levels=c("up",'n.s.', 'down')) ) )+
  geom_point(alpha=0.8,size=0.3)+
  scale_color_manual("", values = c('red','grey','blue'),
                     labels=c('up', 'n.s.','down') )+
  theme_bw()+ylim(-5,5)+
  #xlim(-5,5)+
  geom_hline(yintercept = c(-1,1), linetype=2, color='#555555')+
  #geom_vline(xintercept = c(-0.26,0.26), linetype=2, color='#555555')+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=1)))+
  labs(x="Distal Usage Change of APA level",
       y="log2(Fold Change of mRNA level)",
       title="Sync vs normal HeLa")
dev.off()
#




# RNA as up, down
# define sig by RNA: |log2FC|>log2(2), FDR<0.05;
df_text=newDt[which(newDt$sigRNA!="n.s." & newDt$sigAPA!="ns"),]
dim(df_text) #25 8
head(df_text)

library(ggrepel)
g=ggplot(newDt, aes(apaDelta, rnaFC, color=factor(sigRNA, levels = c("up",'n.s.','down')) ) )+
  geom_point(alpha=0.5,size=1)+
  scale_color_manual("", values = c('red','grey','blue'),
                     labels=c('Up', 'n.s.','Down') )+
  theme_bw()+ylim(-7, 7)+
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
#


# plot again, only threshold no P cutoff.
newDt2=newDt
head(newDt2)
newDt2$RNAthresh="no"
newDt2[which(newDt2$rnaFC>log2(2) & abs(newDt2$apaDelta)>30 ),"RNAthresh"]='up'
newDt2[which(newDt2$rnaFC< (-log2(2)) & abs(newDt2$apaDelta)>30 ),"RNAthresh"]='down'
#
#df_text=newDt2[which(newDt2$RNAthresh!="ns"),] # by threshold
#df_text=newDt2[which(newDt2$sigAPA!="ns" & newDt2$sigAPA!="n.s."),] #by P and threshold
df_text=newDt2[which(newDt2$sigAPA!="ns" & newDt2$RNAthresh!="ns"),] # by RNA threshold

dim(df_text) #200 9
# df_text=df_text[order(-abs(df_text$apaDelta) ),] #order
df_text=df_text[order(-abs(df_text$rnaFC) ),] #order

#
g=ggplot(newDt2, aes(apaDelta, rnaFC, color=factor(RNAthresh, levels = c("up",'no','down')) ) )+
  geom_point(alpha=0.5,size=1)+
  scale_color_manual("", values = c('red','grey','blue'),
                     labels=c('Up', 'n.s.','Down') )+
  theme_bw()+ylim(-15, 10)+
  #xlim(-5,5)+
  geom_hline(yintercept = c(-1,1), linetype=2, color='#555555')+
  geom_vline(xintercept = c(-30,30), linetype=2, color='#555555')+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=3)))+
  labs(x="Defference of distal polyA usage\nsync-normal", #expression( paste(Delta, "DPAU", sep="" ) ),
       y="log2(Fold Change of expression)",
       title="Sync vs normal HeLa"); g
#
g2=g+geom_text_repel(data=df_text[1:15,], aes(x=apaDelta, y=rnaFC, label=gene), alpha=1,
                     colour="black",size=2.5)
pdf('04_foldChangePlot_APA_RNA-3.pdf', width=4,height=3.5)
g2
dev.off()
#
length( rownames( newDt2[which(newDt2$rnaFC>1 & newDt2$apaDelta>30),]) ) #33
length( rownames( newDt2[which(newDt2$rnaFC>1 & newDt2$apaDelta<(-30) ),]) ) #48
length( rownames( newDt2[which(newDt2$rnaFC <(-1) & newDt2$apaDelta>30),]) ) #62
length( rownames( newDt2[which(newDt2$rnaFC<(-1) & newDt2$apaDelta<(-30) ),]) ) #57

#


# left: shorten
barDf=data.frame(
  x=c('up','down'),
  y=c(48,57)
)
barDf$label=paste0(barDf$x, "\n", barDf$y)

p1=ggplot(barDf, aes(x, y, fill=factor(x, levels = c('up','down'))))+
  geom_bar(stat = "identity")+
  scale_fill_manual("",values=c('red','blue'), labels=NULL)+
  labs(title="", x="", y="")+
  theme_classic()+
  geom_text( mapping=aes(x, y=y+7, label=label, angle=0 ), size=2 )+
  scale_y_continuous(expand = c(0,0), limits=c(0,70), breaks = seq(0,70,20))+
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
  y=c(33,62)
)
barDf2$label=paste0(barDf2$x, "\n", barDf2$y)

p2=ggplot(barDf2, aes(x, y, fill=factor(x, levels = c('up','down'))))+
  geom_bar(stat = "identity")+
  scale_fill_manual("",values=c('red','blue'), labels=NULL)+
  labs(title="", x="", y="")+
  theme_classic()+
  geom_text( mapping=aes(x, y=y+10, label=label, angle=0 ), size=2 )+
  scale_y_continuous(expand = c(0,0), limits=c(0,80), breaks = seq(0,60,20))+
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
# width\height fig size; 
# x\y ratio position of subplot in main fig, relative to "just"
print(p1,vp=vp)
#
vp2 <- viewport(width = 0.2, height = 0.4, x = 0.6,y = 0.5,just=c("left","top"))
print(p2,vp=vp2)
#
















## venn plot of RNA and APA: up down
library("VennDiagram")
### 
A <- rownames( dif[which(dif$thresh=="up"),] ); length(A)
B <- rownames( dif[which(dif$thresh=="down"),] ); length(B)
C <- rownames( dif2[which(dif2$thresh=="up"),] ); length(C)
D <- rownames( dif2[which(dif2$thresh=="down"),] ); length(D)
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
  filename =NULL, #Venn_4set_pretty.pdf",
  col = "transparent", 
  fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
  alpha = 0.50, 
  label.col = c("orange", "white", "darkorchid4", "white",
                "white", "white", "white", "white", "darkblue", "white",
                "white", "white", "white", "darkgreen", "white"),
  cex = 1.5,    #label size
  fontfamily = "serif", 
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.cex = 1.5, 
  cat.pos = 0, 
  cat.dist = 0.07, 
  cat.fontfamily = "serif", 
  rotation.degree = 0, 
  margin = 0.2 
);
plot(NULL, ann=F)
grid.draw(venn.plot)

#

#










# test threshold
ct=(function(t1=log2(2)){
  set1=newDt[which(newDt$rnaFC>t1 & newDt$sigAPA=='up'),]; dim( set1) #26 5
  # left: down
  set2=newDt[which(newDt$rnaFC>t1 & newDt$sigAPA=='down'),]; dim(set2) #39 5
  #
  # right: up
  set3=newDt[which(newDt$rnaFC<=-t1 & newDt$sigAPA=='up'),];dim(set3) #51 5
  # right: down
  set4=newDt[which(newDt$rnaFC<=-t1 & newDt$sigAPA=="down"),];dim(set4) #48 5
  
  # switch
  if(1==1){
    # write genes
    writeLines( rownames(set1), "04_set1.gene.txt")
    writeLines( rownames(set2), "04_set2.gene.txt")
    
    writeLines( rownames(set3), "04_set3.gene.txt")
    writeLines( rownames(set4), "04_set4.gene.txt")
  }
  
  c(nrow(set1), nrow(set2), nrow(set3), nrow(set4))
})();ct
#ct=c(26, 39, 51, 48)


#
dt2=data.frame(
  length=c("lengthen","lengthen", '', "shorten","shorten"),
  exp=c('up','down', '', 'up','down'),
  #value=c(54,84,NA,85,83)
  value=c(ct[1:2],NA, ct[3:4])
)
dt2
#ggplot(dt2, aes(length, value, color=exp))+geom_histogram()

#
pdf("05_barplot_geneNumber_RNA_APA.pdf", width=2.8, height=3.8)
par(mgp=c(1.4,0.5,0) )
posX=barplot( rev(dt2$value), col=rev(c("#E82364",'#5BC4BA','white',"#E82364",'#5BC4BA')),
              ylim=c(0,60), border=NA, ylab="Gene number" )
text(x=posX, y=rev(dt2$value)+4, labels=rev(dt2$value) )
#text(x=posX-0.3, y=-6, offset = 1,
#     srt = 30, xpd = TRUE, cex=0.9,
#     labels=rev(dt2$exp) )
# arrows
arrows(3.5,-5,6,-5, length = 0.08,
       col='black', xpd = T, code=2, lwd=2)
text(x=5, y=-12,  xpd = TRUE, cex=0.9,
     labels="Lengthen" )
#
arrows(0,-5,2.5,-5, length = 0.08,
       col='black', xpd = T, code=1, lwd=2)
text(x=1, y=-12,  xpd = TRUE, cex=0.9,
     labels="Shorten" )
#
legend(-2,-15,xpd = TRUE, horiz=T,
       #box.lwd=3,
       text.width=2.5, 
       x.intersp=0.2,
       border=NA, cex=0.8,
       #title="RNA", title.adj=1,
       fill=c('white','#5BC4BA',"#E82364"), legend = c('RNA level','down', 'up'), bty='n' )
dev.off()

