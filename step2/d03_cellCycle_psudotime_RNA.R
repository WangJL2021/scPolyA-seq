# update version

setwd('/data/jinwf/wangjl/apa/20200701Fig/f3/cell_cycle/psudoTimeRNA2/')
getwd()


#1. tools
pP=function(...){  print(paste(...)) }
pP0=function(...){  print(paste0(...)) }
ph=function(...){  print(head(...)) }
pd=function(...){  print(dim(...)) }
pl=function(...){  print(length(...)) }
debug=function(...){
  pd(...)
  ph(...)
}


#2. tools
getDF_fromNamesXX=function(namedXX){
  rs=data.frame(
    name=names(namedXX),
    value=as.numeric(namedXX)
  )
  row.names(rs)=rs$name
  rs=rs[order(-rs$value),]
  print(dim(rs))
  print(head(rs))
  return(rs)
}



################
# step1 load data
################
#(1) load rnaM
rnaM=read.csv('/data/jinwf/wangjl/apa/20200701Fig/f2/BC_HeLa/BC_HeLa.222cells.count.V4.csv', row.names = 1)
dim(rnaM) #18662   222
rnaM[1:5,1:10]

# RNA normalization
rnaM.cpm=apply(rnaM, 2, function(x){
  1e6*x/sum(x)
})
getNormalizedCts <- function ( cts ) {
  #cts <- read.table ( ctsPath , header = T , as.is = T )
  apply ( cts , 2 , function ( x ) { log2 ( ( 10^6 ) * x / sum ( x ) + 1 ) })
}
# rnaM.log2cpm= apply(rnaM.cpm, 2, function(x){
#     log2(x+1)
# }) #
rnaM.log2cpm=getNormalizedCts(rnaM)
dim(rnaM.log2cpm) #18662   222
rnaM.log2cpm[1:5,1:8]



#(2) load cellInfo
cellInfo=read.table('/data/jinwf/wangjl/apa/20200701Fig/f3/cell_cycle/cellInfo.V6.txt',header = T,row.names = 1)
dim(cellInfo)
head(cellInfo)
#              cid cellType0 cluster cellType cellType2 countsPerCell geneNumber cellCycle
#c12ROW03 c12ROW03   unknown       0     BC_0      BC_0       1820666       6659       MG1
#c12ROW04 c12ROW04   unknown       0     BC_0      BC_0       2816130       7930       MG1

cid.HeLa=rownames(cellInfo[which( substr(cellInfo$cellType,1,2)=="He" ),]); head(cid.HeLa); length(cid.HeLa); #57
cid.BC=rownames(cellInfo[which( substr(cellInfo$cellType,1,2)=="BC" ),]); head(cid.BC); length(cid.BC); #165


table(cellInfo$cellType)


#(3) load psudu time
library(pheatmap)
getCellOrder=function(cellType, keyword){
  #cellType for figure
  # keyword for filename
  df0=read.table(paste0('../result/',keyword,'_PhaseRefCor.txt') )
  df=df0[, 1:10]
  head(df)
  #library(pheatmap)
  pheatmap(t(df)[], border_color = NA,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           main=paste(cellType, nrow(df)) )
  return(df0)
}
#cellOrder.HeLa=getCellOrder('HeLa','HeLa2')
cellOrder.BC=getCellOrder('BC','BC')




################
# step2 RNA expression levels along cell cycle timeline
################
getPhase=function(){
  c('G1S',"S","G2M","M","MG1")
}
draw_psudoTime_RNA=function(cellOrder, keyword, gene, color1='red', index=0, 
                            rawColor='white', ### hide the raw lines
                            to1=F, ylim=c(-2,-1),isMain=T,sub="",...){
  if( !(keyword %in% c("BC", "HeLa") ) ){
    stop('Error: keyword must in c("BC", "HeLa")')
  }
  n=nrow(cellOrder) # cell number
  cid.order=rownames(cellOrder)
  cycleNames=as.character(cellOrder$assignedPhase)
  print(length(cid.order))
  #
  dt0=rnaM.log2cpm[, as.character(cid.order)];
  dt0.1=as.data.frame(t( apply(dt0,1,function(x){ (x-min(x)) / (max(x)-min(x)) }) ) )

  # verticle lines
  cycle_df=as.data.frame(table(cycleNames))
  rownames(cycle_df)=cycle_df$cycleNames;
  cycle_posX=cumsum( cycle_df[getPhase(),]$Freq )
  
  getRNA=function(gene){
    dt2=dt0;
    if(to1){ dt2=dt0.1;}
    #
    as.numeric(dt2[gene, ])
  }
  dt=getRNA(gene)
  #
  minY=min(dt);  maxY=max(dt)
  if(ylim[1]>=0) minY=ylim[1];
  if(ylim[2]>0) maxY=ylim[2];
  if( is.na(ylim[3]) ){
    if(maxY<10)maxY=10;
  }
  indexChar=ifelse(index==0,'', paste0("; i=", index))
  main=ifelse(isMain==T, paste0('RNA (',keyword,': ', nrow(cellOrder),' cells), gene: ', gene, indexChar), '')
  plot( dt, type="l", xaxt="n", col=rawColor, ylim=c(minY,maxY),
        yaxs="i", 
        xaxs ="i",# no space between fig and axis
        mgp=c(1.5,0.5,0),
        xlab='Cell phase', ylab="RNA expression level(log2 scale)", 
        main="", 
        sub=sub, ... )
  # bg color
  cycle_colors=c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854")
  arr=c(0, cycle_posX)
  for(i in 1:5){
    start=arr[i]
    end=arr[i+1]
    L=end-start+1
    yy=c(rep(maxY,L), rep(0,L) )
    
    polygon( c(start:end, end:start), yy, col = paste0(cycle_colors[i], "66"), border = NA)
  }

  abline(v=cycle_posX[1:4], lty=2, col="grey") #vertical lines
  text(x= (c(0,cycle_posX[1:4])+cycle_posX)/2, y=par('usr')[3]-0.5, 
       col=c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"),
       labels = getPhase(), xpd=T  )
  
  # 
  addLines=function(y0, color,...){
    n=length(y0);
    #y=y0;
    y=c(y0,y0,y0)
    model2=loess(y ~ seq(1, 3*n), span=0.1) ########## critical parameter!!! as the array lenthenes, the span decreases.
    y2=predict( model2 )
    lines(y2[(n+1): (2*n)], type="l", col=color,...)
    #lines(y2, type="l", col=color,...)
  }
  addLines(dt, color1) #B in M
  return(function(geneName, color,...){
    addLines(getRNA(geneName), color,...)
  })
}
### test: to determin the optimal span
fn8.addLines=draw_psudoTime_RNA(cellOrder.BC, "BC", gene='CCNE2', 'red', to1=F, ylim=c(0,5)) #G1S


#####################
# pick a gene for each phase //use
#####################
pdf("04_RNA_curve_5genes.pdf", width=4.3,height=4.3)
# pic1: origin
fn8.addLines=draw_psudoTime_RNA(cellOrder.BC, "BC", gene='CCNE2', 'red', to1=F, ylim=c(0,13)) #G1S
geneArr5=c("CCNE2","CDK4",  "CCNB1", "CCNB2",'MKI67',"GAPDH")
fn8.addLines(geneArr5[1], 'red', lty=1, lwd=2)
fn8.addLines(geneArr5[2], 'red', lty=2, lwd=2)

fn8.addLines(geneArr5[3], 'blue', lwd=2)
fn8.addLines(geneArr5[4], 'blue', lty=2, lwd=2)

fn8.addLines(geneArr5[5], 'black', lwd=2)
fn8.addLines(geneArr5[6], 'black',lty=2, lwd=2)

legend(x=3, y=par()$usr[4]+2.5, lwd=2,
       legend=geneArr5, 
       lty=c(1,2,1,2,1,2),
       col=c('red','red','blue','blue','black','black'),
       ncol=4,#horiz=T,
       
       x.intersp=0.4, 
       text.width=28, 
       cex=0.9,
       xpd=T, bty='n')

# pic2
fn8.addLines=draw_psudoTime_RNA(cellOrder.BC, "BC", gene='CCNE2', 'red', to1=F, ylim=c(0,13)) #G1S
fn8.addLines('CCNE2', 'red', lty=1, lwd=2)
fn8.addLines('CCNE1', 'red', lty=2, lwd=2)

fn8.addLines('CDC20', 'blue', lwd=2)
fn8.addLines('CCNB1', 'blue', lty=2, lwd=2)

fn8.addLines('PTTG1', 'black', lwd=2)
fn8.addLines('GAPDH', 'black',lty=2, lwd=2)

legend(x=3, y=par()$usr[4]+2.5, lwd=2,
       legend=c("CCNE2", 'CCNE1', "CDC20", "CCNB1","PTTG1","GAPDH"), 
       lty=c(1,2,1,2,1,2),
       col=c('red','red','blue','blue','black','black'),
       ncol=4,#horiz=T,
       
       x.intersp=0.5, 
       text.width=28, 
       cex=0.9, 
       xpd=T, bty='n')
dev.off()



################
# draw all the genes' plots.
for(phase in getPhase()){
  print(phase)
  
  genes=readLines( 
    paste0("/data/jinwf/wangjl/apa/20200701Fig/f3/cell_cycle/result/BC_CycleRelatedGene_RNA_",phase,".txt"))
  
  pdf( paste0("03_All_gene_RNA_curve_",phase,".pdf"), width=6, height=7)
  par(mfrow=c(4,2), mar=c(1,2,3,0))
  for(i in 1:length(genes) ){
    gene=genes[i]
    #if(i >9)break;
    cat(phase, "|", i,'/',length(genes),":\t", gene)
    fn8.addLines=draw_psudoTime_RNA(cellOrder.BC, "BC", gene=gene, 'red', to1=F, ylim=c(0,10))
    title(main=paste(i,': ', gene))
  }
  dev.off()
}
#


#################
fn8.addLines=draw_psudoTime_RNA(cellOrder.BC, "BC", gene='CCNE2', 'red', to1=F, ylim=c(0,13)) #G1S

# G1S: CDCA7L, 
fn8.addLines('CDC6', 'red', lty=2)
fn8.addLines("CASP8AP2", 'red', lty=3)
fn8.addLines("DTL", 'red', lty=4) #OK
fn8.addLines("POLD3", 'red', lty=5)

#S RRM2
fn8.addLines('RRM2', 'blue', lty=1)


# G2M: AURKB, CCNA2, NDC80, UBE2C, CDC20, PLK1, ODF2, 
fn8.addLines('CCNB1', 'blue', lty=1)
fn8.addLines('CCNA2', 'blue', lty=2)
fn8.addLines('CDC20', 'blue', lty=3)
fn8.addLines('PLK1', 'blue', lty=4)
fn8.addLines('ODF2', 'blue', lty=5)


# M: CDC25C, BUB1, CCNB2, CENPF, DLGAP5, DEPDC1, HMMR, MKI67, NUSAP1, TPX2, 
fn8.addLines('CDC25C', 'black', lty=1)
fn8.addLines('MKI67', 'black', lty=2)
fn8.addLines('TPX2', 'black', lty=3)
#fn8.addLines('CDC25A', 'black', lty=1)
#


#####################
# pick a gene for each phase
pdf("04_RNA_curve_5genes-2.pdf", width=5,height=4.5)
fn8.addLines=draw_psudoTime_RNA(cellOrder.BC, "BC", gene='CCNE2', 'red', to1=F, ylim=c(0,13)) #G1S
fn8.addLines('DTL', 'red', lty=2)

fn8.addLines('RRM2', 'blue')
fn8.addLines('PLK1', 'blue', lty=2)

fn8.addLines('CDC25C', 'black')
fn8.addLines('GAPDH', 'black', lty=2)

legend(x=par()$usr[2]/5, y=par()$usr[4]+2, lwd=1,
       
       legend=c("CCNE2", 'DTL', "RRM2", "PLK1","CDC25C",'GAPDH'), 
       lty=c(1,2,1,2,1,2),
       col=c('red','red','blue','blue','black','black'),
       ncol=4,
       
       x.intersp=0.5, 
       text.width=18, 
       cex=0.7, 
       
       xpd=T, bty='n')
dev.off()
















#############
# get and show correlated genes
#############
#similar to gene target gene xx
getCoTranscriptedGenes=function(targetGene, dt, method='spearman', p.cutoff=0.05){
  #targetGene="GAPDH"
  #dt=rnaM.log2cpm[, cid.BC]
  # remove all 0 rows
  pP('targetGene', targetGene, '; method',method)
  keep=apply(dt, 1, sum)>0
  dt=dt[keep,]
  
  # order by RNA expression
  targetGene.exp=dt[targetGene,]
  dt.cor=as.data.frame(t( apply(dt, 1, function(x){
    rs=cor.test(x, targetGene.exp, method=method);
    return( c( as.numeric(rs$estimate), rs$p.value) )
  }) ) )
  colnames(dt.cor)=c('cor','p')
  #dt.cor$gene=rownames(dt.cor)
  #order by cor
  dt.cor=dt.cor[order(-abs(dt.cor$cor)),]
  # fileter by p
  dt.cor$adjp=p.adjust(dt.cor$p, method='fdr')
  #dt.cor=dt.cor[which(dt.cor$adjp<p.cutoff),]
  print(dim(dt.cor))
  #
  dt.cor
}




########################
# in BC, cor with CCNB1
coExp.CCNB1=getCoTranscriptedGenes("CCNB1", rnaM.log2cpm[, cid.BC], method="pearson")
head(coExp.CCNB1)


########## for mean+3*sd
pdf("05_RNA_cor_with_CCNB1.hist.pdf", width=4.5,height=3.5)
# hist
hist(coExp.CCNB1$cor, n=50, main="Histgram of correlation with CCNB1", xlab="Spearman correlation")
thresh=mean(coExp.CCNB1$cor)+sd(coExp.CCNB1$cor)*3; thresh
abline(v=thresh, lty=2, col="red")
#
geneDF.coExp=coExp.CCNB1[which(coExp.CCNB1$cor>thresh),]
dim(geneDF.coExp)
text(x=0.3, y=500, labels = paste0("mean+3*sd\nn=",nrow(geneDF.coExp)), adj=0, col='red')
dev.off()

#save
thresh=mean(coExp.CCNB1$cor)+sd(coExp.CCNB1$cor)*1.96; thresh
x=coExp.CCNB1; genes=rownames(x[x$adjp<0.05 & abs(x$cor)>= thresh,])
print(length(genes))
writeLines( genes, "05_RNA_cor_with_CCNB1.gene.txt" )

# draw co exp genes
pdf("05_RNA_cor_with_CCNB1.linePlot.pdf", width=5,height=4.2)
fn8.addLines=draw_psudoTime_RNA(cellOrder.BC, "BC", gene='CCNB1', #to1=T,ylim=c(0,1,1),
                                color1='red',rawColor='white')
for(i in 2:length(genes)){
  gene=genes[i]
  fn8.addLines(gene, '#aaaaaa')
}
fn8.addLines("CCNB1", 'red', lwd=2)
title("CCNB1 Co-expressed genes")
dev.off()



########## for mean+4*sd
pdf("05_RNA_cor_with_CCNB1-4sd.hist.pdf", width=4.5,height=3.5)
# hist
hist(coExp.CCNB1$cor, n=50, main="Histgram of correlation with CCNB1", xlab="Spearman correlation")
thresh=mean(coExp.CCNB1$cor)+sd(coExp.CCNB1$cor)*4; thresh
abline(v=thresh, lty=2, col="red")
#
geneDF.coExp=coExp.CCNB1[which(coExp.CCNB1$cor>thresh),]
dim(geneDF.coExp)
text(x=0.4, y=500, labels = paste0("mean+4*sd=",round(thresh,2),
                                   "\nn=",nrow(geneDF.coExp)), adj=0, col='red')
dev.off()

#save
x=coExp.CCNB1; genes=rownames(x[x$adjp<0.05 & x$cor>= thresh,])
print(length(genes))
writeLines( genes, "05_RNA_cor_with_CCNB1-4sd.gene.txt" )

# draw co exp genes
pdf("05_RNA_cor_with_CCNB1-4sd.linePlot.pdf", width=5,height=4.2)
fn8.addLines=draw_psudoTime_RNA(cellOrder.BC, "BC", gene='CCNB1', #to1=T,ylim=c(0,1,1),
                                color1='red',rawColor='white')
for(i in 2:length(genes)){
  gene=genes[i]
  fn8.addLines(gene, '#aaaaaa')
}
fn8.addLines("CCNB1", 'red', lwd=2)
title("CCNB1 Co-expressed genes")
dev.off()




########## for 0.4
pdf("05_RNA_cor_with_CCNB1-0.4.hist.pdf", width=4.5,height=3.5)
# hist
hist(coExp.CCNB1$cor, n=50, main="Histgram of correlation with CCNB1", xlab="Spearman correlation")
thresh=0.4 #mean(coExp.CCNB1$cor)+sd(coExp.CCNB1$cor)*4; thresh
abline(v=thresh, lty=2, col="red")
#
geneDF.coExp=coExp.CCNB1[which(coExp.CCNB1$cor>thresh),]
dim(geneDF.coExp)
text(x=0.45, y=500, labels = paste0("thresh=",round(thresh,2),
                                   "\nn=",nrow(geneDF.coExp)), adj=0, col='red')
dev.off()

#save
x=coExp.CCNB1; genes=rownames(x[x$adjp<0.05 & x$cor>= thresh,])
print(length(genes))
writeLines( genes, "05_RNA_cor_with_CCNB1-0.4.gene.txt" )

# draw co exp genes
pdf("05_RNA_cor_with_CCNB1-0.4.linePlot.pdf", width=4.3,height=4.3)
fn8.addLines=draw_psudoTime_RNA(cellOrder.BC, "BC", gene='CCNB1', #to1=T,ylim=c(0,1,1),
                                color1='red',rawColor='white')
for(i in 2:length(genes)){
  gene=genes[i]
  fn8.addLines(gene, '#aaaaaa')
}
fn8.addLines("CCNB1", 'red', lwd=2)
title("CCNB1 Co-expressed genes")
dev.off()
















########################
# in BC, cor with CCNB1
# draw co exp genes
rnaM.log2cpm.BC=rnaM.log2cpm[, as.character(cid.BC)]


drawCoExp=function(gName, method2="pearson", outdir="co-exp/", thresh0=0.4){
  coExp_df=getCoTranscriptedGenes(gName, rnaM.log2cpm.BC, method=method2)
  #head(coExp_df)
  
  pdf( paste0(outdir, "RNA_coExp_",gName,"-gt0.4.linePlot.",method2,".pdf"), width=4.3,height=4.3)
  thresh1=mean(coExp_df$cor)+1.96*sd(coExp_df$cor);thresh1
  thresh=thresh1
  if(thresh<thresh0)thresh=thresh0;
  
  # fig 1
  hist(coExp_df$cor, n=100, main=paste0(gName, "\nthresh Black=", round(thresh1,2), 
                                        ", Red=", round(thresh, 2)))
  abline(v=thresh1, lty=2, col="black")
  abline(v=thresh, lty=2, col="red")
  
  
  x=coExp_df; genes=rownames(x[x$adjp<0.05 & x$cor>= thresh,])
  #x=coExp_df; genes=rownames(x[x$p<0.05 & x$cor>= thresh,])
  cat("gene length=",length(genes), '\n')
  writeLines(genes, paste0(outdir, "RNA_coExp_",gName,"-gt0.4.linePlot.",method2,".gene.txt"))
  #
  # fig 2
  fn8.addLines=draw_psudoTime_RNA(cellOrder.BC, "BC", gene=gName, #to1=T,ylim=c(0,1,1),
                                  color1='red',rawColor='white')
  if(length(genes)>2){
    for(i in 2:length(genes)){
      gene=genes[i]
      fn8.addLines(gene, '#aaaaaa')
    }
  }
  fn8.addLines(gName, 'red', lwd=2)
  title( paste0(gName," Co-expressed genes \n Method:", method2, ", n:", length(genes)) )
  dev.off()
}
#
drawCoExp("CCNB1", "pearson")
#drawCoExp("CCNB1", "pearson")
drawCoExp("CCNB2", "spearman")
#drawCoExp("CCNB1", "spearman")
drawCoExp("CCNE2", "spearman",thresh0=0)
drawCoExp("CCNA2", "spearman",thresh0=0.35)

#
drawCoExp("CDK4", "spearman",thresh0=0.15)
drawCoExp("CDK6", "spearman",thresh0=0.15)
drawCoExp("MKI67", "spearman",thresh0=0.4)
drawCoExp("CDC25C", "spearman",thresh0=0.35)


##
dim(rnaM.log2cpm.BC)
rnaM.log2cpm.BC[1:4,1:4]

coExp_df=as.data.frame(t(apply(rnaM.log2cpm.BC,1, function(x){
  rs=cor.test( as.numeric(rnaM.log2cpm.BC['CCNE2',]), 
               as.numeric(x), method="spearman" )
  c( rownames(x), as.numeric(rs$estimate), rs$p.value)
})))
colnames(coExp_df)=c("cor","p")
coExp_df$p.adj=p.adjust(coExp_df$p, method="fdr")
dim(coExp_df)
head(coExp_df)
#                 cor         p     p.adj
#A1BG-AS1 -0.02783074 0.7227070 0.9827703

hist(coExp_df$p, n=100)
hist(coExp_df$p.adj, n=100)
table(coExp_df$p<0.05)
table(coExp_df$p.adj<0.05) #4
coExp_df[which(coExp_df$p.adj<0.05),]
#











############
# all correlation in BC
dim(rnaM.log2cpm.BC) # 18662   165
rnaM.log2cpm.BC[1:4,1:4]

# filter out all 0 rows in BC
rnaM2=rnaM.log2cpm.BC[rowSums(rnaM.log2cpm.BC)>0,]
dim(rnaM2) # 18240   165
rnaM2=rnaM2[grep("-",rownames(rnaM2),invert=T),]
dim(rnaM2) #17414   165
rnaM2[1:4,1:4]
write.table(rnaM2, "co-exp/RNA_log2cpm_BC.165.df.txt")


## select 1e6 pair of genes
uniqGenes=rownames(rnaM2)
length(uniqGenes)
set.seed(202106)
cor_df0=data.frame(
  g1=sample(uniqGenes, 1e6, replace = T),
  g2=sample(uniqGenes, 1e6, replace = T)
)
head(cor_df0)
dim(cor_df0) #1e4 2

################
# rm duplicate 1
cor_df0 = unique(cor_df0)
dim(cor_df0) #99987     3
# rm duplicate 2
#cor_df0$key=paste(cor_df0$g1,cor_df0$g2, sep="|")
#t1=c(paste(cor_df0$g1,cor_df0$g2, sep="|"), paste(cor_df0$g2,cor_df0$g1, sep="|") )
#head(t1)
#
#df1=c()
#df2=c()

cor_df1=t(apply(cor_df0, 1, function(x){
  tmp=data.frame(
    X=rnaM2[x[1],],
    Y=rnaM2[x[2],]
  )
  tmp=tmp[rowSums(tmp)>0,]
  if(nrow(tmp)>10){
    rs=cor.test(tmp$X, tmp$Y)
    c( as.numeric(rs$estimate), rs$p.value, nrow(tmp))
  }else{
    c(-1, -1, -1)
  }
}))
#1e5: 10:23--> 10:24
#2e5: 10:26--> 10:29
#2e5: 10:34--> 10:37
#2e6: 14:07-->
colnames(cor_df1)=c("R","p","N")
tail(cor_df1)
cor_df2=cbind(cor_df0[,1:2], cor_df1)
head(cor_df2)
# filter out -1 and 1
cor_df2=cor_df2[which(cor_df2$N>30),]
cor_df2$p.adj=p.adjust(cor_df2$p, method="fdr")
dim(cor_df2) #84045     5
head(cor_df2)

pdf("co-exp/BC_gene_gene_cor.spearman.pdf", width=4,height=4)
hist(as.numeric(cor_df2$N), n=1000)
hist(as.numeric(cor_df2$p), n=1000)
hist(as.numeric(cor_df2$p.adj), n=1000)

hist(as.numeric(cor_df2$R), n=1000)
hist(as.numeric(cor_df2$R), n=40)
# fig 1
plot( density(as.numeric(cor_df2$R), n=50), col="blue", mgp=c(2,1,0),
     xlab="Spearman correlation", main="MDA-MB-468 \ngene-gene expression correlation")
hist(as.numeric(cor_df2$R), n=50, freq = F, add=T, border ="#666666", col=NA)
abline(v=0, lty=2, lwd=2, col="red")
# fig 2
hist(as.numeric(cor_df2[which(cor_df2$p.adj<0.05),]$R), n=50,
     xlab="Spearman correlation\n(p.adj<0.05)",
     main="MDA-MB-468 \ngene-gene expression correlation")
abline(v=0, lty=2, lwd=2, col="red")
dev.off()
#cor_df2[which(cor_df2$R>0.5),]
#

#
plot(rnaM2['MAN1C1',], rnaM2['LINC01551',], cex=0.2)
#












############
# all correlation in HeLa
dim(rnaM)
rnaM.H=rnaM[,cid.HeLa]
rnaM.H=rnaM.H[rowSums(rnaM.H)>0,] # filter all 0 rows

rnaM3=apply(rnaM.H, 2, function(x){
  log2(x/sum(x)*1e6+1)
})
dim(rnaM3) #16989    57
rnaM3[1:4,1:4]

write.table(rnaM3, "co-exp/RNA_log2cpm_HeLa.57.df.txt")



getCor_DF=function(rnaM0){
  ## select 1e6 pair of genes
  uniqGenes=rownames(rnaM0)
  #length(uniqGenes)
  set.seed(202106)
  cor_df0=data.frame(
    g1=sample(uniqGenes, 1e6, replace = T),
    g2=sample(uniqGenes, 1e6, replace = T)
  )
  # rm duplicate 1
  cor_df0 = unique(cor_df0)
  cat("rm dup 1: ", dim(cor_df0), "\n" )
  # rm dup 2
  cor_df0=cor_df0 %>% rowwise() %>% 
    mutate(class = paste0(sort(c(g1, g2)), collapse = '')) %>% 
    group_by(class) %>% 
    slice_sample() # %>% data.frame %>% select( -c(class) )
  cor_df0=data.frame(cor_df0)
  cor_df0=cor_df0[,c("g1","g2")]
  cat("rm dup 2: ", dim(cor_df0), "\n" ) #99987     3
  
  cor_df1=t(apply(cor_df0, 1, function(x){
    tmp=data.frame(
      X=rnaM0[x[1],],
      Y=rnaM0[x[2],]
    )
    tmp=tmp[rowSums(tmp)>0,]
    if(nrow(tmp)>10){
      rs=cor.test(tmp$X, tmp$Y)
      c( as.numeric(rs$estimate), rs$p.value, nrow(tmp))
    }else{
      c(-1, -1, -1)
    }
  }))
  
  colnames(cor_df1)=c("R","p","N")
  cor_df2=cbind(cor_df0[,1:2], cor_df1)
  # filter out -1 and 1
  cor_df2=cor_df2[which(cor_df2$N>30),]
  cor_df2$p.adj=p.adjust(cor_df2$p, method="fdr")
  cor_df2
}
cor_df2.HeLa=getCor_DF(rnaM3) # time consuming

dim(cor_df2.HeLa)
head(cor_df2.HeLa)

pdf("co-exp/HeLa_gene_gene_cor.spearman.pdf", width=4,height=4)
# fig 1
plot( density(as.numeric(cor_df2.HeLa$R), n=50), col="blue", mgp=c(2,1,0),
      xlab="Spearman correlation", main="HeLa \ngene-gene expression correlation")
hist(as.numeric(cor_df2.HeLa$R), n=50, freq = F, add=T, border ="#666666", col=NA)
abline(v=0, lty=2, lwd=2, col="red")
# fig 2
hist(as.numeric(cor_df2.HeLa[which(cor_df2.HeLa$p.adj<0.05),]$R), n=50,
     xlab="Spearman correlation\n(p.adj<0.05)",
     main="HeLa \ngene-gene expression correlation")
abline(v=0, lty=2, lwd=2, col="red")
dev.off()
