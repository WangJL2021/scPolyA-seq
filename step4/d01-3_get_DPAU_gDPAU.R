# get DPAU and gDPAU

setwd('/home/wangjl/data/apa/20200701Fig/f4/DPAU/')
getwd()

library(ggplot2)
library(Cairo)
library(pheatmap)

printP=function(...){ print(paste(...)) }
printP0=function(...){ print(paste0(...)) }

#1. tools
pP=function(...){  print(paste(...)) }
pP0=function(...){  print(paste0(...)) }
ph=function(...){  print(head(...)) }
pd=function(...){  print(dim(...)) }
pl=function(...){  print(length(...)) }
debug=function(...){ pd(...); ph(...) }

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

#################
#load data
#################
#(1) cell info
cellInfo=read.table("/data/jinwf/wangjl/apa/20200701Fig/f3/cell_cycle/cellInfo.V6.txt")
dim(cellInfo) #222   8
head(cellInfo)
#              cid cellType0 cluster cellType cellType2 countsPerCell geneNumber cellCycle
#c12ROW03 c12ROW03   unknown       0     BC_0      BC_0       1820666       6659       MG1
#c12ROW04 c12ROW04   unknown       0     BC_0      BC_0       2816130       7930       MG1
table(cellInfo$cellType)
#BC_0        BC_1 HeLa_normal   HeLa_sync 
#92          73          30          27

#(2) load apa site, and apa Matrix
apaSiteFile="/data/jinwf/wangjl/apa/190705PAS/bed/pasPostions_Location_transcriptName-noChrM_noInnerPrime_PY-filterCountCell-motif-Header.bed"
apaMatrixFile="/data/jinwf/wangjl/apa/190705PAS/bed/apa_Matrix_afterFilterByCounts_cell_20222.csv"
#
apaSite=read.table(apaSiteFile, header = T,row.names = 1, stringsAsFactors = F)
row.names(apaSite)=apaSite$PASid
dim(apaSite) #20222 12
apaSite[1:4,]
#                          PASid   chr     pos strand count region transcript                down20   gene
# chr10:320403:-   chr10:320403:- chr10  320403      - 10945   UTR3  DIP2C-201 TGAAAATAGCAGTTTCTTAAT  DIP2
# apa matrix
apaM=read.csv(apaMatrixFile, header = T,row.names = 1, stringsAsFactors = F)
apaM=apaM[, rownames(cellInfo)]
dim(apaM) #20222 222
head(apaM[, 1:5])
#                c12ROW03 c12ROW04 c12ROW05 c12ROW06 c12ROW08
#chr10:320403:-         1      362      358       16        0
#chr10:855485:-        36        1        4        0        0


#################
# get DPAU
#################
#1. get genes with nulti apa sites
gene.apaSiteNumber=(function(){
  rs=getDF_fromNamesXX(table(apaSite$gene))
  colnames(rs)=c('gene', 'apaSite');
  rs=rs[which(rs$apaSite>=2 & rs$gene!='.'),];  print(dim(rs))
  rs
})()
tail(gene.apaSiteNumber)
head(gene.apaSiteNumber) #4103    2
#           gene apaSite
#MALAT1   MALAT1      46
#EGFR       EGFR      40
gene.2sites=rownames( gene.apaSiteNumber[which(gene.apaSiteNumber$apaSite==2),] )
length( gene.2sites ) #2082


#2. get DPAU
DPAU=(function(){
  apa_df_2sites=NULL;
  for(i in 1:length(gene.2sites) ){
    gene=as.character( gene.2sites[i] )
    if(i>10){
      #break;
    }
    if(i%%200==0){
      print(paste('processing ',i) )
    }
    
    #(1) one gene, get 2 apa site
    tmpDF=apaSite[which(apaSite$gene == gene),]
    # order by position
    tmpDF=tmpDF[order(tmpDF$pos),]
    
    #get 2 sites
    apa_2sites=as.character( tmpDF$PASid )
    
    # get the most 3' end apa site;
    if( as.character( tmpDF[1, 'strand'] )=='+'){
      apa_3end=apa_2sites[length(apa_2sites)]
    }else{
      apa_3end=apa_2sites[1]
    }
    
    #(2) get apa matrix for the gene;
    tmpDF2=apaM[tmpDF$PASid,]
    
    #get one row template;
    tmpDF0=apaM[1,];
    row.names(tmpDF0)=gene
    #
    for(j in 1:ncol(tmpDF0) ){
      cell=colnames(tmpDF0)[j]
      
      # filter out sites with less than 5 reads at both sites. #-----------> regard a gene's sum counts less than 5 in a cell as noise.
      if( sum(tmpDF2[,cell]) < 5 ){
        tmpDF0[1, cell]=NA;
      }else{
        # calculate most 3'end read pct;
        DPAU=tmpDF2[apa_3end,cell]/sum(tmpDF2[,cell])*100
        tmpDF0[1, cell]=DPAU;
      }
    }
    # output;
    apa_df_2sites=rbind(apa_df_2sites, tmpDF0)
  } #2min;
  return(apa_df_2sites)
})()

# check
dim(DPAU) #2082  222
DPAU[1:5,1:10]
#         c12ROW03   c12ROW04 c12ROW05 c12ROW06 c12ROW08 c12ROW15 c12ROW19 c12ROW23 c12ROW25  c12ROW27
#AARD           NA   0.000000       NA       NA       NA       NA       NA       NA      100        NA
#AARS            0 100.000000       NA       NA       NA       80      100      100       NA        NA
# save to file
write.table(DPAU, '0_DPAU_2082x222.df.txt')


library(vioplot)
CairoPDF('1_DPAU_QC.pdf',width=12,height=8)
par(mfrow=c(2,3))
# non NA per gene
hist( apply(DPAU, 1, function(x){
  sum(!is.na(x))
}), n=100, xlab="Obs. number per gene")
boxplot(t( DPAU[1:20, ]), main="DPAU per gene", las=2)
vioplot(t( DPAU[1:20, ]), main="DPAU per gene", las=2)
# non NA per cell
hist( apply(DPAU, 2, function(x){
  sum(!is.na(x))
}), n=50, xlab="Obs. number per cell")
boxplot(DPAU[,1:20], main="DPAU per cell", las=2)
vioplot(DPAU[,1:20], main="DPAU per cell", las=2)
dev.off()


(function(){
  rs=apply(DPAU, 1, function(x){
    sum(!is.na(x))
  })
  rs=rs[rs<50]
  writeLines(names(rs), '1_obsLessThan50.gene.txt')
  length(rs) #611
})()






##############
# plot gene DPAU distribution by cell cycle
##############
head(cellInfo)
# plot gene DPAU distribution by cell cycle
drawDPAU_distribution=function(DPAU_, ymax=0.07, main="DPAU distribution of genes within each cell",
                               sub="[One line one cell]", showLegend=T){
  # set colors for cell cycles
  col=RColorBrewer::brewer.pal(n = 5,name = "Set2")
  #
  plot(density(DPAU_[which( !is.na(DPAU_[,1]) ),1]), col="white",
       xlab="DPAU", ylim=c(0, ymax), main=main,sub=sub)
  for(i in seq(1,ncol(DPAU_))){
    cid=colnames(DPAU_)[i]
    cellCycle=cellInfo[cid, 'cellCycle']
    color=switch(cellCycle, G1S=col[1], S=col[2], G2M=col[3], M=col[4], MG1=col[5])
    lines(density(DPAU_[which( !is.na(DPAU_[,i]) ),i]), col=color)
  }
  #
  if(showLegend){
    legend("topleft", legend=c('G1S', 'S', 'G2M', 'M', 'MG1'), col =col, lty = 1,
           xjust = 0, yjust = 1, x.intersp = 1, y.intersp =1, adj = c(0, 0.5), text.width = NULL,xpd=TRUE)
  }
}

###
CairoPDF('1_DPAU_distribution_each_gene_By_cellCycle.pdf',width=24,height=4)
par(mfrow=c(1,6))
## use 0 and 100%
drawDPAU_distribution(DPAU)
# by each cell cycle
for(cc in c('G1S', 'S', 'G2M', 'M', 'MG1')){
  cids=rownames(cellInfo[which(cellInfo$cellCycle==cc),])
  drawDPAU_distribution(DPAU[, cids], main=paste0("DPAU distribution of genes\nphase:", cc, ', cellNumber:',length(cids)),
                        showLegend=F)
}
dev.off()


# find the outlier
outlier=(function(){
  DPAU_=DPAU
  df=NULL;
  for(i in seq(1,ncol(DPAU_))){
    #if(i>2) break;
    cid=colnames(DPAU_)[i]
    #print(cid) #"c12ROW03"
    cellCycle=cellInfo[cid, 'cellCycle']
    dst=density(DPAU_[which( !is.na(DPAU_[,i]) ),i]);
    #
    y0=c(); y100=c();
    n=length(dst$x)
    #print(n) #512
    for(j in 1:n){
      x=dst$x[j]; y=dst$y[j]
      #print(j)
      if(x>-0.5 & x<0.5){ y0=c(y0, y);}
      if(x>99.5 & x<100.5){ y100=c(y100, y);}
    }
    tmp=data.frame(
      cid=cid,
      y0=mean(y0),
      y100=mean(y100)
    )
    df=rbind(df, tmp)
  }
  df
})()
head(outlier)
hist(outlier$y0, n=200)
hist(outlier$y100, n=200)
#
outlier[which(outlier$y0>0.014),] #73 c19ROW01
outlier[which(outlier$y100>0.03),] 
#        cid         y0       y100
#73 c19ROW01 0.02084412 0.06828576
#74 c19ROW03 0.01184019 0.04188290
#89 c19ROW36 0.01003147 0.03210570
cid.outlier=as.character(outlier[which(outlier$y100>0.03),'cid'] ); cid.outlier #"c19ROW01" "c19ROW03" "c19ROW36"
cellInfo[cid.outlier, ] #all BC_0
#             cid cellType0 cluster cellType cellType2 countsPerCell geneNumber cellCycle
#c19ROW01 c19ROW01   unknown       0     BC_0      BC_0        226252       3483       G2M
#c19ROW03 c19ROW03   unknown       0     BC_0      BC_0        199351       3945         S
#c19ROW36 c19ROW36   unknown       0     BC_0      BC_0        566442       5519       MG1
head(cellInfo[order(cellInfo$geneNumber),]) #c19ROW01, least gene;
head(cellInfo[order(cellInfo$countsPerCell),]) #c19ROW03, least counts;
(function(){
  rs=getDF_fromNamesXX( apply(DPAU, 2, function(x){sum(!is.na(x))}) )
  colnames(rs)=c('cid', 'obs')
  rs=rs[order(rs$obs),]
  print( grep('c19ROW36', rownames(rs)) ) #28位
  head(rs, n=10)
})()











##############
# plot gene DPAU distribution by cell type
##############
drawDPAU_distribution_ByCellType=function(DPAU_, ymax=0.07, main="DPAU distribution of genes", 
                                          sub="[One line one cell]", showLegend=T){
  # set colors for cell cycles
  col=c("#F8766D", "#00BFC4")
  #
  plot(density(DPAU_[which( !is.na(DPAU_[,1]) ),1]), col="white",
       xlab="DPAU", ylim=c(0, ymax), main=main,sub=sub)
  for(i in seq(1,ncol(DPAU_))){
    cid=colnames(DPAU_)[i]
    cellType=substr( cellInfo[cid, 'cellType'], 0, 2)
    color=switch(cellType, BC=col[1], He=col[2])
    
    lines(density(DPAU_[which( !is.na(DPAU_[,i]) ),i]), col=color)
  }
  #
  if(showLegend){
    legend("topleft", legend=c('BC', 'HeLa'), col =col, lty = 1,
           xjust = 0, yjust = 1, x.intersp = 1, y.intersp =1, adj = c(0, 0.5), text.width = NULL,xpd=TRUE)
  }
}

###
CairoPDF('1_DPAU_distribution_each_gene_By_cellType.pdf',width=12,height=8)
par(mfrow=c(2,3))
drawDPAU_distribution_ByCellType(DPAU)

for(ct in c("BC","He")){
  cids=rownames(cellInfo[which( substr( cellInfo$cellType,1,2 )==ct),])
  drawDPAU_distribution_ByCellType(DPAU[, cids], 
                                   main= paste0("DPAU distribution of genes\n",ct,", cellNumber:", length(cids)), showLegend = F )
}

for(ct in c("HeLa_normal",  "HeLa_sync", "HeLa_syncMix")){
  cids=rownames(cellInfo[which( cellInfo$cellType2 ==ct),])
  drawDPAU_distribution_ByCellType(DPAU[, cids], 
                                   main= paste0("DPAU distribution of genes\n",ct,", cellNumber:", length(cids)), showLegend = F )
}
dev.off()








#################
# a gene a line
#################
drawDPAU_distribution_byGene=function(DPAU_, geneSymbol,ymax=18, main="DPAU distribution of genes", 
                                      sub="[One line one gene]", showLegend=F){
  
  dataRow1=DPAU_[geneSymbol, which( !is.na(DPAU_[i,]) )];
  
  if( ncol(dataRow1)<1 ){
    printP('All NA, gene=',geneSymbol, '; ncol=',ncol(dataRow1))
    #print(dim(dataRow1) )
    return();
  }
  
  #
  par(mfrow=c(1,3))
  hist( as.numeric(dataRow1), col='red', main=paste("gene:",geneSymbol),xlab="", sub=sub, n=225, border=NA)
  plot(density( as.numeric(dataRow1)  ), col='red', main=paste("gene:",geneSymbol),
       xlim=c(-10,110),ylim=c(0,0.2), sub=sub)
  plot(density( as.numeric(dataRow1)  ), col='red', main=paste("gene:",geneSymbol),
       xlim=c(-10,110), sub=sub)
}
# test
i=51; drawDPAU_distribution_byGene(DPAU, as.character( rownames(DPAU)[i]), sub=i)
i=52; drawDPAU_distribution_byGene(DPAU, as.character( rownames(DPAU)[i]), sub=i)
i=53; drawDPAU_distribution_byGene(DPAU, as.character( rownames(DPAU)[i]), sub=i)
pdf('2-densityOfGeneDPAU_1m.pdf', width=20,height=4)
for(i in 1:nrow(DPAU)){
  if(i>200){break;}
  drawDPAU_distribution_byGene(DPAU, as.character( rownames(DPAU)[i]), sub=i)
}
dev.off()



# all genes
color=rgb(255, 0, 0, 25, maxColorValue=255)
drawDPAU_distribution2=function(DPAU_, ymax=18, main="DPAU distribution of genes", 
                                sub="[One line one gene]", showLegend=F){
  # set colors for cell cycles
  #col=c("#F8766D", "#00BFC4")
  #
  dataRow1=DPAU_[1, which( !is.na(DPAU_[1,]) )];
  plot(density( as.numeric(dataRow1) ) , col="white",
       xlab="DPAU", ylim=c(0, ymax), main=main, sub=sub)
  for(i in seq(1,nrow(DPAU_))){
    genename=rownames(DPAU_)[i]
    
    dataRow1=DPAU_[i, which( !is.na(DPAU_[i,]) )];
    
    if(length(dataRow1)<2){
      print( paste('Discard(length<2): ',i,genename) )
      next;
    }
    
    lines(density( as.numeric(dataRow1)  ), col=color)
  }
  #
  if(showLegend){
    legend("topleft", legend=c('BC', 'HeLa'), col =col, lty = 1,
           xjust = 0, yjust = 1, x.intersp = 1, y.intersp =1, adj = c(0, 0.5), text.width = NULL,xpd=TRUE)
  }
  print(i)
}
CairoPDF('2-densityOf2sites_GeneDPAU_all.pdf', width=10,height=5)
par(mfrow=c(1,2))
drawDPAU_distribution2(DPAU, ymax=0.2)
drawDPAU_distribution2(DPAU, ymax=15)
dev.off()






#############
##get each gene's non-NA cell number, mean, sd of DPAU
#############
geneAPAInfo=(function(){
  rsDF=NULL
  for(i in 1:length(gene.2sites)){
    #if(i>10){break;}
    geneSymbol=gene.2sites[i]
    #printP(i, ' ', geneSymbol)
    dataRow1=DPAU[geneSymbol, which( !is.na(DPAU[geneSymbol,]) )];
    tmp=as.numeric(dataRow1)
    
    if(length(tmp)<2){
      print(paste('Discard(obs<2): ', i, ' ', geneSymbol) )
      next;
    }
    
    tmpDF=data.frame(
      gene=geneSymbol,
      DPAUcellNumber=length(tmp),
      meanDPAU=mean(tmp),
      sdDPAU=sd(tmp)
    )
    rsDF=rbind(rsDF, tmpDF)
  } #10seconds
  #
  rownames(rsDF)=rsDF$gene
  print(dim(rsDF)) #2082 4
  print( head(rsDF) )
  #    gene cellNumber meanDPAU     sdDPAU
  #1   RPL8        225 99.86553  0.3462266
  #2   RPL3        225 99.73653  0.9470744
  #3  RPS23         19 82.75423 36.9176673
  #View(rsDF)
  #
  plot(rsDF$mean, rsDF$sd, xlab="Mean of DPAU", ylab="sd of DPAU", main="[One gene one point]" )
  rsDF
})()
plot(geneAPAInfo$mean, geneAPAInfo$sd, xlab="Mean of DPAU", ylab="sd of DPAU", main="[One gene one point]" )
#             gene DPAUcellNumber meanDPAU   sdDPAU
#AARD         AARD             27 53.21499 48.17641
#AARS         AARS            113 44.73042 43.32136





# add RNA read counts for each genes
rnaM=read.csv('/data/jinwf/wangjl/apa/20200701Fig/f2/BC_HeLa/BC_HeLa.222cells.count.V4.csv', row.names = 1)
dim(rnaM) #18662 225
rnaM[1:10,1:10]
#
geneAPAInfo2=(function(){
  rnaM.cpm=apply(rnaM, 2, function(x){x/sum(x)*1e6})
  countsPerGene=getDF_fromNamesXX( apply(rnaM.cpm, 1, sum)  )
  pd(countsPerGene)
  #
  gene.common=intersect(rownames(geneAPAInfo), rownames(rnaM.cpm))
  countsPerGene=countsPerGene[gene.common, ]
  geneAPAInfo=geneAPAInfo[gene.common,]
  #
  pd(countsPerGene) #2036 2
  ph(countsPerGene)
  
  geneAPAInfo$cpm=countsPerGene$value
  geneAPAInfo$log2cpm=log2( 1 + countsPerGene$value/ncol(rnaM) )
  head(geneAPAInfo)
  #        gene DPAUcellNumber meanDPAU     sdDPAU RNACounts
  #RPL8     RPL8        225 99.86553  0.3462266    712259
  #RPL3     RPL3        225 99.73653  0.9470744    542531
  geneAPAInfo
})()
pd(geneAPAInfo) #2082 4
pd(geneAPAInfo2) #2036 6
head(geneAPAInfo2)
#             gene DPAUcellNumber meanDPAU   sdDPAU       cpm  log2cpm
#AARD         AARD             27 53.21499 48.17641 14908.435 13.86394
#AARS         AARS            113 44.73042 43.32136 16194.540 13.98331
#save to file
write.table(geneAPAInfo2, '3_gene_DPAUcellNumber_meanDPAU_sdDPAU_cpm_log2cpm.df.txt')


(function(){
  rsDF=geneAPAInfo2
  
  CairoPDF("3_1_gene_DPAUcellNumber_meanDPAU_sdDPAU_RNACounts.pdf", width=4.2, height=3.2)
  g1=ggplot(rsDF, aes(meanDPAU, sdDPAU, color=DPAUcellNumber))+
    geom_point(size=0.4) + 
    scale_colour_gradientn('Cell number', colours=c('blue', 'white', 'red') )+
    labs(x="Mean of gene' DPAU", y="Standard deviation of gene' DPAU")+
    theme_bw()
  print(g1)
  ##小结: 
  #
  g2=ggplot(rsDF, aes(log2cpm, sdDPAU, color=DPAUcellNumber))+
    geom_point(size=0.4) + 
    scale_colour_gradientn('Cell number', colours=c('blue', 'white', 'red') )+
    labs(x="RNA expression level [log2(cpm+1)]", y="Standard deviation of gene' DPAU")+
    geom_vline(aes(xintercept=6.5), colour="#000000", linetype="dashed")+
    theme_bw();#g2
  print(g2)
  #table(rsDF$log2cpm>6)
  #
  g3=ggplot(rsDF, aes(DPAUcellNumber, log2cpm, color=meanDPAU))+
    geom_point(size=0.4) + 
    scale_colour_gradientn('Mean DPAU', colours=c('blue', 'white', 'red') )+
    labs(x="Observation number",y="RNA expression level [log2(cpm+1)]")+
    theme_bw()
  print(g3)
  #
  g4=ggplot(rsDF, aes(meanDPAU,log2cpm, color=DPAUcellNumber))+
    geom_point(size=0.4) + 
    scale_colour_gradientn('Cell number', colours=c('blue', 'white', 'red') )+
    labs(x="Mean of gene' DPAU",y="RNA expression level [log2(cpm+1)]")+
    theme_bw()
  print(g4)
  
  # hist plot of a gene's DPAU
  plotDPAU_hist_byGene=function(gene, n=0){
    tmp=as.numeric(DPAU[gene, !is.na(DPAU[gene,])])
    print(length(tmp))
    if(n==0)n=length(tmp);
    hist( tmp, n=n, xlab="DPAU", main=paste('gene:',gene) )
  }
  
  ## gene's mean DPAU around 50
  print( rsDF[which(rsDF$mean>49.8 & rsDF$mean<50),] )
  plotDPAU_hist_byGene('FANCA', 50)
  plotDPAU_hist_byGene('NUDT15', 30)
  #hist( as.numeric(DPAU['NCAPD3', !is.na(DPAU['NCAPD3',])]), n=20 )
  #
  ## gene's mean DPAU around 0
  print( rsDF[which(rsDF$mean>0 & rsDF$mean<1),] )
  plotDPAU_hist_byGene('NDUFB10', 30)
  plotDPAU_hist_byGene('PHPT1', 30)
  #
  ## gene's mean DPAU around 100
  print( rsDF[which(rsDF$mean>99.9 & rsDF$mean<100),] )
  plotDPAU_hist_byGene('UBC', 30)
  plotDPAU_hist_byGene('RPL7', 30)
  #
  dev.off()
})()

# save gene list for 2 groups
(function(){
  rsDF=geneAPAInfo2
  head(rsDF)
  dim(rsDF) #2036    6
  #
  geneSet1=rsDF[which(rsDF$log2cpm<=6.5),]
  writeLines( rownames(geneSet1), "3_2_lowRNA_sdDPAU.gene.txt" ) #1515
  geneSet2=rsDF[which(rsDF$log2cpm>6.5),]
  writeLines( rownames(geneSet2), "3_2_highRNA_sdDPAU.gene.txt" ) #521
})()
#
table(rsDF$log2cpm>6.5)
#
















#############################
# get gDPAU
#############################
getGeneralRatio_Long_vs_Short=function(Debug=F){
  gene.set1=unique(apaSite$gene)
  pl(gene.set1) #9455
  
  ###############
  # step1 filter and get gene list with multiple polyA sites.
  gene_df=getDF_fromNamesXX(table(apaSite$gene)) #9455    2
  # remove NA gene name
  colnames(gene_df)=c('gene', 'apaSites')
  gene_df=gene_df[gene_df$gene!='.',] # 9454    2
  #remove 1 site
  gene_df2=gene_df[which(gene_df$apaSites>1),] #4103    2
  
  ###############
  #step2 get gDPAU
  #w=1
  gDPAU=NULL;
  gPPAU=NULL;
  gRatio=NULL;
  for(gene in rownames(gene_df2)){
    #w=w+1
    #if(w>30)break;
    
    #(1) use gene name to get its polyA sites
    polyA.sites=apaSite[which(apaSite$gene==gene),]
    #(2) order the polyA site
    if(polyA.sites[1,'strand']=="+"){
      polyA.sites=polyA.sites[order(polyA.sites$pos),]
    }else{
      polyA.sites=polyA.sites[order(-polyA.sites$pos),]
    }
    #print(polyA.sites[, -c(1,7,8,10)])
    polyA.list=rownames(polyA.sites)
    #(3) get apaM of this gene by polyA.list
    apaM.gene=apaM[polyA.list, ]
    #(4) calculate a gene's gDPAU, gPPAU in each cell, filter genes whose read count at all polyA sites of this gene is less than 5
    # rename as NA, or 0? 0 will be better I think
    gD.tmp=NULL;
    gP.tmp=NULL;
    gR.tmp=NULL;
    for(cid in colnames(apaM.gene)){
      apaCounts.cell=as.numeric(apaM.gene[,cid]);
      gD=0;gP=0;gR=0;
      if(sum(apaCounts.cell)>=5){ #less than 5, gD, gP=NA will not chanage;
        # for genes whose apa reads larger than 5, add one reads for all polyA sites to eliminate 0 division later
        #apaCounts.cell=apaCounts.cell+1
        # counts to pct
        apaCounts.cell.pct=apaCounts.cell/sum(apaCounts.cell);
        n=length(apaCounts.cell.pct);
        for(i in 1:length(apaCounts.cell.pct)){
          gD=gD+ (i-1)*apaCounts.cell.pct[i]
          gP=gP+ (n-i)*apaCounts.cell.pct[i]
        }
        gD=gD/(n-1);
        gP=gP/(n-1);
        gR=gD/gP;
      }else{
        gD=NA;gP=NA;gR=NA;
      }
      gD.tmp=c(gD.tmp, gD);
      gP.tmp=c(gP.tmp, gP);
      gR.tmp=c(gR.tmp, gR);
    }
    #(5) name gDPAU, and combine to data.frame
    nameArr=function(arr){
      names(arr)=colnames(apaM)
      arr=as.data.frame( t(arr) )
      rownames(arr)=gene;
      return(arr)
    }
    gDPAU=rbind(gDPAU, nameArr(gD.tmp));
    gPPAU=rbind(gPPAU, nameArr(gP.tmp));
    gRatio=rbind(gRatio, nameArr(gR.tmp));
  }
  # check
  if(Debug){
    print(gDPAU[,1:5])
    print(gPPAU[,1:5])
    print(gRatio[,1:5])
  }
  
  return(list('gDPAU'=gDPAU, 'gPPAU'=gPPAU, 'gRatio'=gRatio))
}

### 
GeneralRatio_Long_vs_Short=getGeneralRatio_Long_vs_Short() #time consuming: 10:31 - 10:35



# check
length(GeneralRatio_Long_vs_Short) #3
gDPAU=GeneralRatio_Long_vs_Short[['gDPAU']]
dim(gDPAU) #4103 222
#
gPPAU=GeneralRatio_Long_vs_Short[['gPPAU']]
dim(gPPAU)
#
gRatio=GeneralRatio_Long_vs_Short[['gRatio']]
dim(gRatio)

head(gRatio[,1:10])
#         c12ROW03  c12ROW04  c12ROW05   c12ROW06  c12ROW08  c12ROW15  c12ROW19  c12ROW23  c12ROW25  c12ROW27
#MALAT1  0.4982261 0.9531091 0.4999605  0.5672896 0.4510395 0.6260630 0.6668382 0.7426471 0.6175132 0.4588433
#EGFR    8.1820388 3.3139629 3.0766551 13.6462882 1.4642260 0.9905847 1.0324693 1.7500000 2.9803009 3.7085926

#
table(is.na(gDPAU)) #FALSE   TRUE #475240 435626
table(is.na(gPPAU)) #FALSE   TRUE #475240 435626 
table(is.na(gRatio)) #FALSE   TRUE #475240 435626 

## filter out all NA genes
keep=apply(gRatio, 1, function(x){
  sum(!is.na(x))
})
table(keep<10)
#FALSE  TRUE 
#4096     7
table(keep==0) #no all NA genes
#FALSE 
#4103


#save to files
fn_saveListToFile=function(sourcelist, output){
  write.table(sourcelist[['gDPAU']], paste0(output, "01_generalDPAU.txt"), quote = F)
  write.table(sourcelist[['gPPAU']], paste0(output,"02_generalPPAU.txt"), quote = F)
  write.table(sourcelist[['gRatio']], paste0(output,"03_generalRatio_long_vs_short.txt"), quote = F)
}
fn_saveListToFile(GeneralRatio_Long_vs_Short, output="Matrix_")


# test DPAU and gDPAU is the same, except the range
(function(){
  gene.common=intersect(rownames(DPAU), rownames(gDPAU))
  pd(DPAU)
  pd(gDPAU)
  pl(gene.common);
  
  #
  gene=gene.common[3]
  getNonNA=function(arr){
    arr[!is.na(arr)]
  }
  a1=getNonNA( as.numeric( gDPAU[gene,] ) ) #range [0,1]
  a2=getNonNA( as.numeric( DPAU[gene,] ) ) #range [0,100]
  
  print(a2)
  print(a1*100-a2)
})()


#### plot the Non-NA pct per gene
(function(i){
  #i=5
  arr=gDPAU[i,]
  arr=as.numeric(arr)
  arr=arr[!is.na(arr)] *100
  hist(arr, n=40, main=paste( 'gDPAU of gene: ', rownames(gDPAU)[i] ), 
       xlim=c(0,100), sub=paste('Observation: ', length(arr)) )
})(9)
#








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




##################
# compare: sync vs normal
##################
cid.A=rownames( cellInfo[which(cellInfo$cellType == "HeLa_sync"),] ); length(cid.A) #27
cid.B=rownames( cellInfo[which(cellInfo$cellType == "HeLa_normal"),] ); length(cid.B) #30
#(1)
gDPAU.sync_normal=compareA_B(gDPAU, cid.A, cid.B)
debug(gDPAU.sync_normal) #3551   12
delta_gDPAU_df=(function(){
  df=gDPAU.sync_normal
  df$delta=df$delta*100
  #plot(log2(df$cvA), log2(df$cvB) )
  #plot(df$meanA, df$meanB)
  hist(df$delta, n=100)
  #df=df[which( df$adj.p<0.05  ),]
  ph(df)

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
  debug(tmp)
  tmp
})()
debug(delta_gDPAU_df) #122 12 #126 12
writeLines(rownames(delta_gDPAU_df), '5-gDPAU-delta.gene.txt')


#(2)
DPAU.sync_normal=compareA_B(DPAU, cid.A, cid.B)
debug(DPAU.sync_normal) #1977   12
(function(){
  df=DPAU.sync_normal
  #plot(log2(df$cvA), log2(df$cvB) )
  #plot(df$meanA, df$meanB)
  #hist(df$delta, n=100)
  plot(density(df$delta))
  df=df[which( df$adj.p<0.05  ),]
  
  rownames(df)
})()





##### test on reads
chisqTestA_B=function(cid.A, cid.B){
  # genes with more than 2 apa site
  apaGenes=(function(){
    # filter out 1 site genes
    rs=getDF_fromNamesXX( table(apaSite$gene) )
    rs=rs[order(rs$value),]
    rs=rs[which(rs$value>=2),]
    rs=rs[which(rs$name!='.'),]
    #print(head(rs) )
    genes=rownames(rs); pl(genes); #4103
    genes
  })()
  print(length(apaGenes)) #4103
  #
  i=0
  result=NULL
  gene.all0oneside=c();
  for(gene in apaGenes){
    i=i+1
    #if(i>10){ break; }
    # get apa id
    coord=apaSite[which( apaSite$gene==gene ),'PASid']
    # get matrix
    apa_df.A=apply( apaM[coord, cid.A], 1, sum)
    apa_df.B=apply( apaM[coord, cid.B], 1, sum)
    # filter out all 0 at leat in one side
    if(sum(apa_df.A)==0 | sum(apa_df.B)==0){
      gene.all0oneside=c(gene.all0oneside, gene);
      print(gene); next;
    }
    #
    apa_df=rbind(apa_df.A,apa_df.B)
    # X^2 test
    test=chisq.test(apa_df)
    # result
    rs_df=data.frame(
      gene=gene,
      countsA=sum(apa_df.A),
      countsB=sum(apa_df.B),
      p=test$p.value
    )
    result=rbind(result, rs_df)
  }
  # a adj and order
  result$padj=p.adjust(result$p, method='fdr')
  result=result[order(result$padj), ]
  rownames(result)=result$gene
  return( list(result, gene.all0oneside) )
} # about 30 second
rs=chisqTestA_B(cid.A, cid.B)
length(rs[[2]]) #68
dim(rs[[1]]) #4035    5
#
chisqResult=(function(){
  df=rs[[1]]
  df=df[which(df$padj<0.05), ] #3409
  pd(df)
  head(df)
  df
})()
debug(chisqResult)
#[1] 3409    5
#             gene countsA countsB p padj
#GRB10       GRB10    3556     610 0    0
#HIST3H2A HIST3H2A    3906     390 0    0




###############
# chisqP<0.05 and abs(delta)>20
##############
debug(DPAU.sync_normal)
#[1]  1977   12
#         nA nB     meanA    meanB      sdA       sdB            p       cvA        cvB        adj.p         fc     delta
#GRB10    21 11 96.324719 20.38576 15.743215 39.61901856 3.308477e-06 0.16343899 1.9434650819 0.005349807 4.72509756  75.938956
#HIST3H2A 27  6  8.661464 98.46650 21.18527  3.250418 1.769333e-19 2.4459226 0.03301040 2.861012e-16 0.08796356 -89.80504
deltaDPAU_df=(function(){
  df=DPAU.sync_normal;
  pd(df) #1617 12
  n=length(df$delta) #1617
  m1=mean(df$delta) #0.4719
  s1=sd(df$delta) #0.4719
  pP(n, m1, s1) #38 -37
  #
  t1=m1+1.96*s1
  t2=m1-1.96*s1
  pP(t1,t2)
  plot(density(df$delta))
  abline(v= t1, col='red', lty=2)
  abline(v= t2, col='red', lty=2)
  #
  df2=df[which(df$delta>t1 | df$delta<t2), ]
  pd(df2)
  ph(df2)
  writeLines(rownames(df2), '4-z-score-DPAU.gene.txt')
  write.table(df2, '4-z-score-DPAU.df.txt')
  df2
})()
debug(deltaDPAU_df) #137 12
View(deltaDPAU_df)



(function(){
  df0=chisqResult
  df1=DPAU.sync_normal
  #
  gene.common=intersect(rownames(df0), rownames(df1))
  df0=df0[gene.common, ]
  df1=df1[gene.common, ]
  #
  df2=cbind(df0[,c(1:5)], df1[, c(1:6, 8,9,11,12)])
  #ph(df2)
  #
  n=length(df2$delta) #1216
  m1=mean(df2$delta) #0.449712
  s1=sd(df2$delta) #21.0657106
  pP('n=',n,'; mean=', m1, "; sd=",s1)
  print(m1+2*s1)#42.58113
  print(m1-2*s1)
  df2.sigB=df2[which(df2$padj<0.01 & df2$delta>20 ), ]
  df2.sigS=df2[which(df2$padj<0.01 & df2$delta<(-20) ), ]
  
  #rownames(df)
  debug( rbind(df2.sigB, df2.sigS) ) #1216   15
  writeLines( rownames(df2.sigB), '4-sync_normal_hela.pos.gene.txt')
  writeLines( rownames(df2.sigS), '4-sync_normal_hela.neg.gene.txt')
})()
#



(function(){
  df0=chisqResult
  df1=gDPAU.sync_normal
  df1$delta=df1$delta*100
  #
  gene.common=intersect(rownames(df0), rownames(df1))
  df0=df0[gene.common, ]
  df1=df1[gene.common, ]
  #
  df2=cbind(df0[,c(1:5)], df1[, c(1:6, 8,9,11,12)])
  #ph(df2)
  #
  n=length(df2$delta) #1216
  m1=mean(df2$delta) #0.449712
  s1=sd(df2$delta) #21.0657106
  pP('n=',n,'; mean=', m1, "; sd=",s1)
  print(m1+2*s1)#42.58113
  print(m1-2*s1)
  df2.sigB=df2[which(df2$padj<0.01 & df2$delta>30 ), ]
  df2.sigS=df2[which(df2$padj<0.01 & df2$delta<(-30) ), ]
  
  #rownames(df)
  debug( rbind(df2.sigB, df2.sigS) ) #1216   15
  writeLines( rownames(df2.sigB), '5-gDPAU_sync_normal_hela.pos.gene.txt')
  writeLines( rownames(df2.sigS), '5-gDPAU_sync_normal_hela.neg.gene.txt')
  writeLines( rownames(rbind(df2.sigB, df2.sigS)), '5-gDPAU_sync_normal_hela.all.gene.txt')
})()
#




writeLines( rownames(rsDF[which(rsDF$log2cpm>6 & rsDF$sdDPAU>40),]), "6-sdDPAUhigh_mRNAhigh.gene.txt")

writeLines( rownames(rsDF[which( rsDF$sdDPAU>48),]), "7-sdDPAUhigh.gene.txt")
