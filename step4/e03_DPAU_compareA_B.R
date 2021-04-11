# compare A_B

setwd('/home/wangjl/data/apa/20200701Fig/f4/cycles/')
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
dim(cellInfo) #222 8
head(cellInfo)
#              cid cellType0 cluster cellType cellType2 countsPerCell geneNumber cellCycle
#c12ROW03 c12ROW03   unknown       0     BC_0      BC_0       1820666       6659       MG1
#c12ROW04 c12ROW04   unknown       0     BC_0      BC_0       2816130       7930       MG1
table(cellInfo$cellType)
#BC_0        BC_1 HeLa_normal   HeLa_sync 
#92          73          30          27
cid.BC=row.names(cellInfo[which(substr(cellInfo$cellType,1,2)=='BC'),]);length(cid.BC)
head(cid.BC)

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


#(3) DPAU
DPAU=read.table('/home/wangjl/data/apa/20200701Fig/f4/DPAU/0_DPAU_2082x222.df.txt', row.names = 1)
dim(DPAU) #2082  222
DPAU[1:5,1:5]
#         c12ROW03   c12ROW04 c12ROW05 c12ROW06 c12ROW08
#AARD           NA   0.000000       NA       NA       NA
#AARS            0 100.000000       NA       NA       NA

#(4) gDPAU
gDPAU=read.table('/home/wangjl/data/apa/20200701Fig/f4/DPAU/Matrix_01_generalDPAU.txt', row.names = 1)
dim(gDPAU) # 4103  222
gDPAU[1:5,1:5]
#         c12ROW03  c12ROW04  c12ROW05  c12ROW06  c12ROW08
#MALAT1  0.3325440 0.4879958 0.3333158 0.3619558 0.3108389
#EGFR    0.8910917 0.7681946 0.7547009 0.9317233 0.5941931









######### example 


####
# compare delta of DPAU or gDPAU: sync vs normal
# normalize gRatio to (0,1)
#   method 1(log, (x-min)/(max-min) ),
#   method 2(original raw ratio)
#   method 3: x/(x+1)
compareDeltaA_B=function(dt, cid.A, cid.B){
  dt=dt[, c(cid.A, cid.B)]
  #filter out all NA rows
  keep=apply(dt, 1, function(x){ sum( !is.na(x) )})>0
  dt=dt[keep,]
  pd(dt)
  
  #
  dt=as.data.frame(dt)
  genes= rownames(dt);    
  df=NULL;
  j=0
  for(i in 1:length(genes)){
    gene=genes[i]
    arrA=as.numeric(dt[gene, cid.A])
    arrB=as.numeric(dt[gene, cid.B])
    # filter out NA
    arrA=arrA[!is.na(arrA)]
    arrB=arrB[!is.na(arrB)]
    if(length(arrA)==0 | length(arrB)==0){
      #pP('>>>Discard(either lenght is 0):',i, ' ', gene)
      next;
    }
    p=1
    #if(length(arrA)>=2 & length(arrB)>=2){
      #if(sd(arrA)!=0 & sd(arrB)!=0){
        p=wilcox.test(arrA, arrB)$p.value;
        j=j+1
      #}
    #}
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
  df=df[which(!is.na(df$p)),]
  df$adj.p=p.adjust(df$p, method ='fdr')
  df=df[order(df$adj.p),]
  df$fc=df$meanA / df$meanB
  df$delta=df$meanA - df$meanB
  df
}






##################
# compare delta: sync vs normal
##################
cid.A=rownames( cellInfo[which(cellInfo$cellType == "HeLa_sync"),] ); length(cid.A) #27
cid.B=rownames( cellInfo[which(cellInfo$cellType == "HeLa_normal"),] ); length(cid.B) #30
#(1)
DPAU.sync_normal=compareDeltaA_B(DPAU, cid.A, cid.B) #1982
debug(DPAU.sync_normal) #1982   12
View(DPAU.sync_normal)
(function(){
  df=DPAU.sync_normal
  plot(density(df$delta))
  df=df[which( df$adj.p<0.05  ),]
  
  pd(df)
  rownames(df)
})()

drawHist=function(df){
  arr=df$delta; m1=mean(arr); s1=sd(arr); t1=m1-2*s1; t2=m1+2*s1;
  pP(m1, s1, ' || ', t1, t2)
  plot(density(arr));  abline(v=c(t1,t2), col='red', lty=2)
  df[which(df$delta<t1 | df$delta>t2), ]
}
tmp=drawHist(DPAU.sync_normal)
writeLines(rownames(tmp), '01-deltaDPAU.gene.txt')

#(2)
gDPAU.sync_normal=compareDeltaA_B(gDPAU, cid.A, cid.B)
gDPAU.sync_normal$delta=gDPAU.sync_normal$delta*100
debug(gDPAU.sync_normal) #3940   12
#       nA nB     meanA     meanB        sdA         sdB            p        cvA         cvB        adj.p        fc      delta
#EIF3K  27 30 0.7508131 0.9992413 0.30973699 0.002001580 4.957500e-11 0.41253539 0.002003100 1.953255e-07 0.7513832 -24.842814
View(gDPAU.sync_normal)
tmp2=drawHist(gDPAU.sync_normal)
writeLines(rownames(tmp2), '01-delta_gDPAU.gene.txt')
# q value of GO >0.05, maybe combine p value will be better.







##################
# compare X^2 p value: sync vs normal
##################
##### test on reads
chisqTestA_B=function(cid.A, cid.B){
  pP('Cell Number: ', length(cid.A), length(cid.B))
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
      #print(gene); 
      next;
    }
    # X^2 test
    apa_df=rbind(apa_df.A, apa_df.B)
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
  print(paste('Abandon gene number(low counts):', length(gene.all0oneside) ) )
  return( list(result, gene.all0oneside) )
} # about 30 second

rs=chisqTestA_B(cid.A, cid.B)
length(rs[[2]]) #68
dim(rs[[1]]) #4035    5
#
chisqResult=(function(){
  df=rs[[1]]
  df=df[order(df$padj), ]
  #df=df[which(df$padj<0.05), ] #3409
  pd(df) #4035    5
  head(df)
  df
})()
debug(chisqResult)


# combine p value, and delta or fc
volcanoPlot=function(df1, df2, keyword, cid.A,cid.B, output=''){
  #df1=gDPAU.sync_normal;
  #df2=rs[[1]];
  gene.common=intersect(rownames(df1), rownames(df2))
  df1=df1[gene.common, ]
  df2=df2[gene.common, ]
  #
  df=cbind(df1[, -c(7)], df2[, -c(4)])
  debug(df) #3940   16
  df$sig="ns"
  df[which( (df$delta)>30  & df$padj<0.05), ]$sig="up"
  df[which( (df$delta)<(-30)  & df$padj<0.05), ]$sig="down"
  # save to file
  write.table(df, paste0(output, '02_volcano_gDPAU_chisqP_',keyword,'.txt'))
  writeLines( rownames(df[which(df$sig=='up'),]), paste0(output, '02_volcano_gDPAU_chisqP_',keyword,'.up.gene.txt'))
  writeLines( rownames(df[which(df$sig=='down'),]), paste0(output, '02_volcano_gDPAU_chisqP_',keyword,'.down.gene.txt'))
  writeLines( rownames(df[which(df$sig!='ns'),]), paste0(output, '02_volcano_gDPAU_chisqP_',keyword,'.All.gene.txt'))
  #
  tb=table(df$sig); print(tb)
  CairoPDF(paste0(output, '02_volcano_gDPAU_chisqP_',keyword,'.pdf'),width=3.5,height=4)
  g=ggplot(df, aes(delta, -log10(padj), color=factor(sig, levels=c('down', 'ns', 'up')) ) )+
    geom_point(size=0.1, alpha=0.5)+theme_bw()+
    labs( title=paste("gDPAU:", keyword, length(cid.A),length(cid.B)), x='Change of gDPAU', y="-log10(adj.p)" ) + theme_bw()+
    theme(legend.box = "horizontal", # 图例，水平，标到底部
          legend.key.size=unit(6,"pt"), #图例方块的高度
          legend.position="bottom") +
    scale_color_manual('', labels=c( paste0("shorten(",tb[1],")"),paste0('n.s.(',tb[2],')'),paste0("lengthen(",tb[3],")") ), #图例文字
                       values=c("blue", "#dddddd",'red') ) #自定义颜色
  
  dd_text=rbind(
    (function(){
      df2=df[which(df$sig=='up'),]
      head(df2[order(df2$padj), ])
    })(),
    (function(){
      df2=df[which(df$sig=='down'),]
      head(df2[order(df2$padj), ])
    })()
  );pd(dd_text)
  g2=g+geom_text_repel(data=dd_text, aes(x=delta, y=-log10(padj), label=rownames(dd_text)),
                       color="black",size=3,alpha=0.6)+
    guides(colour = guide_legend(override.aes = list(alpha = 1,size=2)))#放大图例的点
  
  print(g2)
  dev.off()
  return(df)
}
#df1=gDPAU.sync_normal;
#df2=rs[[1]];
#volcanoPlot(gDPAU.sync_normal, rs[[1]],'sync_vs_normal_HeLa',cid.A,cid.B,)
#



############
# combine to 1 function
compare_A_B=function(cid.A, cid.B, keyword='sync_vs_normal_HeLa', output="."){
  #1. delta
  print('############### part 1: get delta')
  gDPAU.A_B=compareDeltaA_B(gDPAU, cid.A, cid.B)
  gDPAU.A_B$delta=gDPAU.A_B$delta*100
  
  #2. p value
  print('############### part 2: get p value')
  rs=chisqTestA_B(cid.A, cid.B)
  print(length(rs[[2]])) #68
  pd(rs[[1]]) #4035    5
  
  #3. volcanoPlot
  print('############### part 3: plot')
  volcanoPlot(gDPAU.A_B, rs[[1]],keyword,cid.A,cid.B, output)
}

#(1) test: sync vs normal HeLa
cid.A=rownames( cellInfo[which(cellInfo$cellType == "HeLa_sync"),] ); length(cid.A) #27
cid.B=rownames( cellInfo[which(cellInfo$cellType == "HeLa_normal"),] ); length(cid.B) #30
chisqResult=compare_A_B(cid.A, cid.B, keyword="sync_normal", output='sync_normal/')
head(chisqResult)
#




#(0) test: BC vs HeLa
cid.A=rownames(cellInfo[which( substr(cellInfo$cellType,1,2)=="BC"), ]); length(cid.A) #165
cid.B=rownames(cellInfo[which( substr(cellInfo$cellType,1,2)=="He"),]); length(cid.B) #57
pP('cell number:',length(cid.A), length(cid.B)) #165 57
#
chisqResult=compare_A_B(cid.A, cid.B, keyword="BC_HeLa", output='gDPAU_BC_HeLa/')
head(chisqResult)
tail(chisqResult)
#

# BC vs HeLa_normal
table(cellInfo$cellType)
cid.A=rownames(cellInfo[which( substr(cellInfo$cellType,1,2)=="BC"), ]); length(cid.A) #165
cid.B=rownames(cellInfo[which( cellInfo$cellType=="HeLa_normal"),]); length(cid.B) #30
pP('cell number:',length(cid.A), length(cid.B)) #165 30
#
chisqResult=compare_A_B(cid.A, cid.B, keyword="BC_HeLaNormal", output='gDPAU_BC_HeLaNormal/')
head(chisqResult)
tail(chisqResult)
#















#
# (2) BC, BC_0 vs BC_1
cid.A=rownames(cellInfo[which(cellInfo$cellType=="BC_0"),])
cid.B=rownames(cellInfo[which(cellInfo$cellType=="BC_1"),])
pP('cell number:',length(cid.A), length(cid.B))
#
chisqResult=compare_A_B(cid.A, cid.B, keyword="BC0_vs_BC1", output='gDPAU_BC/')
head(chisqResult)



#(3) BC, cell cycle pair-wise compare
(function(){
  cycles=rev( c('G1S', 'S', 'G2M', 'M', 'MG1') )
  cycles=c(cycles, cycles)
  print(cycles)
  for(i in 1:5){
    pP('==============>>>>', i, cycles[i], cycles[(i+1)])
    p1=cycles[i]; p2=cycles[(i+1)]
    #
    cid.A=rownames(cellInfo[which( substr(cellInfo$cellType,1,2)=="BC" & cellInfo$cellCycle==p1),])
    cid.B=rownames(cellInfo[which( substr(cellInfo$cellType,1,2)=="BC" & cellInfo$cellCycle==p2),])
    #
    pP('-------------->>>>cell number:',length(cid.A), length(cid.B))
    chisqResult=compare_A_B(cid.A, cid.B, keyword=paste0('BC_',p1,"_vs_",p2), output='gDPAU_BC/')
  }
})()


#(4) HeLa, cell cycle pair-wise compare
(function(){
  cycles=rev( c('G1S', 'S', 'G2M', 'M', 'MG1') )
  cycles=c(cycles, cycles)
  print(cycles)
  for(i in 1:5){
    pP('==============>>>>', i, cycles[i], cycles[(i+1)])
    p1=cycles[i]; p2=cycles[(i+1)]
    #
    cid.A=rownames(cellInfo[which( substr(cellInfo$cellType,1,2)=="He" & cellInfo$cellCycle==p1),])
    cid.B=rownames(cellInfo[which( substr(cellInfo$cellType,1,2)=="He" & cellInfo$cellCycle==p2),])
    #
    pP('-------------->>>>cell number:',length(cid.A), length(cid.B))
    chisqResult=compare_A_B(cid.A, cid.B, keyword=paste0('HeLa_',p1,"_vs_",p2), output='gDPAU_HeLa/')
  }
})()





##############
## combine genes
# load all genes
geneList.BC=list()
cycles=rev( c('G1S', 'S', 'G2M', 'M', "MG1") )
cycles2=c(cycles,cycles)
for(i in 1:length(cycles)){
  keyword=paste0(cycles2[i], '_vs_',cycles2[i+1])
  fileName=paste0("gDPAU_BC/02_volcano_gDPAU_chisqP_BC_",keyword,".All.gene.txt")
  #
  printP(i, keyword, fileName)
  geneList.BC[[keyword]]= c( readLines(fileName) )
}
sapply(geneList.BC, length) 
#geneList.BC

# load all genes
geneList.HeLa=list()
for(i in 1:length(cycles)){
  keyword=paste0(cycles2[i], '_vs_',cycles2[i+1])
  fileName=paste0("gDPAU_HeLa/02_volcano_gDPAU_chisqP_HeLa_",keyword,".All.gene.txt")
  #
  printP(i, keyword, fileName)
  geneList.HeLa[[keyword]]= c( readLines(fileName) )
}
sapply(geneList.HeLa, length) 
#geneList.HeLa




##########################
# pool all genes: BC
getUniq_genes_VennPlot=function(geneList, output, keyword){
  #geneList=geneList.BC
  #output="gDPAU_BC/"
  #keyword="BC"
  
  ##
  pP('Total gene number:', sum(sapply(geneList, length)) )
  
  genesDF=data.frame(
    'group'=names(unlist(geneList)),
    'genes'=unname(unlist(geneList))
  )
  dim(genesDF)
  head(genesDF)
  
  # uniq genes
  genes=unique(as.character( genesDF$genes) )
  pP("Uniq gene number:",length(genes)) #834
  ph(genes)
  
  
  ################
  # Venn plot
  cols=RColorBrewer::brewer.pal(n = 5,name = "Set2")
  #
  library("VennDiagram")
  grid.newpage()
  n=1800
  venn.plot <- venn.diagram(
    x =geneList,
    filename = NULL,
    #filename = paste0(output, "Venn/Cycle_gene_APA_level_Venn_5set.png" ), #dpi=300,
    resolution =200, 
    imagetype = "png", units = "px",height = n, width = n, 
    
    col = "black",
    fill = cols, #c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    alpha = 0.50,
    cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
            1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
    cat.col = cols,#c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    cat.cex = 1.5,
    cat.fontface = "bold",
    margin = 0.5,
    main=paste0("Cycle related genes of ",keyword," at APA level")
  );
  
  pdf(file=paste0(output, "Venn/Cycle_gene_APA_level_Venn_5set.pdf" ), width=5, height=5 )
  grid.draw(venn.plot)
  dev.off()
  
  #########
  writeLines(genes, paste0(output, keyword,'_cellCycle_related_gene_APAlevel_uniq.gene.txt'))
  return(genes);
}

cycleGenes.BC=getUniq_genes_VennPlot(geneList.BC, "gDPAU_BC/", "BC")

cycleGenes.HeLa=getUniq_genes_VennPlot(geneList.HeLa, "gDPAU_HeLa/", "HeLa")



## intersect genes
(function(){
  g1=cycleGenes.BC
  g2=cycleGenes.HeLa
  #     # union
  #     cycleGenes.union=union( g1, g2)
  #     length(cycleGenes.union)
  #     head(cycleGenes.union) # 568
  #     #writeLines(cycleGenes.union, 'BC_HeLa_CycleRelatedGene_union_568.txt')
  
  pl(g1)
  pl(g2)
  #
  g3=intersect(g1, g2) #cycleGenes.intersect
  pl(g3)
  writeLines(g3, 'BC_HeLa_CycleRelatedGene_intersect.txt')
  
  library (VennDiagram)
  #grid.newpage()
  pP( length(g1), length(g2), length(g3) )
  venn.plot <- draw.pairwise.venn( length(g1), length(g2), length(g3), 
                                   c("BC", "HeLa"), 
                                   
                                   fill = c('#F81082', '#619CFF'),
                                   col = "transparent",
                                   cex = 1.5,cat.cex=1.5,
                                   alpha = 0.50,
                                   main="Cycle related genes at APA level",
                                   
                                   cat.pos = c(-30, 10), #angle of text, dgree=0 means top of cycle
                                   cat.dist = 0.05,   #distance of text to border, can be negative
                                   
                                   cat.col = c('#F81082', '#619CFF'),
                                   margin=0.3,
                                   scaled = T);
  pdf(file="Venn/BC_HeLa_CycleRelatedGene_intersect_venn.pdf", width=4, height=4)
  grid.draw(venn.plot)
  dev.off()
  print('==end==')
})()
#













##################################
# pseudo-time plot
##################################
## step1: get total gene's mean gDPAU at cell cycle time point
get_timeDPAU=function(df, output, keyword){
  if(!( keyword %in% c("BC", "HeLa")) ){
    stop('Error: keyword must in c("BC", "HeLa")');
  }
  #df=gDPAU.BC[cycleGenes.BC, cid.BC]
  pd(df) 
  
  ###############
  ## cell cycle time point, mean DPAU, sd DPAU, num
  cycles=c('G1S', 'S', 'G2M', 'M', "MG1")
  timeDPAU=NULL;
  
  num=1
  for(gene in rownames(df)){
    #if(num>20){break;}
    
    tmp_DF=NULL;
    for(cycle in cycles){
      num=num+1
      #if(num>20){break;}
      #
      cid.tmp=rownames(cellInfo[which( substr(cellInfo$cellType,1,2)==substr(keyword,1,2) & cellInfo$cellCycle==cycle ), ] )
      #printP(cycle, length(cid.tmp)) #check
      DPAU.tmp=as.numeric( df[gene, cid.tmp] )
      DPAU.tmp=DPAU.tmp[!is.na(DPAU.tmp)]
      #
      tmp=data.frame(
        'a'=mean(DPAU.tmp),
        'b'=sd(DPAU.tmp),
        'c'=length(DPAU.tmp)
      );
      
      
      colnames(tmp)=paste0(cycle, c('_mean', '_sd', '_n') )
      #
      if(is.null(tmp_DF)){tmp_DF=tmp}else{
        tmp_DF=cbind(tmp_DF, tmp );
      }        
    }
    tmp_DF$gene=gene;
    timeDPAU=rbind(timeDPAU, tmp_DF)
  }
  ##
  rownames(timeDPAU)=timeDPAU$gene
  pd(timeDPAU);ph(timeDPAU, n=2)
  #save to file
  write.table(timeDPAU, paste0(output, keyword, '_mean_gDPAU_at_cellCycle_timePoint.df.txt'))
  pP('saved at:', output, keyword)
  return(timeDPAU)
}
##
timeDPAU.BC0=get_timeDPAU(gDPAU[cycleGenes.BC, ], "gDPAU_BC/", 'BC')
timeDPAU.HeLa0=get_timeDPAU(gDPAU[cycleGenes.HeLa, ], "gDPAU_HeLa/", 'HeLa')




##################################
# step2: filter out genes whose observation are less than 2 at any time point;
filter_timeDPAU=function(df0, output, keyword){
  #df0=timeDPAU.BC
  pP('1. input: ', nrow(df0))
  
  # filter by cell number, at least 2 per time point
  df_n=df0[, seq(3,15,3)]
  lowCellID=apply(df_n, 1, function(x){
    sum(x<2)
  })
  print(table(lowCellID))
  #return(keep)
  df0=df0[!lowCellID,]
  pP("2. after filter out by cell number(<2):", nrow(df0))
  
  #get mean DPAU at each phase
  df=df0[,seq(1,15,3)]
  pd( df[complete.cases(df), ] )
  pP('3. after filter out NA:', nrow(df) )
  
  # simplify col names
  colnames(df)=gsub('_mean','',colnames(df))
  debug(df)
  # pheatmap( as.matrix( df ) )
  #save to file
  write.table(df, paste0(output, keyword, '03_mean_gDPAU_at_cellCycle_timePoint.filterGenesLessThan2cell.df.txt'))
  pP('saved at:', output, keyword)
  return(df)
}

#
timeDPAU.BC=filter_timeDPAU(timeDPAU.BC0, "gDPAU_BC/", 'BC'); dim(timeDPAU.BC)
timeDPAU.HeLa=filter_timeDPAU(timeDPAU.HeLa0, "gDPAU_HeLa/", 'HeLa'); dim(timeDPAU.HeLa)
# overlap nenn plot: 569 257 555  
(function(){
  g1=rownames(timeDPAU.BC)
  g2=rownames(timeDPAU.HeLa)
  g3=intersect( g1, g2 )
  
  pP( length(g1), length(g2), length(g3) )
  writeLines(g3, '03_BC_HeLa_CycleRelatedGene_intersect.filterLessThan2cells.txt')
  venn.plot <- draw.pairwise.venn( length(g1), length(g2), length(g3), 
                                   c("BC", "HeLa"), 
                                   fill = c('#F81082', '#619CFF'),
                                   col = "transparent",
                                   cex = 1.8,cat.cex=1.8,alpha = 0.50,
                                   main="Cycle related genes at APA level, \n
                            after filter out genes whose observation are less than 2 cells at any cell phase time point",
                                   
                                   cat.pos = c(-30, 10), #angle of text, dgree=0 means top of cycle
                                   cat.dist = 0.05,   #distance of text to border, can be negative
                                   
                                   cat.col = c('#F81082', '#619CFF'),
                                   margin=0.3,
                                   scaled = T);
  pdf(file="Venn/BC_HeLa_CycleRelatedGene_intersect_filtered_venn.pdf", width=4, height=4)
  grid.draw(venn.plot)
  dev.off()
  print('==end==')
})()
#






######
# order gene by correlation
# cluster genes accoring to their correlation across cell cycle time point
#   - gene correlation heatmap
drawHM01=function(df, output, keyword){
  #df=timeDPAU.BC
  #output="gDPAU_BC/"
  #keyword="BC"
  pd(df)
  #
  clusterNum=6 ##----------todo: set parameters
  
  
  ## using correlation
  df.cor=cor( t(as.matrix( df ) ), use = 'pairwise.complete.obs', 
              method="spearman"
              #method="pearson"
  )
  #diag(df.cor)=NA;
  pHM=pheatmap( df.cor, 
                clustering_method = 'ward.D2', na_col = 'grey',
                #clustering_distance_rows = "correlation",
                #scale='row',
                show_rownames=F,show_colnames=F,
                
                cutree_cols = clusterNum,
                cutree_rows = clusterNum,
                #treeheight_row = 100,
                
                main=paste0("Correlation of cell cycle related genes \nAt APA level(", keyword,") ", nrow(df),' genes' ),
                border=FALSE)
  colors.cluster=c("#E41E25", "#FBD800", "#208A41", "#446CFF", "#36AFE5", "#AF4A9C")
  
  #dev.off()
  CairoPDF( paste0(output, keyword, '04_heatmap_cycleGene.pdf'), width=4.5, height=4.5)
  # pic1
  clust_ward <- hclust( dist(df.cor) , method = 'ward.D2')
  plot(clust_ward)
  rect.hclust(clust_ward, k = clusterNum, border =colors.cluster )
  # pic2
  grid.newpage()
  print(pHM)
  dev.off()
  
  pP('draw end: ', output, keyword)
  return(pHM)
}
pheatmapRS.BC=drawHM01(timeDPAU.BC, "gDPAU_BC/", "BC" )


# order from heatmap
getNewOrder=function(pheatmapRS, df, output="gDPAU_BC/",keyword="BC"){
  clusterNum=6
  # get gene order after cluster
  row_cluster = cutree(pheatmapRS$tree_row,k=clusterNum)
  # reorder data frame, add Cluster column
  newOrder = df[pheatmapRS$tree_row$order,]
  newOrder[,ncol(newOrder)+1]= row_cluster[match(rownames(newOrder),names(row_cluster))]
  colnames(newOrder)[ncol(newOrder)]="Cluster"
  #order by cluster
  newOrder=newOrder[order(newOrder$Cluster),]
  #head(newOrder)
  
  # how many clusters?
  print(unique(newOrder$Cluster))
  # how many genes each cluster?
  print(table(newOrder$Cluster))
  # 1   2   3   4   5  6
  #116 189 206  99  80 122 
  #
  # write df to file
  newOrder$Cluster = paste0("Cluster",newOrder$Cluster)
  write.table(newOrder, paste0(output, keyword,"05_gDPAU_cycleRelatedGenes.pheatmap.cluster.df.txt"),
              sep="\t",quote = F,row.names = T,col.names = T)
  return(newOrder)
}

#
newOrder.BC=getNewOrder(pheatmapRS.BC, timeDPAU.BC )
dim(newOrder.BC) #812 6
head(newOrder.BC)



# plot heatmap for newOrdered data with Cluster number - 
#    long heatmap
plotHM02=function(newOrder, output, keyword, is_cluster_rows=F){
  # anno for rows
  annotation_row=data.frame(
    gene=rownames(newOrder),
    Cluster=newOrder$Cluster, 
    row.names=1
  )
  #head(annotation_row)
  
  # define color list
  colors.cluster=c("#E41E25", "#FBD800", "#208A41", "#446CFF", "#36AFE5", "#AF4A9C")
  col=colors.cluster # RColorBrewer::brewer.pal(n = 8,name = "Set2")
  #
  num=length(unique(newOrder$Cluster))
  Cluster=col[1:num]
  names(Cluster)=paste0('Cluster',1:num)
  #
  ann_colors = list( Cluster = Cluster )
  #print(head(ann_colors))
  #
  gaps_row =(table(newOrder$Cluster))
  gaps_row=cumsum(gaps_row)
  #
  ######################
  ## plot gene's DPAU across cycles
  ######################
  fileName=paste0(output,keyword,'06_heatmap_cycle_gDPAU_long.pdf')
  
  CairoPDF(fileName, width=ifelse(is_cluster_rows, 5,3), height=7)
  p1.DPAU=pheatmap( newOrder[,1:5], clustering_method = 'ward.D2', scale='row',
                    clustering_distance_rows = "correlation",
                    annotation_row =annotation_row,
                    annotation_colors = ann_colors,
                    #
                    gaps_row=gaps_row,
                    #scale='row',
                    show_rownames=F, #show_colnames=F,
                    cluster_rows=is_cluster_rows,  #--------------<<< 是否对基因聚类
                    cluster_cols=F,
                    cutree_rows = 6,
                    treeheight_row = 100,
                    main=paste0(" gDPAU of cell cycle related \ngenes at APA level: ", keyword),
                    angle_col=90,
                    border=FALSE)
  dev.off()
  #
  printP('heatmap end: ',fileName)
}
plotHM02(newOrder.BC,output="gDPAU_BC/",keyword="BC", is_cluster_rows=F)
#





##############
# plot gene by gene cluster
##############
# add gene number after cluster number
addGeneNumber2Cluster=function(newOrder){
  #print(head(newOrder))
  mapCN=sapply(split(newOrder, newOrder$Cluster), nrow)
  #print( mapCN )
  
  for(C in names(mapCN)){
    N=unname(mapCN[C])
    C2=sub('Cluster', 'Cluster ', C)
    newOrder[which(newOrder$Cluster==C), ]$Cluster=paste0(C2, " (",N, " genes)" )
  }
  newOrder
}
# test
head( addGeneNumber2Cluster(newOrder.BC) )


# plot genes by cluster
# v0.2
plotLineByCluster=function(newOrder, output="gDPAU_BC/", keyword='BC'){
  # change from [0,1] to [0,100]
  newOrder.2=as.data.frame( apply(newOrder[,1:5],2, function(x){x*100}) )
  newOrder.2$Cluster=newOrder[,6]
  newOrder=newOrder.2
  
  #add gene column
  newOrder$gene = rownames(newOrder)
  
  # add gene numbers per Cluster
  newOrder=addGeneNumber2Cluster(newOrder)
  
  #change factor to number
  cycleName=colnames(newOrder.BC)
  colnames(newOrder.BC)=c(1:5,"Cluster")
  
  #wide to long
  library(reshape2)    
  data_new = melt(newOrder)
  
  #dev.off()
  CairoPDF( paste0(output, keyword, '07_clusterLinesByCycle.pdf'), width=3, height=7)
  # line plot for each cluster
  colors.cluster=c("#E41E25", "#FBD800", "#208A41", "#446CFF", "#36AFE5", "#AF4A9C")
  p=ggplot(data_new,aes( as.numeric( variable ), value, group=gene)) + geom_line(color="gray90",size=0.4) + 
    geom_hline(yintercept =50,linetype=2, size=0.1) +
    geom_vline(xintercept =c(2,3,4),linetype=2, size=0.01, color="#999999") +
    #stat_summary(aes(group=1, color=Cluster),fun.y=mean, geom="line", size=0.8, color="#F81082" )+
    stat_summary(aes(group=1, color=Cluster),fun.y=mean, geom="line", size=0.8, #color="#c51b7d"
    ) + scale_color_manual('Gene cluster',values = colors.cluster)+
    facet_wrap(Cluster~., ncol=1) +
    labs(x="Cell cycle", y='Mean gDPAU', title=paste0("Gene cluster: ", keyword))+
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text = element_text(size=8, face = "bold"),
          axis.text.x=element_text(angle=60, hjust=1), #x标签旋转60度
          legend.position = 'none', 
          strip.text = element_text(size = 8, face = "bold"))+
    scale_x_discrete(expand = c(0,0), limits=1:5, labels = cycleName)
  print(p)
  dev.off()
  
  
  ## output gene list for each cluster
  #head(newOrder)
  num=length(unique(newOrder$Cluster))
  for(i in 1:num){
    CName=unique(newOrder$Cluster)[i]
    gene.tmp=rownames( newOrder[which(newOrder$Cluster==CName ),] )
    printP(i, length(gene.tmp))
    writeLines(gene.tmp, paste0(output, keyword, "07_cycleCluster",i,".gene.txt") )
  }
}
## 
plotLineByCluster( newOrder.BC, output="gDPAU_BC/", keyword='BC')
writeLines(rownames(newOrder.BC), 'gDPAU_BC/BC07_cycleClusterAll.gene.txt')
#





################
# for HeLa: cell number too low, so not reliable.
################
pheatmapRS.HeLa=drawHM01(timeDPAU.HeLa, "gDPAU_HeLa/", "HeLa" )
#
newOrder.HeLa=getNewOrder(pheatmapRS.HeLa, timeDPAU.HeLa )
#  1   2   3   4   5   6 
#160 208  90  82 136 150
dim(newOrder.HeLa) #826 6
head(newOrder.HeLa)
#
plotHM02(newOrder.HeLa,output="gDPAU_HeLa/",keyword="HeLa", is_cluster_rows=F)
#
plotLineByCluster(newOrder.HeLa,  output="gDPAU_HeLa/", keyword='HeLa')
writeLines(rownames(newOrder.HeLa), 'gDPAU_HeLa/HeLa07_cycleClusterAll.gene.txt')
#








########## heatmap
tmp=as.matrix(gDPAU.BC)*100
tmp[is.na(tmp)]=-20
#View(tmp)

# col annotation info
annotation_col = data.frame(
  cellType = as.character(cellInfo[cid.BC,]$cellType),
  cellCycle=factor(cellInfo[cid.BC,]$cellCycle, levels = c('G1S', 'S', 'G2M', 'M', 'MG1') )
)

rownames(annotation_col) = cellInfo[cid.BC,]$cid
head(annotation_col)
str(annotation_col)
#
col=RColorBrewer::brewer.pal(n = 5,name = "Set2")
#barplot(1:5,col=col)

# colors
ann_colors = list(
  cellType = c('BC_0'="#FF9ECE", 'BC_1'="#F81082", 'HeLa_normal'="#005FFF", 'HeLa_sync'="#619CFF"),
  cellCycle=c('G1S'=col[1], 'S'=col[2], 'G2M'=col[3],   'M'=col[4], 'MG1'=col[5])
)
head(ann_colors)

# load cycle phase
order.BC=read.table('~/data/apa/20200701Fig/f3/cell_cycle/result/BC_PhaseRefCor.txt', row.names = 1)
dim(order.BC)
head(order.BC)
#View(order.BC)


#
CairoPDF('gDPAU_BC/BC_gDPAU_heatmap.pdf', width=5, height=7)
pheatmap(tmp[,rownames(order.BC)], #scale='row',
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         main = '1',
         annotation_col =annotation_col,
         annotation_colors = ann_colors,
         color = c(colorRampPalette(colors = c("grey","grey"))(20),colorRampPalette(colors = c("blue",'white',"red"))(100) ),
         clustering_method = 'ward.D2')
dev.off()
#
dim(tmp)#4103
dim(newOrder.BC) #812 6
CairoPDF('gDPAU_BC/BC_gDPAU_heatmap.smallSet-cluster_cols.pdf', width=5, height=4)
pheatmap(tmp[ rownames(newOrder.BC), rownames(order.BC)], #scale='row',
         cluster_cols = T,
         show_rownames = F,
         show_colnames = F,
         main = '1',
         annotation_col =annotation_col,
         annotation_colors = ann_colors,
         color = c(colorRampPalette(colors = c("grey","grey"))(20),colorRampPalette(colors = c("blue",'white',"red"))(100) ),
         clustering_method = 'ward.D2')
dev.off()
#






##########################
# get intersect genes in 5 phases!
(function(){
  # load all genes
  geneList=geneList.BC
  
  cycles=rev( c('G1S', 'S', 'G2M', 'M', "MG1") )
  cycles2=c(cycles,cycles)
  geneCommon=geneList[['MG1_vs_M']]
  for(i in 1:length(cycles)){
    keyword=paste0(cycles2[i], '_vs_',cycles2[i+1])
    geneCommon=intersect(geneCommon, geneList[[keyword]])
  }
  #sapply(geneList, length)
  geneCommon
})()
## BC: [3] "DENND3" "UBE3C"  "GNPTAB"
## HeLa: [15] "RBL2"  "VPS33A"  "AL390719.1" "CCNY" "CEP63" "CTSC"  "DDX10"     
# "GTPBP8" "KDM4A" "MAPKAP1"  "NFRKB" "NIPSNAP3A" "PFKFB2"  "SH3KBP1" "TFB1M"

#
gDPAU.BC=gDPAU[,rownames(cellInfo[which(substr(cellInfo$cellType,1,2)=="BC"), ])]
dim(gDPAU.BC)

plotgDPAU_ByGene=function(df0, g){
  # gDPAU.BC['DENND3',!is.na(gDPAU.BC['DENND3',])]
  #df0=gDPAU.BC
  df=df0[g,]
  df=df[g, !is.na(df[g,])]
  df=as.data.frame(t(df*100))
  colnames(df)=c('gDPAU')
  # add cycle
  df$phase=cellInfo[rownames(df), ]$cellCycle
  df=df[order(df$phase),]
  p=ggplot(df, aes(phase, gDPAU))+
      geom_boxplot()+geom_jitter()+theme_bw()+
      labs(title=g)+scale_x_discrete(limits=c('G1S', 'S', 'G2M', 'M', "MG1"))
  plot(p)
  df
}
plotgDPAU_ByGene(gDPAU.BC,'DENND3')
plotgDPAU_ByGene(gDPAU.BC,'UBE3C')
plotgDPAU_ByGene(gDPAU.BC,'GNPTAB')
plotgDPAU_ByGene(gDPAU.BC,'HIST3H2A')
#
plotgDPAU_ByGene(DPAU[,cid.BC],'HIST3H2A')
plotgDPAU_ByGene(DPAU[,cid.BC],'AGPS')


## check
tmp=gDPAU['HIST3H2A', cid.BC]
tmp[,!is.na(tmp)]*100
#
tmp2=apaSite[which(apaSite$gene=='HIST3H2A'),]
rownames(tmp2)
apaM[rownames(tmp2),'c16ROW37']
