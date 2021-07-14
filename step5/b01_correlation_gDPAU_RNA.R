# library(here)
setwd("/data/jinwf/wangjl/apa/20200701Fig/f6/result/")
getwd()

library(ggplot2)
library(Cairo)

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
# 92          73          30          27 
cid.BC=rownames(cellInfo[which( substr(cellInfo$cellType,1,2)=='BC'),]); length(cid.BC) #165
cid.HeLa=rownames(cellInfo[which( substr(cellInfo$cellType,1,2)=='He'),]); length(cid.HeLa) #57
#
#(2) add RNA read counts for each genes
rnaM=read.csv('/data/jinwf/wangjl/apa/20200701Fig/f2/BC_HeLa/BC_HeLa.222cells.count.V4.csv', row.names = 1)
dim(rnaM) #18662 222
rnaM[1:10,1:10]
#
rnaM.logcpm=as.data.frame(apply(rnaM, 2, function(x){
  log2(x/sum(x)*1e6+1)
}))
dim(rnaM.logcpm)
rnaM.logcpm[1:10,1:10]
#




#(3) gDPAU
gDPAU=read.table('/home/wangjl/data/apa/20200701Fig/f4/DPAU/Matrix_01_generalDPAU.txt', row.names = 1)
dim(gDPAU) # 4103  222
gDPAU[1:5,1:5]




###################
# process1: get apa RNA relation
###################
source("../src/f1_correlaton_gDPAU_RNA.src.R")
df_cor.BC=get_apa_rna_cor(cid.BC)
#df_cor.HeLa=get_apa_rna_cor(cid.HeLa)
#df_cor.all=get_apa_rna_cor( c(cid.BC, cid.HeLa) ) #time consuming:1min
#
#df_cor.all['CSTF2',]



###################
# process2: plot apa RNA point plot
###################
source("../src/f1_correlaton_gDPAU_RNA.src.R")
#save_plot_cor(df_cor.BC, 'BC', 'cor_BC/')
#save_plot_cor(df_cor.HeLa, 'HeLa', 'cor_HeLa/')
#save_plot_cor( df_cor.all, 'all', 'cor_all/')
#


save_plot_cor2=function(df_cor, keyword="BC", output='cor_BC/'){
  #define significance
  print(dim(df_cor)) 
  #order
  df_cor=df_cor[order(-df_cor$cor),]
  # filter out too few cell items
  df_cor=df_cor[which(df_cor$n>=10),]
  print(dim(df_cor))
  #
  # (1) positive cor
  df_cor_pos=df_cor[which(df_cor$cor>0 & df_cor$adj.p<0.05),]
  df_cor_min=df_cor[which(df_cor$cor<0 & df_cor$adj.p<0.05),]
  #
  df_cor$sig='n.s.'
  df_cor$sig[which(  df_cor$gene %in% df_cor_pos$gene )]='positive'
  df_cor$sig[which(  df_cor$gene %in% df_cor_min$gene )]='negative'
  #
  tb=table(df_cor$sig);
  df_cor$sig=factor(df_cor$sig, levels=c('positive','n.s.', 'negative'));
  print(tb)
  
  # save
  #df_cor.high=df_cor[which( (df_cor$cor<0.3 & (df_cor$cor> -0.3)) & df_cor$adj.p<0.01 ),]; dim(df_cor.high)
  write.table(df_cor, paste0(output, keyword, '_01-volcano_Spearman_cor_gDPAU_RNA.df.txt') )
  # df_cor=read.table("cor_BC/MDA-MB-468_01-volcano_Spearman_cor_gDPAU_RNA.df.txt")
  table(df_cor$sig)
  #n.s. negative positive 
  #3181      595      222 
  
  # save
  writeLines( rownames(df_cor[which( df_cor$cor> -0.3 & df_cor$cor< 0.3 & df_cor$adj.p<0.05 ),]), 
              paste0(output, keyword, '_02-Spearman_cor_gDPAU_RNA.low.gene.txt'))
  
  writeLines( rownames(df_cor[which(df_cor$sig=='positive'),]), paste0(output, keyword, '_02-Spearman_cor_gDPAU_RNA.pos_high.gene.txt')) #222
  writeLines( rownames(df_cor[which(df_cor$sig=='negative'),]), paste0(output, keyword, '_02-Spearman_cor_gDPAU_RNA.neg_high.gene.txt')) #595
  writeLines( rownames(df_cor[which( df_cor$cor> -0.02 & df_cor$cor< 0.02 ),]), 
              paste0(output, keyword, '_02-Spearman_cor_gDPAU_RNA.abs_lt0.02.gene.txt'))
  #
  g1=ggplot(df_cor, aes(cor, -log10(p), color=sig))+
    geom_point(size=0.1)+theme_bw()+
    scale_color_manual('Correlation', values=c('red', 'grey80', 'blue'), limits=c('positive', 'n.s.', 'negative') )+
    #scale_color_manual('Correlation', values=c('grey80', 'grey80', 'purple', 'red'))+
    labs(x='Spearman correlation\n(gDPAU and gene expression)', y="-log10(p value of log-rank test)",
         title=keyword);
  library(ggrepel)
  df3=df_cor;
  dd_text=df3[which( df3$adj.p<1e-13 ),];
  print('==> add text, point number: ')
  print( dim(dd_text) )
  g2=g1+geom_text_repel(data=dd_text, aes(cor, -log10(p), label=rownames(dd_text)),
                        color="black",size=3,alpha=0.8)+ 
    guides(colour = guide_legend(override.aes = list(alpha = 1,size=2)))#
  
  CairoPDF(file=paste0(output, keyword, "_01-volcano_Spearman_cor_gDPAU_RNA.pdf"), width=4.3, height=3.5) #width=3.7, height=3
  print(g2)
  dev.off()

  print(paste('Saved: ', output, keyword))
}

save_plot_cor2(df_cor.BC, 'MDA-MB-468', 'cor_BC/')

df_cor.BC['PPP1CA',]
df_cor.BC['UNC50',]
#

