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
dim(rnaM) #18662 225
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

# df_cor.pearson=df_cor
df_cor.BC=get_apa_rna_cor(cid.BC)
df_cor.HeLa=get_apa_rna_cor(cid.HeLa)
df_cor.all=get_apa_rna_cor( c(cid.BC, cid.HeLa) ) #time consuming:1min
#
df_cor.all['CSTF2',]
#df_cor['ITGAE', ]


# check
g1=ggplot(df_cor.HeLa, aes(cor, -log10(p), color=n))+
  geom_point(size=0.2)+theme_bw()+
  #scale_color_manual('Correlation', values=c('blue', 'grey80', 'purple', 'red'))+
  #scale_color_manual('Correlation', values=c('grey80', 'grey80', 'purple', 'red'))+
  labs(x='Spearman correlation between gDPAU and RNA');g1




###################
# process2: plot apa RNA point plot
###################
source("../src/f1_correlaton_gDPAU_RNA.src.R")
save_plot_cor(df_cor.BC, 'BC', 'cor_BC/')
save_plot_cor(df_cor.HeLa, 'HeLa', 'cor_HeLa/')
save_plot_cor( df_cor.all, 'all', 'cor_all/')
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
  # low      mid neg_high pos_high 
  #213     3544       69       46 
  
  #df_cor.high=df_cor[which( (df_cor$cor<0.3 & (df_cor$cor> -0.3)) & df_cor$adj.p<0.01 ),]; dim(df_cor.high)
  write.table(df_cor, paste0(output, keyword, '_01-volcano_Spearman_cor_gDPAU_RNA.df.txt') )
  
  
  #writeLines( rownames(df_cor[which(df_cor$sig=='low'),]), paste0(output, keyword, '_02-Spearman_cor_gDPAU_RNA.low.gene.txt'))
  writeLines( rownames(df_cor[which( df_cor$cor> -0.3 & df_cor$cor< 0.3 & df_cor$adj.p<0.05 ),]), paste0(output, keyword, '_02-Spearman_cor_gDPAU_RNA.low.gene.txt'))
  
  
  writeLines( rownames(df_cor[which(df_cor$sig=='positive'),]), paste0(output, keyword, '_02-Spearman_cor_gDPAU_RNA.pos_high.gene.txt'))
  writeLines( rownames(df_cor[which(df_cor$sig=='negative'),]), paste0(output, keyword, '_02-Spearman_cor_gDPAU_RNA.neg_high.gene.txt'))
  #
  g1=ggplot(df_cor, aes(cor, -log10(p), color=sig))+
    geom_point(size=0.1)+theme_bw()+
    scale_color_manual('Correlation', values=c('red', 'grey80', 'blue'), limits=c('positive', 'n.s.', 'negative') )+
    #scale_color_manual('Correlation', values=c('grey80', 'grey80', 'purple', 'red'))+
    labs(x='Spearman correlation\n(gDPAU and gene expression)', y="-log10(p value of log-rank test)",
         title=keyword);
  library(ggrepel)
  df3=df_cor; 
  #dd_text=df3[which( ((df3$cor >0.8 | df3$cor < (-0.8)) & df3$adj.p<1e-10) | df3$adj.p<1e-18 ),];
  dd_text=df3[which( df3$adj.p<1e-13 ),];
  print('==> add text, point number: ')
  print( dim(dd_text) )
  g2=g1+geom_text_repel(data=dd_text, aes(cor, -log10(p), label=rownames(dd_text)),
                        color="black",size=3,alpha=0.8)+ #xlim(-1,1)+
    guides(colour = guide_legend(override.aes = list(alpha = 1,size=2)))#放大图例的点
  #
  CairoPDF(file=paste0(output, keyword, "_01-volcano_Spearman_cor_gDPAU_RNA.pdf"), width=4.3, height=3.5) #width=3.7, height=3
  print(g2)
  dev.off()
  #scale_color_gradient2(low="navy", mid='lightblue', high='red', midpoint = 120)
  # dim(df_cor.high) #90  5
  #
  print(paste('Saved: ', output, keyword))
}

save_plot_cor2(df_cor.BC, 'MBA-MD-468', 'cor_BC/')

df_cor.BC['PPP1CA',]
df_cor.BC['UNC50',]
#



# delta of cor
head(df_cor.BC)
df_cor.BC_HeLa=(function(){
  gene.common=intersect(rownames(df_cor.BC), rownames(df_cor.HeLa))
  length(gene.common) #3954
  df_cor.BC_HeLa=df_cor.BC[gene.common, ]
  colnames(df_cor.BC_HeLa)=c('gene', 'nBC', 'corBC','pBC', 'adj.pBC')
  df_cor.BC_HeLa=df_cor.BC_HeLa[, -c(4)]
  df_cor.BC_HeLa$nHeLa=df_cor.HeLa[gene.common, ]$n;
  df_cor.BC_HeLa$corHeLa=df_cor.HeLa[gene.common, ]$cor;
  df_cor.BC_HeLa$adj.pHeLa=df_cor.HeLa[gene.common, ]$adj.p;
  #
  delta=c()
  for(gene in gene.common){
    delta=c(delta, df_cor.BC[gene,'cor']-df_cor.HeLa[gene, 'cor'])
  }
  df_cor.BC_HeLa$deltaCor=delta
  df_cor.BC_HeLa=df_cor.BC_HeLa[order(df_cor.BC_HeLa$delta),]
  # filter
  df_cor.BC_HeLa=df_cor.BC_HeLa[which(df_cor.BC_HeLa$adj.pBC<0.05 &df_cor.BC_HeLa$adj.pBC<0.05),]
  df_cor.BC_HeLa=df_cor.BC_HeLa[which(df_cor.BC_HeLa$nBC>=156*0.2 & df_cor.BC_HeLa$nHeLa>=57*0.2),]
  df_cor.BC_HeLa
})()
#
head(df_cor.BC_HeLa)
tail(df_cor.BC_HeLa)
dim(df_cor.BC_HeLa) #681 9

CairoPDF(file="cor_BC/BC03_cor_changes_BC_vs_HeLa.pdf", width=4, height=4)
hist(df_cor.BC_HeLa$delta, n=80, 
     main="delta cor(gDPAU, RNA)\nBC and HeLa", 
     xlab="Correlation of apa and RNA changes\n for each gene in BC and HeLa")
dev.off()
#



# cor changed
# correlation changed gene
rs1=df_cor.BC_HeLa[which(df_cor.BC_HeLa$deltaCor>0.5), ]; rs1
#        gene nBC     corBC         pBC    adj.pBC nHeLa    corHeLa   adj.pHeLa    delta
# ZMIZ1 ZMIZ1  33 0.4854158 0.004190201 0.02573105    24 -0.7179156 0.003979174 1.203331
rs2=df_cor.BC_HeLa[which(df_cor.BC_HeLa$deltaCor< -0.5), ]; rs2
#             gene nBC      corBC          pBC      adj.pBC nHeLa   corHeLa    adj.pHeLa     delta
#HSD17B12 HSD17B12  97 -0.3987445 5.223759e-05 8.328993e-04    30 0.8480446 1.190864e-06 -1.246789
#NT5E         NT5E  78 -0.4998122 3.179597e-06 8.574242e-05    17 0.5352258 1.692939e-01 -1.035038
corChanged=rbind(rs1,rs2)
dim(corChanged) #28 8
write.table(corChanged, "cor_BC/BC04_corChanged.df.txt")
writeLines(corChanged$gene,"cor_BC/BC04_corChanged.gene.txt")

# cor Not changed
corNotChanged=df_cor.BC_HeLa[which( abs(df_cor.BC_HeLa$deltaCor)<0.1), ]
dim(corNotChanged)
write.table(corNotChanged, "cor_BC/BC04_corNotChanged.df.txt")
writeLines(corNotChanged$gene,"cor_BC/BC04_corNotChanged.gene.txt")








##########
# gene list of RNA expresion level: high
gene.high=getRNAHighGene()
head(gene.high)
length(gene.high) #7819


# corPosHigh + rnaHigh
head(df_cor.BC, n=10)
gene.set1=intersect(readLines('cor_BC/BC_02-Spearman_cor_gDPAU_RNA.pos_high.gene.txt'), 
                    gene.high)
length(gene.set1) #76
writeLines(gene.set1, 'cor_BC/BC_03_posHigh_rnaHigh.gene.txt')
head(gene.set1, n=15)
#[1] "CD59"    "NUTF2"   "PTPMT1"  "WDR77"   "SNRPD3"  "RBM34"   "MTHFS"   "SPINT1"  "COPS3"  
#[10] "FBXL5"   "XPNPEP3" "PLAUR"   "TARDBP"  "RPL22"   "OIP5"   


# corNegHigh + rnaHigh
head(df_cor.BC, n=10)
gene.set2=intersect(readLines('cor_BC/BC_02-Spearman_cor_gDPAU_RNA.neg_high.gene.txt'), 
                    gene.high)
length(gene.set2) #135
writeLines(gene.set2, 'cor_BC/BC_03_negHigh_rnaHigh.gene.txt')
head(gene.set2, n=15)
#[1] "BLOC1S5"  "URGCP"    "MRPL43"   "SFT2D1"   "LTV1"     "COG8"     "MAPKAPK5" "TIMM10"  
#[9] "ISCU"     "ASNSD1"   "KIF2A"    "SNAPC5"   "NCBP2"    "ACTR1A"   "NOL12"




# use
###################
# process3: apa RNA point plot
###################
CairoPDF(file="cor_BC/BC03_point_apa_RNA-.pdf", width=3, height=3)
# pos
plot_apa_RNA_by_gene('NUTF2', cid.BC)
plot_apa_RNA_by_gene('CD59', cid.BC, x=70, y=6)
plot_apa_RNA_by_gene('CDKN2B', cid.BC, y=6)
plot_apa_RNA_by_gene('PHF10', cid.BC,y=8)
plot_apa_RNA_by_gene('WDR77', cid.BC,y=8)
#

# neg
plot_apa_RNA_by_gene('BLOC1S5', cid.BC, y=1)
plot_apa_RNA_by_gene('URGCP', cid.BC, x=65, y=2, xmin=55)
plot_apa_RNA_by_gene('KRAS', cid.BC, x=26, y=2.5)
# plot_apa_RNA('BST2')
plot_apa_RNA_by_gene('CDK6', cid.BC,x=26, y=3 )
# plot_apa_RNA('CCNB2',y=9)
plot_apa_RNA_by_gene('MRPL43', cid.BC, y=1)
plot_apa_RNA_by_gene('ISCU', cid.BC, y=1)
#
# low
# plot_apa_RNA('S100A4', xmin=95, x=96)
# plot_apa_RNA('S100A9', xmin=95,x=96)
plot_apa_RNA_by_gene('NOP56', cid.BC, x=45, y=6.5)
plot_apa_RNA_by_gene('GRB2', cid.BC, y=3)

### cor changed gene
plot_apa_RNA_by_gene(rs1$gene[7], cid.BC)
plot_apa_RNA_by_gene(rs1$gene[7], cid.HeLa)
#
plot_apa_RNA_by_gene(rs2$gene[1], cid.BC)
plot_apa_RNA_by_gene(rs2$gene[1], cid.HeLa)

plot_apa_RNA_by_gene('CD47', cid.BC, x=10, y=9 )
plot_apa_RNA_by_gene('CD47', cid.HeLa, x=10, y=1 )
#
dev.off()
#









###################
# process4: cor(apa, RNA) low, and apa changed gene across cell cycle
# (其apa与RNA表达量不相关，且apa变动的基因)
###################
#(1) RNA, in cycle
gene.BC.setRNA=(function(){
  setAll=c()
  for(phase in c('G1S', 'S', 'G2M', 'M', 'MG1')){
    set1=readLines( paste0('../../f3/cell_cycle/result/BC_CycleRelatedGene_RNA_',phase,'.txt') )
    setAll=c(setAll, set1)
  }
  return( unique(setAll) )
})()
pL(gene.BC.setRNA, 'setRNA') #290
writeLines( gene.BC.setRNA, 'geneSet_BC/BC00_RNA_inCycle.gene.txt')
#(2) apa, in cycle
gene.BC.setAPA=readLines('../../f4/cycles/gDPAU_BC/BC07_cycleClusterAll.gene.txt')
pL(gene.BC.setAPA, 'setAPA') #812
writeLines( gene.BC.setAPA, 'geneSet_BC/BC00_APA_inCycle.gene.txt')

#(3) low Cor(apa, RNA)
gene.BC.setLowCor=readLines('cor_BC/BC_02-Spearman_cor_gDPAU_RNA.low.gene.txt')
pL(gene.BC.setLowCor, 'setLowCor') #212


(function(){
  ######## both change in RNA & APA;
  gene.BC.setRNA_APA=intersect(gene.BC.setRNA, gene.BC.setAPA)
  pL(gene.BC.setRNA_APA, 'setRNA_APA') #39
  writeLines( gene.BC.setRNA_APA, 'geneSet_BC/BC01_RNA_APA.gene.txt')
  
  ## only change in APA, not in RNA
  gene.BC.setAPA_notRNA=setdiff(gene.BC.setAPA, gene.BC.setRNA)
  pL(gene.BC.setAPA_notRNA, 'setRNA_APA') #773
  writeLines( gene.BC.setAPA_notRNA, 'geneSet_BC/BC02_APA_notRNA.gene.txt')
})()
