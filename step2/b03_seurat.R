#01_seurat.R

## Aim: determine BC and HeLa

pP=function(...){print(paste(...))}
pP0=function(...){print(paste0(...))}
pd=function(...){print(dim(...))}
pl=function(...){print(length(...))}
ph=function(...){print(head(...))}
pt=function(...){print(tail(...))}
debug=function(...){
  pd(...)
  ph(...)
}

#
setwd("/home/wangjl/data/apa/20200701Fig/f2/BC_HeLa/")
getwd()
##
library(Matrix)
library(Seurat)
library(dplyr)
library(Cairo)



############
# Load dataset
loadData=function(){
  getRNA=function(filePath, cellType){
    pP(filePath, cellType)
    rna <- read.csv(filePath,row.names = 1)
    colnames(rna)=gsub("_","", colnames(rna))
    cellInfo=data.frame(
      cid=colnames(rna),
      cellType0=cellType
    )
    pd(rna)
    return(list(rna, cellInfo))
  }
  #
  e15=getRNA("/home/wangjl/data/apa/190530Mix/expression_HTseq.15Mixed.csv", "HeLa_syncMix")
  e13=getRNA("/home/wangjl/data/apa/190528L/expression_HTseq.13L.csv", "HeLa_sync")
  e207=getRNA("/home/wangjl/data/apa/190517R/expression_HTseq.207R.csv", "unknown")
  # check gene symbol
  t1= identical(rownames(e15),rownames(e13)) #T
  t2= identical(rownames(e15),rownames(e207)) #T
  pP(t1, t2)
  #
  rnaM=cbind(e15[[1]], e13[[1]], e207[[1]])
  cellInfo=rbind(e15[[2]], e13[[2]], e207[[2]])
  #
  rownames(cellInfo)=cellInfo$cid;  
  return(list(rnaM, cellInfo))
}
rs=loadData()
rnaM=rs[[1]]
cellInfo=rs[[2]]

#check
dim(rnaM) #37143   235
rnaM[1:5,1:10]
dim(cellInfo) #235   2
cellInfo[1:5,]
#
write.csv( rnaM, file="0_expression_HTseq-15mix_13syncHeLa_207Right-withPreColnames.235.csv" )


###########
# filter by counts and gene: median+-3*mad
filterByGeneAndCounts=function(df){
  pP('init cell number:',ncol(df) ) #37143   235
  # by counts
  rs=as.data.frame( apply(df, 2, sum) )
  colnames(rs)=c("counts")
  rs$cid=rownames(rs)
  #
  th1=median(rs$counts)-3*mad(rs$counts)
  th2=median(rs$counts)+3*mad(rs$counts)
  #
  rs=rs[which(rs$counts>=th1 & rs$counts<=th2),] #filter out cells
  pP("filter out cell number by total counts(median+-3*mad):", nrow(rs) ) #230 2
  df=df[,rownames(rs)]
  #
  #by gene number
  rs=as.data.frame( apply(df>0, 2, sum) );
  colnames(rs)=c("geneNumber")
  rs$cid=rownames(rs)
  #
  th1=median(rs$geneNumber)-3*mad(rs$geneNumber)
  th2=median(rs$geneNumber)+3*mad(rs$geneNumber)
  #
  rs=rs[which(rs$geneNumber>=th1 & rs$geneNumber<=th2),] #filter out cells
  pP("filter out cell number by gene number(median+-3*mad):", nrow(rs) ) #230 2
  df=df[,rownames(rs)]
  #
  df
}
rs=filterByGeneAndCounts(rnaM)
head(rs[,1:10])
dim(rs)#
rnaM2=rs

'c15ROW02' %in% colnames(rnaM2)
'c15ROW02' %in% colnames(rnaM) 

#
#########
#from rnaM2
pbmc <- CreateSeuratObject(counts = rnaM2, min.cells = 2, min.features = 200, 
                           project = "BC_HeLa")
pbmc #18710 features across 230 samples within 1 assay
#

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc@meta.data, 5)

# Visualize QC metrics as a violin plot
CairoPDF(file="03-nFeature-nCount-percent.MT.pdf",width=10,height=5)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


###########
#we visualize QC metrics, and use these to filter cells.
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CairoPDF(file="04-percent.MT-nFeature.pdf",width=10,height=5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

# filter out cells by mt pct
pbmc2 <- subset(pbmc, subset = nFeature_RNA > 2000 & nFeature_RNA < 10000 & percent.mt < 40)
pbmc2 #18542 features across 225 samples within 1 assay  

#filter out who?
rs=pbmc@meta.data[pbmc@meta.data$percent.mt>=40,];rs
#         orig.ident nCount_RNA nFeature_RNA percent.mt
#c13ROW17    BC_HeLa     264633         2592   47.13207
#c13ROW18    BC_HeLa     182124         2544   42.88726
#c14ROW33    BC_HeLa      58749         2239   77.47877
#c15ROW10    BC_HeLa     976401         5170   62.41155
#c16ROW06    BC_HeLa     262003         4038   51.85971


########
# save RNA count data
rnaM3=rnaM2[, setdiff( colnames(rnaM2)  , rownames(rs) ) ]
dim(rnaM3)
write.csv(rnaM3,file="1_BC_HeLa.225cells.count.csv")


####
#
#
# Normalizing the data: log( expr/totalExpr * 10000)
pbmc3 <- NormalizeData(pbmc2, normalization.method = "LogNormalize", scale.factor = 1e6)
#or  pbmc3 <- NormalizeData(pbmc2)


######## HVG
# Identification of highly variable features (feature selection)
pbmc4 <- FindVariableFeatures(pbmc3, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc4), 10);top10
#[1] "NEAT1"  "COL1A1" "ACTA2"  "MGP"    "S100A7" "LGALS1" "SPINK6" "C4BPB"  "DDIT3"  "LAMA1" 

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc4)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE,xnudge=0,ynudge=0)+labs('BC_HeLa')
#CombinePlots(plots = list(plot1, plot2)) #error
CairoPDF(file="05-FindVariableGenes.pdf",width=7,height=5)
#CombinePlots(plots = list(plot1, plot2))
plot2
dev.off()


# save HVG
x = pbmc4@assays$RNA@var.features
length(x) # 2000
x[1:10]
write.csv(x,"05_BC_HeLa.225cells.variableGenes.csv")




#
#############
#Scaling the data
# rm(br6,br7,br8,br9,rs,rs2,tmp,top10,top5)
# rm(br.2deg, br.2marker, dense.size, soarse,size)

# The results of this are stored in pbmc[["RNA"]]@scale.data
all.genes <- rownames(pbmc4)
pbmc5 <- ScaleData(pbmc4, features = all.genes)






####### PCA
#Perform linear dimensional reduction
pbmc6 <- RunPCA(pbmc5, features = VariableFeatures(object = pbmc5),
                nfeatures.print = 15,
                npcs=60)
#
#Seurat provides several useful ways of visualizing both cells and features that define the PCA, 
#including VizDimReduction,  DimPlot, and DimHeatmap
#1. Examine and visualize PCA results a few different ways
print(pbmc6[["pca"]], dims = 1:3, nfeatures = 50)
#2
CairoPDF(file="06-VizPCA.pdf",width=8,height=8)
VizDimLoadings(pbmc6, dims = 1:2, reduction = "pca",nfeatures = 50)
dev.off()
#3
CairoPDF(file="07-PCAPlot.pdf",width=15,height=5)
p1=DimPlot(pbmc6, reduction = "pca", dims=c(1,2))
p2=DimPlot(pbmc6, reduction = "pca", dims=c(1,3))
p3=DimPlot(pbmc6, reduction = "pca", dims=c(2,3))
CombinePlots(list(p1,p2,p3), ncol=3)
dev.off()







############
# color sync HeLa in PCA
############
tmp_pc=as.data.frame(pbmc6@reductions$pca@cell.embeddings)
dim(tmp_pc) #[1] 225  60
head(tmp_pc[,1:5])
#
head(cellInfo)
dim(cellInfo[rownames(tmp_pc),])
#
tmp_pc$cellType=cellInfo[rownames(tmp_pc),]$cellType0
head(tmp_pc[,c(1:3, ncol(tmp_pc))])
#
table(tmp_pc$cellType)
write.table(tmp_pc, '07-PCAPlot-HeLa_BC_255cells.txt')

#define function
PCPlot=function(pcX,pcY){
  ggplot(tmp_pc, aes(tmp_pc[,pcX], tmp_pc[,pcY], col=factor(cellType) ))+
    geom_point(size=0.5)+
    labs(title="", x=pcX, y=pcY)+
    guides(color=guide_legend(title=NULL)) +
    theme_bw()+
    scale_color_manual(values=c('black',"#F8766D",'grey'))
}
#
Do_PCPlot=function(){
  p1=PCPlot("PC_1",'PC_2')
  p2=PCPlot("PC_1",'PC_3')
  p3=PCPlot("PC_2",'PC_3')
  #
  CombinePlots(plots = list(p1, p2, p3),ncol=3)
}
CairoPDF(file=paste0("07-PCAPlot-syncHeLa-syncHeLa_inRed.pdf"),width=15,height=4)
Do_PCPlot()
dev.off()





############
CairoPDF(file="08-PCHeatmap.pdf",width=8,height=8)
DimHeatmap(pbmc6, dims = 1, cells = 225, balanced = TRUE)
#for each PC
DimHeatmap(pbmc6, dims = 1:12, cells = 225, balanced = TRUE)
dev.off()



######### Select PC
#Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
pbmc7 <- JackStraw(pbmc6, num.replicate = 100,dims =50)
pbmc8 <- ScoreJackStraw(pbmc7, dims = 1:50)
#
CairoPDF(file="09-JackStraw.pdf",width=7,height=5)
JackStrawPlot(pbmc8, dims = 1:20)
dev.off()

# br8=br7
CairoPDF(file="10-ElbowPlot.pdf",width=5,height=5)
ElbowPlot(pbmc8, ndims=30 ) #+ theme_gray() #类似碎石图
dev.off()



######### tSNE, and UMAP
pbmc9 <- FindNeighbors(pbmc8, k.param=20, dims = 1:8)
pbmc10 <- FindClusters(pbmc9, resolution = 0.6,algorithm=1) #resolution between 0.4-1.2 

head(Idents(pbmc10), 15) #cluster and cid
table(Idents(pbmc10))
#0  1  2 ###
#96 72 57 when 1e6, now I use;
#87 82 56 when 1e4

pbmc11 <- RunUMAP(pbmc10, dims = 1:8)
pbmc11 <- RunTSNE(pbmc11, dims = 1:8)
#
CairoPDF(file="11-DimPlot-UMAP-tSNE.pdf",width=4,height=3.5)
DimPlot(pbmc11, reduction = "umap")
DimPlot(pbmc11, reduction = "tsne")
dev.off()

#save to file
saveRDS(pbmc11, file = "BC_HeLa.225cells_umap_tSNE-pbmc11_version5.rds")

table(pbmc11@active.ident)
#0  1  2 
#93 75 57
#




############
# get cell id and cluster number
head(pbmc11@reductions$umap@cell.embeddings)
tmp_df=pbmc11@reductions$umap@cell.embeddings
dim(tmp_df[tmp_df[,1]<0,]) #59 2
## get cluster number
tmp_df2=data.frame(
  cellid=names(Idents(pbmc10)),
  cluster=unname( Idents(pbmc10) )
)
rownames(tmp_df2)=tmp_df2$cellid
dim(tmp_df2)
head(tmp_df2)
cellInfo=cellInfo[rownames(tmp_df2),]
cellInfo$cluster=tmp_df2$cluster
head(cellInfo)
table(cellInfo$cellType0, cellInfo$cluster)
#              0  1  2
#HeLa_sync     0  0 12
#HeLa_syncMix  0  0 15
#unknown      93 75 30
#
cellInfo$cellType='BC_0';
cellInfo[which(cellInfo$cluster==1),]$cellType="BC_1"
cellInfo[which(cellInfo$cluster==2 & cellInfo$cellType0=='unknown'),]$cellType='HeLa_normal'
cellInfo[which(cellInfo$cluster==2 & cellInfo$cellType0 %in% c('HeLa_sync','HeLa_syncMix') ),]$cellType='HeLa_sync'
#
cellInfo$cellType2=cellInfo$cellType
cellInfo[which( cellInfo$cellType0 == 'HeLa_syncMix' ),]$cellType2='HeLa_syncMix'
#
table(cellInfo$cellType)
table(cellInfo$cellType2)
table(cellInfo$cellType, cellInfo$cellType2)
#
head(cellInfo)

#
#
# Error cluster cells: in cluster0 or 1 area but UMAP_1>5
rs0=intersect(row.names(tmp_df[tmp_df[,1]>5,]), rownames(cellInfo[which(cellInfo$cluster==0),]));
rs0 #"c13ROW23"
rs1=intersect(row.names(tmp_df[tmp_df[,1]>5,]), rownames(cellInfo[which(cellInfo$cluster==1),]));
rs1 #"c14ROW23" "c15ROW23"
#
cellInfo2=cellInfo[!(rownames(cellInfo) %in% c(rs0,rs1) ) ,]
dim(cellInfo2) #222 5
head(cellInfo2)
#
table(cellInfo2$cellType)
table(cellInfo2$cellType2)
table(cellInfo2$cellType, cellInfo$cellType2)
write.table(cellInfo2, '11_cellInfo.v5.0_222cells.txt')
#
cellInfo2=read.table('11_cellInfo.v5.0_222cells.txt')
dim(cellInfo2)
head(cellInfo2)

#
pbmc12=subset(pbmc11, cells=rownames(cellInfo2))
pbmc12 #18710 features across 222 samples within 1 assay


############
# plot sync HeLa in UMAP
#save to file
saveRDS(pbmc12, file = "BC_HeLa.222cells_umap_tSNE-pbmc12_forUMAP_diy.rds")
# pbmc12=readRDS("BC_HeLa.222cells_umap_tSNE-pbmc12_forUMAP_diy.rds")

cellInfo2=read.table('/data/jinwf/wangjl/apa/20200701Fig/f3/cell_cycle/cellInfo.V6.txt', row.names = 1)
head(cellInfo2)



##################
# add cellType info to seurat object, instead of cluster 0, 1, 2
##################
# add new column as cell identity
pbmc12@meta.data$cellType=cellInfo2[rownames( pbmc12@meta.data), ]$cellType
table(pbmc12@meta.data$seurat_clusters, pbmc12$cellType)
#   BC_0 BC_1 HeLa_normal HeLa_sync
#0   92    0           0         0
#1    0   73           0         0
#2    0    0          30        27
Idents(pbmc12) ="cellType"
sortIdent = factor(pbmc12$cellType,levels = c("BC_0","BC_1", "HeLa_normal",   "HeLa_sync"))
pbmc12@active.ident = sortIdent
levels(pbmc12@active.ident)
#[1] "BC_0"        "BC_1"        "HeLa_normal" "HeLa_sync"
#save to file
saveRDS(pbmc12, file = "BC_HeLa.222cells_umap_tSNE-pbmc12.1_forUMAP_diy_addCellType.rds")
# pbmc12=readRDS("BC_HeLa.222cells_umap_tSNE-pbmc12.1_forUMAP_diy_addCellType.rds")



############# tSNE and UMAP plot
write.csv(pbmc12@reductions$umap@cell.embeddings, '11-data-UMAP12-222cells.csv')
write.csv(pbmc12@reductions$tsne@cell.embeddings, '11-data-tsne12-222cells.csv')
write.csv(pbmc12@reductions$pca@cell.embeddings, '11-data-pca60-222cells.csv')
# test
(function(){
  gene='CCNE2'
  df1=as.data.frame(pbmc12@reductions$pca@cell.embeddings)
  #ggplot(df1, aes(PC_1, PC_2))+geom_point()
  ggplot(df1, aes(PC_1, PC_2, color=rnaM.log2cpm[gene, rownames(df1)]))+geom_point()+
    scale_color_gradient2(gene,low = 'navy', mid='white', high='red', midpoint = 5)+
    theme_bw()
})()

#
df1=(function(){
  tmp_df=pbmc12@reductions$umap@cell.embeddings
  #tmp_df=pbmc12@reductions$tsne@cell.embeddings
  tmp_df=as.data.frame(tmp_df)
  pd(tmp_df) #222 2
  ph(tmp_df)
  #
  #set
  tmp_df$type=cellInfo2[rownames(tmp_df),'cellType0']
  print( table(tmp_df$type) )
  # HeLa_sync HeLa_syncMix      unknown 
  #12           15          195
  # set 2: more details
  tmp_df$type2=cellInfo2[rownames(tmp_df),'cellType']
  tmp_df$type2=factor(tmp_df$type2, levels=c('BC_0', 'BC_1', 'HeLa_normal','HeLa_sync'))
  
  #get
  print(tmp_df[1:5,])
  
  #CairoPDF(file="11-UMAPPlot-syncHeLaTag.pdf",width=4,height=2.8)
  CairoPDF(file="11-tSNEPlot-syncHeLaTag.pdf",width=4,height=2.8)
  g=ggplot(tmp_df, aes(UMAP_1, UMAP_2, color=factor(type) ))+
  #g=ggplot(tmp_df, aes(tSNE_1, tSNE_2, color=factor(type) ))+
    geom_point(size=0.2)+
    guides(color=guide_legend(title=NULL)) +
    theme_bw()+
    scale_color_manual(values=c('#98BEFD',"#D1E2FF",'grey'))
    #scale_color_manual(values=c('#619CFF',"#98BEFD",'grey'))
  print(g)
  #
  g2=ggplot(tmp_df, aes(UMAP_1, UMAP_2, color=type2 ))+
  #g2=ggplot(tmp_df, aes(tSNE_1, tSNE_2, color=type2 ))+
    geom_point(size=0.5)+
    guides(color=guide_legend(title=NULL)) +
    theme_bw()+
    scale_color_manual('Cell type',values=c("#FF9ECE","#F81082", '#005FFF', '#98BEFD'))
  print(g2)
  dev.off()
  tmp_df
})()
dim(df1);head(df1)



write.table(df1, '11-UMAPPlot-syncHeLaTag_TextOnPic.df.txt')
# tag on figure
df1$type2=as.character(df1$type2)
df1[which( substring(df1$type2,1,2)=="BC" ),]$type2="MDA-MB-468"
df1$type2=factor(df1$type2, levels=c("MDA-MB-468", "HeLa_normal", "HeLa_sync" ))
levels(df1$type2)
#
CairoPDF(file="11-UMAPPlot-syncHeLaTag_TextOnPic-2.pdf",width=3,height=3)
g3=ggplot(df1, aes(UMAP_1, UMAP_2, color=type2 ))+
  #g2=ggplot(tmp_df, aes(tSNE_1, tSNE_2, color=type2 ))+
  geom_point(size=0.5)+
  guides(color=guide_legend(title=NULL)) +
  annotate("text", x=-2, y=-1, label='MDA-MB-468', size=4, color="#F81082")+
  #annotate("text", x=1.9, y=1, label='BC_0', size=4, color="#FF9ECE")+
  #annotate("text", x=-5.8, y=-0.8, label='BC_1', size=4, color="#F81082")+
  #
  annotate("text", x=4.8, y=-4.9, label='HeLa_normal', size=4, color="#005FFF")+
  annotate("text", x=4.5, y=-3, label='HeLa_sync', size=4, color="#98BEFD")+
  theme_classic()+
  # scale_color_manual('Cell type',values=c("#FF9ECE", "#F81082", '#005FFF', '#98BEFD'))+
  scale_color_manual('Cell type',values=c( "#F81082", '#005FFF', '#98BEFD'))+
  theme(
    legend.position = "none"
  )
print(g3)
dev.off()















#############
# find marker, heatmap
## add claster for HeLa_sync
head( pbmc12@meta.data )

# find all markers distinguishing cluster 5 from clusters 0,1 and 2

# BC vs HeLa
BC_HeLa.markers <- FindMarkers(pbmc12, ident.1 = c("BC_0","BC_1"), ident.2 = c("HeLa_normal",   "HeLa_sync"), min.pct = 0.25)
BC_HeLa.markers$gene=rownames(BC_HeLa.markers)
head(BC_HeLa.markers, n = 5)
BC_HeLa.markers['KRT19',]
# 
rs1=BC_HeLa.markers %>% filter(pct.1>0.7 | pct.2>0.7 ) %>% filter(p_val_adj<0.05 & abs(avg_logFC)>1); dim(rs1)
writeLines((rs1[which(rs1$avg_logFC>1),]$gene), '14_BC_vs_HeLa_plus.gene')
writeLines((rs1[which(rs1$avg_logFC<1),]$gene), '14_BC_vs_HeLa_minus.gene')
hist(rs1$pct.1, n=100)
hist(rs1$pct.2, n=100)
#
drawHM=function(rs){
  rs=rs %>% filter(p_val_adj<1e-5 )
  gset1= rs %>% top_n(n = 25, wt = avg_logFC) #only in BC
  gset2= rs %>% top_n(n = 25, wt = -avg_logFC) #only in HeLa
  
  rs_1=rs[order( -rs$avg_logFC ),]
  # top 100 and top 20 pct
  t1=head(rs_1,n=100)
  t1=t1[order(-t1$pct.1), ]
  #
  t2=tail(rs_1,n=100)
  t2=t2[order(t2$pct.2), ]
  n=25
  rs_2=rbind(t1[1:n, ], tail(t2,n=n) )
  
  p=DoHeatmap(pbmc12, 
            # features = rbind(gset1, gset2)$gene, 
            features = rs_2$gene, 
            group.colors= c("#FF9ECE","#F81082", '#005FFF', '#98BEFD') )
  print(p)
  rs_2$gene;
};
CairoPDF(file="14-heatmap_BC_HeLa.pdf",width=5,height=8)
mk.gene.BC_HeLa=drawHM(rs1)
dev.off()



# BC0 vs BC1
BC0_BC1.markers <- FindMarkers(pbmc12, ident.1 = "BC_0", ident.2 ="BC_1", test.use = 'wilcox',min.pct = 0.25)
BC0_BC1.markers$gene=rownames(BC0_BC1.markers)
head(BC0_BC1.markers)
rs1.bc0_bc1=BC0_BC1.markers %>% filter(pct.1>0.6 | pct.2>0.6) %>% filter(p_val_adj<0.05 & abs(avg_logFC)>log2(1.5)); debug(rs1.bc0_bc1)
hist(rs1.bc0_bc1$pct.2, n=100)
#
FeaturePlot(pbmc12,pt.size = 0.3, features = rs1.bc0_bc1$gene[1:16] )
FeaturePlot(pbmc12,pt.size = 0.3, features = rev(rs1.bc0_bc1$gene)[1:16] )
DoHeatmap(pbmc12, features = rs1.bc0_bc1[order(-rs1.bc0_bc1$pct.2),]$gene, 
          group.colors= c("#FF9ECE","#F81082", '#005FFF', '#98BEFD') )
drawHM(rs1.bc0_bc1)
mk.gene.BC0_BC1=rs1.bc0_bc1[order(-rs1.bc0_bc1$pct.2),]$gene
length(mk.gene.BC0_BC1)
#[1] 104
#

# BC0,BC1 vs normal HeLa
BC_HeLa_normal.markers <- FindMarkers(pbmc12, ident.1 = c("BC_0","BC_1"), ident.2 ="HeLa_normal", test.use = 'wilcox',min.pct = 0.25)
BC_HeLa_normal.markers$gene=rownames(BC_HeLa_normal.markers)
head(BC_HeLa_normal.markers)
dim(BC_HeLa_normal.markers) #6704
#
rs1.bc_helaN=BC_HeLa_normal.markers %>% filter(pct.1>0.5 | pct.2>0.5) %>% filter(p_val_adj<0.05 & abs(avg_logFC)>log2(1.5)); debug(rs1.bc_helaN)
head(rs1.bc_helaN[order(rs1.bc_helaN$avg_logFC), ])
head(rs1.bc_helaN[order(-rs1.bc_helaN$avg_logFC), ])
hist(rs1.bc_helaN$pct.2, n=100)
#
FeaturePlot(pbmc12,pt.size = 0.3, features = rs1.bc_helaN[order(-rs1.bc_helaN$avg_logFC), ]$gene[1:16] ) #BC high
FeaturePlot(pbmc12,pt.size = 0.3, features = rs1.bc_helaN[order(rs1.bc_helaN$avg_logFC), ]$gene[1:16] ) #HeLa high
#

# BC1 vs all
bc1.marker=(function(){
  tmp <- FindMarkers(pbmc12, ident.1 = c("BC_1"),only.pos = F, ident.2 =c('BC_0',"HeLa_normal", 'HeLa_sync'), 
                    test.use = 'wilcox',min.pct = 0.25)
  tmp$gene=rownames(tmp)
  tmp1=tmp %>% filter(pct.1>0.5 | pct.2>0.5) %>% filter(p_val_adj<0.05 & abs(avg_logFC)>log2(1.5));
  tmp1=tmp1[order(-tmp1$avg_logFC),]
  tmp1
})()
head(bc1.marker)
tail(bc1.marker)
FeaturePlot(pbmc12,pt.size = 0.3, features = bc1.marker$gene[1:16] )


# BC1 vs BC0
bc1.markerI=(function(){
  tmp <- FindMarkers(pbmc12, ident.1 = c("BC_1"),only.pos = F, ident.2 =c('BC_0'), 
                     test.use = 'wilcox',min.pct = 0.25)
  tmp$gene=rownames(tmp)
  tmp1=tmp %>% filter(pct.1>0.5 | pct.2>0.5) %>% filter(p_val_adj<0.05 & abs(avg_logFC)>log2(1.5));
  tmp1=tmp1[order(-tmp1$avg_logFC),]
  tmp1
})()
head(bc1.markerI)
tail(bc1.markerI)
FeaturePlot(pbmc12, pt.size = 0.3, features = bc1.markerI$gene[1:16] )
FeaturePlot(pbmc12, pt.size = 0.3, features = rev(bc1.markerI$gene)[1:16] )
#

# HeLa normal vs sync
helaN.markerI=(function(){
  tmp <- FindMarkers(pbmc12, ident.1 = c("HeLa_normal"),only.pos = F, ident.2 =c('HeLa_sync'), 
                     test.use = 'wilcox',min.pct = 0.25)
  tmp$gene=rownames(tmp)
  tmp1=tmp %>% filter(pct.1>0.5 | pct.2>0.5) %>% filter(p_val_adj<0.05 & abs(avg_logFC)>log2(1.5));
  tmp1=tmp1[order(-tmp1$avg_logFC),]
  tmp1
})()
head(helaN.markerI)
tail(helaN.markerI)
FeaturePlot(pbmc12, pt.size = 0.3, features = helaN.markerI$gene[1:16] )
FeaturePlot(pbmc12, pt.size = 0.3, features = rev(helaN.markerI$gene)[1:16] )
#

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc12, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
rs0=pbmc.markers %>% group_by(cluster) %>% filter(p_val_adj<0.05 & abs(avg_logFC)>log2(2) ) 
#rs=rbind( rs0 %>% top_n(n = 10, wt = avg_logFC ), rs0 %>% top_n(n = 10, wt = -avg_logFC ))
rs= rs0 %>% top_n(n = 10, wt = avg_logFC )
rs=rs[!duplicated(rs$gene),]
rs=rs[order(rs$cluster, -rs$avg_logFC),]
rs
#    p_val avg_logFC pct.1 pct.2   p_val_adj cluster  gene   
#    <dbl>     <dbl> <dbl> <dbl>       <dbl> <fct>    <chr>  
#1 2.12e-11     -1.02 0.75  0.977  3.96e- 7 BC_0        LIMCH1 
#2 1.82e- 7      1.04 0.989 0.938  3.40e- 3 BC_0        KRT17  
#3 2.12e- 7     -1.00 0.598 0.823  3.96e- 3 BC_0        GFPT1  
#4 1.48e- 7      2.25 0.384 0.094  2.77e- 3 BC_1        ZNF551 
#5 2.19e- 7      2.52 0.411 0.128  4.10e- 3 BC_1        ZNF761 
#6 2.26e- 7      2.05 0.384 0.094  4.22e- 3 BC_1        TTN-AS1
#7 6.23e- 9      3.05 0.7   0.297  1.17e- 4 HeLa_normal SDR16C5
#8 1.72e- 8      2.77 0.467 0.089  3.23e- 4 HeLa_normal CDH10  
#9 2.65e- 6      2.93 0.4   0.104  4.96e- 2 HeLa_normal SACS   
#10 4.00e-20      5.09 1     0.349  7.49e-16 HeLa_sync   ACTA2  
#11 3.17e-18      4.38 1     0.436  5.93e-14 HeLa_sync   RPS4XP2
#12 1.05e-16      4.40 1     0.636  1.97e-12 HeLa_sync   COL1A1 

DoHeatmap(pbmc12, features = c( rs$gene ), 
          group.colors= c("#FF9ECE","#F81082", '#005FFF', '#98BEFD') )
#
CairoPDF(file="14-heatmap_all4cluster-raw.data.pdf",width=5,height=8)
features=unique(c(
  rs1.bc_helaN[order(-rs1.bc_helaN$avg_logFC), ]$gene[1:16], #BC high
  mk.gene.BC_HeLa[1:8], 'KRT17', 'S100A14', #BC high
                 rs$gene[11:20], # BC1 high
                 rs$gene[1:10], #BC0 low, HeLa high
                 'BST2','PDE4D','TGFBR1','CD55','EMP3', #'hela high'
                 rev(mk.gene.BC_HeLa)[1:5],  #HeLa high
  rs1.bc_helaN[order(rs1.bc_helaN$avg_logFC), ]$gene[1:16], #HeLa high
                 rs$gene[21:40] ) )
features=features[-grep('\\-', features) ]
DoHeatmap(pbmc12, features = features, #slot='counts',
          group.colors= c("#FF9ECE","#F81082", '#005FFF', '#98BEFD') )
dev.off()


write.table(pbmc12@assays$RNA@scale.data, '14-heatmap_all4cluster.scale.data.txt')
write.table(pbmc12@assays$RNA@data, '14-heatmap_all4cluster.data.txt')










# main figure: heatmap
(function(){
  #dt=read.table('14-heatmap_all4cluster.data.txt', row.names = 1)
  dt=pbmc12@assays$RNA@data
  pd(dt)
  ph(dt[,1:5])
  #
  # add cell info
  annotation_col = data.frame(
    cellType = pbmc12@meta.data$cellType,
    id=1
  )
  rownames(annotation_col) = colnames(dt)
  annotation_col=annotation_col[order(annotation_col$cellType),]
  #
  cellInfo=cellInfo[rownames(annotation_col),]
  annotation_col2=data.frame(
    #cellType2=cellInfo$cellType2,
    cellType=annotation_col$cellType
    )
  rownames(annotation_col2) = rownames(annotation_col)
  
  genelist.hm=c(
    #BC
    "EGFR",  "PI3", "MGP",  "EPCAM",  "KRT19",  "RIPK4", 'MMP7', 'S100A9',
    "KRT17",  "LINC00665", 
             
             #HeLa
             "TGFBR1", "VIM",
             "CD55", 
             "EMP3", 
             "BST2", 
             "C4BPB", 
             "PDE4D", 
             "BASP1", 
             "CTAG2", 
             "TRABD2A", 
             "CPS1", 
             "PRSS21", "FSTL3", 
             
             #BC1
             "EXPH5", "LYPD6B",
             "ZNF551", 
             "SLFN11", 
             "LEAP2", 
             "FBXO33", #"IMPAD1", 
             'ZNF781','APOBEC3F','ZNF800',
             # BC1
             "WDFY3",
             "KAT6A",
             "SECISBP2L",
             "ANKFY1",
             "ZNF800",
             "PHLDB2",

             #hela normal
            
             
             # HeLa sync
             "ACTA2", 
             "COL1A1", "COX5BP1" ,"COL1A2", 
    
             "RPS4XP2", 
             "RPS2P28", 
             
             "PRRX1", 
             "PLEKHA5", 
             "COL3A1", 
             "ACTG2",
    
     "RPS2P28" , "COL1A1" ,  "PRRX1",    "COL3A1" ,  
     "EIF1P5" ,  "CYCSP1"  ,"SPARC" ,   "RPS3AP12", "RPS7P10",  "COL5A2",   "ZFHX4"
  )
  genelist.hm=unique(genelist.hm)
  
  dt=dt[genelist.hm, rownames(annotation_col)]
  writeLines(rownames(dt), 'tmp.txt')
  write.table(dt, '14-heatmap_all4cluster.drawHM.data.txt')
  #
  # define colors
  ann_colors = list(
    cellType = c('BC_0'="#FF9ECE", 'BC_1'="#F81082", 'HeLa_normal'="#005FFF", 'HeLa_sync'="#98BEFD") #D1E2FF
    #cellType2 = c('BC_0'="#FF9ECE", 'BC_1'="#F81082", 'HeLa_normal'="#005FFF", 'HeLa_sync'="#98BEFD", 'HeLa_syncMix'='#D1E2FF') #
  )
  #head(ann_colors)
  #gap
  gaps=cumsum(as.numeric( table(annotation_col2$cellType) ) )
  bk <- seq(-1.5,1.5,0.001) #c(seq(-3,-0.1,by=0.01),seq(0,3,by=0.01))
  CairoPDF(file="14-heatmap_all4cluster-data-remake-small.pdf",width=5, height=6.5)
  pheatmap(dt, scale='row', 
           clustering_method = 'ward.D2',
           border_color = NA,
           show_colnames = F, 
           cluster_cols = F, cluster_rows = F,
           annotation_col = annotation_col2,
           annotation_colors = ann_colors,
           gaps_col = gaps,
           #gaps_row = c(8, 21, 31),
           #color
           color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
           legend_breaks=seq(-3,3,1), fontsize_row = 8,
           breaks=bk)
  dev.off()
  print('end')
})()
#

#
#
pbmc.markers.roc=FindAllMarkers(pbmc12, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
dim(pbmc.markers.roc) #863 7
head(pbmc.markers.roc)
#       myAUC  avg_diff power pct.1 pct.2 cluster   gene
#RPL28  0.782 0.5342063 0.564     1 1.000       0  RPL28
#RPS11  0.769 0.2944425 0.538     1 1.000       0  RPS11
#KRT19  0.759 0.6678969 0.518     1 0.838       0  KRT19
#CHCHD2 0.752 0.4275734 0.504     1 0.992       0 CHCHD2
#MT-TT  0.751 0.6135705 0.502     1 0.877       0  MT-TT
#S100A9 0.750 0.8479954 0.500     1 0.808       0 S100A9


VlnPlot(pbmc12, pt.size = 0.3, features ='KRT19' )
FeaturePlot(pbmc12,pt.size = 0.3, features = 'KRT19' )
FeaturePlot(pbmc12,pt.size = 0.3, features = 'KRT17' )
#
FeaturePlot(pbmc12,pt.size = 0.3, features = 'EGFR' )
FeaturePlot(pbmc12,pt.size = 0.3, features = 'BST2' )
FeaturePlot(pbmc12,pt.size = 0.3, features = 'ZNF761' )
FeaturePlot(pbmc12,pt.size = 0.3, features = 'LEAP2' )
FeaturePlot(pbmc12,pt.size = 0.3, features = 'TGFBR1' )
FeaturePlot(pbmc12,pt.size = 0.3, features = 'VIM' )
FeaturePlot(pbmc12,pt.size = 0.3, features = 'ITGA1' )+theme(legend.position="none")

FeaturePlot(pbmc12,pt.size = 0.3, features = 'LY6K' )
VlnPlot(pbmc12, pt.size = 0.3, features ='LANCL1' )
VlnPlot(pbmc12, pt.size = 0.3, features ='CDC25C' )
#
VlnPlot(pbmc12, pt.size = 0.3, features ='COL1A1' )
#
FeaturePlot(pbmc12,pt.size = 0.3, features = 'CCNB1' )
FeaturePlot(pbmc12,pt.size = 0.3, features = 'CCNB2' )
####### plot
CairoPDF(file="11-2-VlnPlot-FeaturePlot.pdf",width=10,height=20)
#'CHCHD2','RPS11',
fts=c('KRT19','EGFR','EPCAM', 'PI3','MGP','KRT4',"S100A9", 'RPS12',"KRT17",'S100A14', # HeLa sync -
      "ZNF761", 'DEFB1', 'EXPH5', 'LYPD6B', # BC_1 +
      'BST2','PDE4D','TGFBR1','CD55','EMP3','VIM','BASP1' ,'C4BPB', #HeLa +
      'ACTA2','COL1A1','RPS4XP2', 'RPS2P28', 'COL1A2', 'PRRX1','MTUS2','PLEKHA5','COL3A1','ACTG2'); length(fts)# Hela Sync +;
#1
VlnPlot(pbmc12, pt.size = 0.3, features =fts )
#2
FeaturePlot(pbmc12,pt.size = 0.3, features = fts )#+theme(legend.position="none")
dev.off()
#

#3
top10 <- pbmc.markers %>% group_by(cluster)%>% filter(p_val_adj<0.05 & abs(avg_logFC)>log2(1.5) )  %>% top_n(n = 10, wt = avg_logFC)
dim(top10)
top10$gene

CairoPDF(file="11-3-DoHeatmap.pdf",width=5,height=6)
DoHeatmap(pbmc12, features = top10$gene,
          group.colors= c("#FF9ECE","#F81082", '#005FFF', '#98BEFD') ) 
dev.off()

#4.dot plot
DotPlotPDF=function(genelist,prefix="",width=26,height=10,pbmc=pbmc12){
  CairoPDF( paste0(prefix,"_DotPlot.pdf") , width = width, height = height)
  marker_dot_plot <- DotPlot(object = pbmc, features = genelist,
                             cols = c("white","red"), #c("black","red")
                             dot.scale = 10)
  marker_dot_plot <- marker_dot_plot + theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
  print(marker_dot_plot)
  dev.off()
}
# check
levels(pbmc12@active.ident)
table(pbmc12@active.ident)
#
DotPlotPDF( fts, "12-fts",width=14,height=5)
#
tmp0=pbmc.markers %>% group_by(cluster) %>% filter(p_val_adj<0.05 & abs(avg_logFC)>1 );dim(tmp0)
tmp=tmp0 %>% top_n(n = 20, wt = avg_logFC)
dim(tmp)
DotPlotPDF(rev(unique(tmp$gene)),"12-adjP0.05",width=20,height=5)
#
top5=tmp %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
top5
DotPlotPDF(rev(top5$gene),"13-top5",width=10,height=4.5)


# cluster specific gene list
# get sync HeLa specific gene list
tmp.HeLa_sync=pbmc.markers %>% filter(cluster=='HeLa_sync' & pct.1>0.7) %>% filter(p_val_adj<0.05 & abs(avg_logFC)>1 ); dim(tmp.HeLa_sync)
tmp.HeLa_sync=tmp.HeLa_sync[order(-tmp.HeLa_sync$avg_logFC),]
hist(tmp.HeLa_sync$pct.1, n=100)
writeLines(tmp.HeLa_sync$gene, '00_specific.HeLa_sync.gene.txt')
write.csv(tmp.HeLa_sync, '00_specific.HeLa_sync.df.csv')



############
cellInfo2=cellInfo2[order(cellInfo2$cellType),]
annotation_col=data.frame(
  row.names = rownames(cellInfo2),
  cellType = cellInfo2$cellType,
  cellType2 = cellInfo2$cellType2
  #UMAP_Cluster=rep(c(0,1,2),c(length(c0), length(c1), length(c2)))
)
annotation_col$cellType=as.character(annotation_col$cellType)
head(annotation_col)
#
ann_colors = list(
  cellType = c('BC_0'="#FF9ECE", 'BC_1'="#F81082", 'HeLa_normal'="#005FFF", 'HeLa_sync'="#619CFF"),
  cellType2 = c('BC_0'="#FF9ECE", 'BC_1'="#F81082", 'HeLa_normal'="#005FFF", 'HeLa_sync'="#619CFF",
                'HeLa_syncMix'="#D1E2FF")
)
head(ann_colors)
#
matrix2=pbmc12@assays$RNA@scale.data[tmp0$gene, ]
dim(matrix2)
#
bk <- c(seq(-3,-0.1,by=0.01),seq(0,3,by=0.01)) #[-2, 6]
CairoPDF(file="14-heatmap.pdf",width=5,height=6)
pheatmap(matrix2, border_color = NA, scale='row', #22:08 - 
         clustering_method = 'ward.D2',
         show_colnames = F, show_rownames = F,
         annotation_col = annotation_col, annotation_colors = ann_colors,
         #color = colorRampPalette( c('navy','white','firebrick3') )(100),
         color = c(colorRampPalette(colors = c("#B63FFF","black"))(length(bk)/2),
                   colorRampPalette(colors = c("black","yellow"))(length(bk)/2)),
         legend_breaks=seq(-3,3,2),
         breaks=bk,
         main="(RNA) BC and HeLa")
dev.off()
#



#########
# save rna matrix
getRNAM=function(cellInfo2){
  df=rnaM3 # 37143   225
  
  # filter by cell id
  df=df[, rownames(cellInfo2)]
  pd(df) #37143   222
  
  # filter out all 0 genes
  keep=apply(df, 1, sum)>0
  df=df[keep,]
  pd(df) #37143   222
  
  # keep genes at least in 2 cells
  keep=apply(df>0, 1, sum)>=2
  df=df[keep,]
  pd(df)
  
  ##########
  # get counts per cell
  counts.cell=as.data.frame( apply(df, 2, function(x){
    sum(x)
  }) )
  colnames(counts.cell)="counts"
  counts.cell$cid=rownames(counts.cell)
  counts.cell=counts.cell[rownames(cellInfo2),]
  #ph(counts.cell)
  cellInfo2$countsPerCell=counts.cell$counts
  
  # get gene number per cell
  geneNum.cell=as.data.frame( apply(df, 2, function(x){
    sum(x>0)
  }) )
  colnames(geneNum.cell)="geneNumber"
  geneNum.cell$cid=rownames(geneNum.cell)
  geneNum.cell=geneNum.cell[rownames(cellInfo2),]
  ph(geneNum.cell)
  cellInfo2$geneNumber=geneNum.cell$geneNumber
  
  #
  pd(cellInfo2)
  ph(cellInfo2)
  return( list(df, cellInfo2) )
}
rs=getRNAM(cellInfo2)
rnaM4=rs[[1]]
cellInfo3=rs[[2]]
dim(rnaM4) #18662   222
dim(cellInfo3) #222 7
write.csv(rnaM4, "BC_HeLa.222cells.count.V4.csv")
write.table(cellInfo3, "cellInfo.V5.txt")

# read from file
rnaM=read.csv('/data/jinwf/wangjl/apa/20200701Fig/f2/BC_HeLa/BC_HeLa.222cells.count.V4.csv', row.names = 1)
rnaM.log2cpm=apply(rnaM, 2, function(x){
  log2(x/sum(x)*1e6+1)
})


##########
# box plot for RNA
library("ggsignif")
library('ggplot2')
CairoPDF('15-boxplot_counts_geneNumber_perCell.pdf', width=3.5,height=4)
my_comparisons <- list(c("BC_0","BC_1"), c("HeLa_normal","HeLa_sync") )
g1=ggplot(cellInfo3, aes(cellType, countsPerCell/1e6, color=cellType))+
  geom_boxplot()+geom_jitter(size=0.5)+theme_bw()+
  labs(x='',y="Counts per cell (Million reads)")+
  geom_signif(comparisons = my_comparisons,
              step_increase = 0.1,
              map_signif_level = F,
              test = t.test, size=1, textsize = 5)+
  scale_color_manual('cell Type',values=c("#FF9ECE", "#F81082", "#005FFF","#619CFF"))+
  theme(axis.text.x=element_text(angle=60, hjust=1, size=10), 
        legend.position="none")+
  ylim(0,4.3)
print(g1)
#
g2=ggplot(cellInfo3, aes(cellType, geneNumber, color=cellType))+
  geom_boxplot()+geom_jitter(size=0.5)+theme_bw()+
  labs(x='',y="Gene number per cell")+
  geom_signif(comparisons = my_comparisons,
              step_increase = 0.1,
              map_signif_level = F,
              test = t.test, size=1, textsize = 5)+
  scale_color_manual('cell Type',values=c("#FF9ECE", "#F81082", "#005FFF","#619CFF"))+
  theme(axis.text.x=element_text(angle=60, hjust=1, size=10), 
        legend.position="none") +ylim(2.7e3,10.5e3)
print(g2)
dev.off()










# RNA expression level of genes ( selected by compare DPAU across groups )
##########
(function(){
  #dt=read.table('14-heatmap_all4cluster.data.txt', row.names = 1)
  dt=pbmc12@assays$RNA@data
  pd(dt)
  ph(dt[,1:5])
  #
  # add cell info
  annotation_col = data.frame(
    cellType = pbmc12@meta.data$cellType,
    id=1
  )
  rownames(annotation_col) = colnames(dt)
  annotation_col=annotation_col[order(annotation_col$cellType),]
  #
  cellInfo=cellInfo[rownames(annotation_col),]
  annotation_col2=data.frame(
    cellType2=cellInfo$cellType2,
    cellType=annotation_col$cellType
  )
  rownames(annotation_col2) = rownames(annotation_col)
  
  genelist.hm=readLines('/home/wangjl/data/apa/20200701Fig/f4/DPAU/5-gDPAU_sync_normal_hela.all.gene.txt')
  genelist.hm=unique(genelist.hm)
  #
  ttt=setdiff(genelist.hm, rownames(dt))
  print(ttt)
  genelist.hm=intersect(genelist.hm, rownames(dt))
  pl(genelist.hm)
  
  
  dt=dt[genelist.hm, rownames(annotation_col)]
  #dt=dt[1:200,]
  #dt=dt[200:nrow(dt),]
  writeLines(rownames(dt), 'tmp.txt')
  #
  # define colors
  ann_colors = list(
    cellType = c('BC_0'="#FF9ECE", 'BC_1'="#F81082", 'HeLa_normal'="#005FFF", 'HeLa_sync'="#98BEFD"), #D1E2FF
    cellType2 = c('BC_0'="#FF9ECE", 'BC_1'="#F81082", 'HeLa_normal'="#005FFF", 'HeLa_sync'="#98BEFD", 'HeLa_syncMix'='#D1E2FF') #
  )
  #head(ann_colors)
  #gap
  gaps=cumsum(as.numeric( table(annotation_col2$cellType) ) )
  bk <- seq(-5,5,0.001)
  #bk <- seq(-1.5,1.5,0.001) #c(seq(-3,-0.1,by=0.01),seq(0,3,by=0.01))
  CairoPDF(file="14-heatmap_all4cluster-data-remake-syncVSHeLa.pdf",width=10,height=50)
  #CairoPDF(file="14-heatmap_all4cluster-data-remake-syncVSHeLa-small.pdf",width=9,height=8)
  pheatmap(dt, #scale='row', 
           clustering_method = 'ward.D2',
           border_color = NA,
           show_colnames = F, 
           show_rownames = T, #gene names 
           cluster_cols = F, cluster_rows = T,
           annotation_col = annotation_col2,
           annotation_colors = ann_colors,
           gaps_col = gaps,
           #gaps_row = c(8, 21, 31),
           #color
           color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
           legend_breaks=seq(-5,5,1),
           #legend_breaks=seq(-3,3,1),
           main = "399 genes' RNA expression level across \n4 clusters (genes selected by compare \nDPAU sync_vs_normal_HeLa)",
           breaks=bk)
  dev.off()
  print('end')
})()


######important version info
#>>
# rnaM origin 37143   235
# rnaM2 rm 5 cells, counts too high, 230
# rnaM3 rm 5 cells, MT too high，225
# rnaM4 rm 3 UMAP mis cluster, remain genes expressed in at least 2 cells: 18662   222

#>>
# cellInfo #235   2  cid    cellType0
# cellInfo2 #222 5: cid cellType0 cluster cellType cellType2
# cellInfo3  # 222   7: cid cellType0 cluster cellType cellType2 countsPerCell geneNumber
