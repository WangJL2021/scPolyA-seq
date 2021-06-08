# b1_RE_
## 696 row

setwd("/home/wangjl/data/apa/20200701Fig/f2/BC_HeLa2/")
getwd()

pd=function(...){
  print(dim(...))
}
ph=function(...){
  print(head(...))
}


# load data
library(Seurat)
pbmc12=readRDS("../BC_HeLa/BC_HeLa.222cells_umap_tSNE-pbmc12.1_forUMAP_diy_addCellType.rds")

cellInfo=read.table('/data/jinwf/wangjl/apa/20200701Fig/f3/cell_cycle/cellInfo.V6.txt', row.names = 1)
head(cellInfo)

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
  # there are no sig difference between breast cancer cell line.
  annotation_col[grep("BC_", annotation_col$cellType),]$cellType="MDA-MB-468"
  #
  #
  cellInfo=cellInfo[rownames(annotation_col),]
  
  annotation_col2=data.frame(
    #cellType2=cellInfo$cellType2,
    cellType=factor(annotation_col$cellType, levels=c("MDA-MB-468", "HeLa_normal", "HeLa_sync"))
  )
  table(annotation_col2$cellType)
  rownames(annotation_col2) = rownames(annotation_col)
  
  genelist.hm=c(
    #BC
    "EGFR",  "MGP", "KRT19",'S100A9',
    "PI3", 
     "EPCAM",    "RIPK4", 'MMP7', 
    #"KRT17",  #"LINC00665", 
    
    #HeLa
    "VIM",
    "CD55", "C4BPB", 
    "EMP3", 
    "TGFBR1", 
    "BST2", 
    
    "PDE4D", 
    #"BASP1", 
    #"CTAG2", 
    #"TRABD2A", 
    #"CPS1", 
    #"PRSS21", "FSTL3", 
    
    #BC1
    #"EXPH5", "LYPD6B",
    #"ZNF551", 
    #"SLFN11", 
    #"LEAP2", 
    #"FBXO33", #"IMPAD1", 
    #'ZNF781','APOBEC3F','ZNF800',
    # BC1
    #"WDFY3",
    #"KAT6A",
    #"SECISBP2L",
    #"ANKFY1",
    #"ZNF800",
    #"PHLDB2",
    
    #hela normal
    
    
    # HeLa sync
    "ACTA2", 
    "COL1A1", #"COX5BP1" ,"COL1A2", 
    
    #"RPS4XP2", 
    #"RPS2P28", 
    
    #"PRRX1", 
    "PLEKHA5", 
    "COL3A1", 
    "ACTG2",
    
    "RPS2P28" , 
    #"COL1A1" ,  #"PRRX1",    #"COL3A1" ,  
    "EIF1P5"   #"CYCSP1"  ,"SPARC" ,   "RPS3AP12", "RPS7P10",  "COL5A2",   "ZFHX4"
  )
  genelist.hm=unique(genelist.hm)
  
  dt=dt[genelist.hm, rownames(annotation_col2)]
  writeLines(rownames(dt), 'tmp-2.txt')
  write.table(dt, '14-heatmap_all4cluster.drawHM.data-2.txt')
  #
  # define colors
  ann_colors = list(
    #cellType = c('BC_0'="#FF9ECE", 'BC_1'="#F81082", 'HeLa_normal'="#005FFF", 'HeLa_sync'="#98BEFD") #D1E2FF
    cellType = c( 'MDA-MB-468'="#F81082", 'HeLa_normal'="#005FFF", 'HeLa_sync'="#98BEFD")
    #cellType2 = c('BC_0'="#FF9ECE", 'BC_1'="#F81082", 'HeLa_normal'="#005FFF", 'HeLa_sync'="#98BEFD", 'HeLa_syncMix'='#D1E2FF') #
  )
  #head(ann_colors)
  #gap
  gaps=cumsum( as.numeric( table(annotation_col2$cellType)) )
  bk <- seq(-1.5,1.5,0.001) #c(seq(-3,-0.1,by=0.01),seq(0,3,by=0.01))
  CairoPDF(file="14-heatmap_all4cluster-data-remake-small-2.pdf",width=8, height=3.5)
  par(mar= c(5, 5, 4, 2) + 0.1 )
  pheatmap(t(dt), scale='column', 
           clustering_method = 'ward.D2',
           border_color = NA,
           show_rownames = F,
           
           cluster_cols = F, cluster_rows = F,
           annotation_row = annotation_col2,
           annotation_colors = ann_colors,
           gaps_row = gaps,
           angle_col=45,
           #gaps_row = c(8, 21, 31),
           #color
           color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
           legend_breaks=seq(-3,3,1), fontsize_row = 8,
           breaks=bk)
  dev.off()
  print('end')
})()
#



