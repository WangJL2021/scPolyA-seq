#

setwd("/data/jinwf/wangjl/apa/20200701Fig/f3/cell_cycle/cycle_plot/")
getwd()

library(ggplot2)


###########################
# load data after cell cycle assignment
###########################
loadData=function(fName){
  phases=c("G1S", "S", "G2M", "M", "MG1")
  getOrderedDF=function(fName){
    ranked0=read.table(fName)
    ranked=ranked0[, phases]
    ranked$phase=ranked0$assignedPhase
    ranked
  }
  
  norPhaseScore=getOrderedDF(fName)
  head(norPhaseScore)
  # check
  library(pheatmap)
  pheatmap(t(norPhaseScore[,1:5]), cluster_rows = F, cluster_cols = F)
  return(norPhaseScore)
}
# BC
fName="/data/jinwf/wangjl/apa/20200701Fig/f3/cell_cycle/result/BC_PhaseRefCor.txt"
norPhaseScore.BC=loadData(fName)
head(norPhaseScore.BC)
norPhaseScore.BC$cellType="MDA-MB-468"

# HeLa
fName="/data/jinwf/wangjl/apa/20200701Fig/f3/cell_cycle/result/HeLa_PhaseRefCor.txt"
norPhaseScore.HeLa=loadData(fName)
head(norPhaseScore.HeLa)
norPhaseScore.HeLa$cellType="HeLa"

# conbine to 1
norPhaseScore=rbind(norPhaseScore.BC, norPhaseScore.HeLa)
head(norPhaseScore)
dim(norPhaseScore) #222 7

# add cell info
cellInfo=read.table("/data/jinwf/wangjl/apa/20200701Fig/f3/cell_cycle/cellInfo.V6.txt")
head(cellInfo)
cellInfo=cellInfo[rownames(norPhaseScore),]
norPhaseScore$cellType2=cellInfo$cellType
#
norPhaseScore$phase=factor(norPhaseScore$phase, levels=c('G1S','S', 'G2M', 'M', 'MG1'))
#
main=paste0("Cell cycle score", "(", nrow(norPhaseScore) ," cells)")
g2=ggplot(norPhaseScore, aes( G1S+S, G2M+M, color=phase, shape=cellType2 ))+
  geom_point(size=1.8) + 
  theme_classic()+
  scale_shape_discrete("Cell cluster")+
  scale_color_manual("Phase", values=c(
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3","#A6D854"
  ))+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=2)));g2 # best

pdf(paste0("03_cycleScore_cyclePlot_AllCells-large.pdf"), width=3.8, height=3.4)
print(g2)
dev.off()
