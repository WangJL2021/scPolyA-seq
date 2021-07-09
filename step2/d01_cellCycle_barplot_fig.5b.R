# 02_cellCycle_visual.R
# version 2: add ns and ***

setwd('/data/jinwf/wangjl/apa/20200701Fig/f3/cell_cycle/')
getwd()

###############
## set cell cycle tag to cellInfo, save to file
cellInfo=(function(){
  df1=read.csv('result/BC_cellCycle_phase.csv', header = T, row.names = 1)
  df2=read.csv('result/HeLa2_cellCycle_phase.csv', header = T, row.names = 1)
  
  df=rbind(df1, df2)
  df$cid=rownames(df)
  print( dim(df) )
  # head(df)
  
  #
  cellInfo=read.table('/data/jinwf/wangjl/apa/20200701Fig/f2/BC_HeLa/cellInfo.V5.txt',header = T,row.names = 1)
  df=df[rownames(cellInfo),]
  cellInfo$cellCycle=df$val
  ## 
  write.table(cellInfo, 'cellInfo.V6.txt')
  cellInfo
})()
dim(cellInfo)
head(cellInfo)



#########
table(cellInfo$cellCycle)
#
cellInfo[which(substr(cellInfo$cellType,1,2)=="BC"),]$cellType="MDA-MB-468"
tbl1=table(cellInfo$cellType, cellInfo$cellCycle)
tbl1
#            G1S G2M  M MG1  S
#HeLa_normal  10   3  3  11  3
#HeLa_sync     3   6  7   0 11
#MDA-MB-468   47  29 25  39 25

tbl1=tbl1[ ,c('G1S','S', 'G2M', 'M', 'MG1')]
tbl1=tbl1[c(3,1,2),c(5,4,3,2,1)] #reorder columns
tbl1

#####
# test
p1=chisq.test(tbl1[1:2,]); p1 #p-value = 0.4465
p2=chisq.test(tbl1[2:3,]); p2 #p-value = 0.0002153=2.15e-4

p2=formatC(p2$p.value, format = "e", digits = 2)
p2


#percentage
tbl2=apply(tbl1,1,function(x){ 100*x/sum(x)})
cellNames=colnames(tbl2)
colnames(tbl2)=NULL
tbl2

# barplot
library(RColorBrewer)
col=RColorBrewer::brewer.pal(n = 5,name = "Set2");
col=rev(col)
barplot(1:5,col=col)
#plot
library(Cairo)


CairoPDF(file='01_cellCycle_barplot.pct-withP.pdf', width=2.8, height=4)
par(mar=c(5, 4, 5, 5) + 0.1)
posX=barplot(as.matrix(tbl2), col=col, mgp=c(2.5,1,0),
             names.arg=NULL,
             space=0.2, # bar space
             ylab="Percentage",
             border=NA)
# x,y range
usr=par('usr'); usr
# legend
legend(usr[2], usr[4]*0.9, border=NA,
       fill=rev(col), legend=rev(rownames(tbl2)),
       box.col="white", inset=-0.35,bty="n", xpd=T)

#add x axis label
text(posX, usr[3]*5.5, labels=cellNames, adj=1, srt=45, xpd=TRUE)
# p value
segments(posX[1], usr[4]*1.03, posX[2], xpd=T, lwd=2)
segments(posX[2], usr[4]*1.1, posX[3], xpd=T, lwd=2)
#
text( mean(posX[1:2]), usr[4]*1.08, labels="n.s.", xpd=T, cex=1.3)
text( mean(posX[2:3]), usr[4]*1.16, labels=paste0("***"), xpd=T, cex=1.3)
#box()
dev.off()








#############
#1 cell_number_at_each_phase
outPath=''
library(RColorBrewer)
col6=RColorBrewer::brewer.pal(n = 6,name = "Set2")

tb2=table(cellInfo$cellType, cellInfo$cellCycle)
tb2

library(ggplot2)
theme_blank=theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

CairoPDF( paste0('02_cell_number_at_each_phase_coloredSe2.pdf'),height=5,width=5)
par(mfrow=c(2,2), mai=c(0,0,0.5,0))
for(i in 1:4){
  tb2.tmp=tb2[i,]; #tb2.tmp
  tb2.tmp2=data.frame(
    "phase"=names(tb2.tmp),
    "cellNumber"=as.numeric(tb2.tmp)
  )
  #order
  tb2.tmp2$phase=as.character(tb2.tmp2$phase)
  rownames(tb2.tmp2)=tb2.tmp2$phase
  tb2.tmp2=tb2.tmp2[c('G1S',"S","G2M","M","MG1"),]
  tb2.tmp2$phase=factor(tb2.tmp2$phase,levels=c('G1S',"S","G2M","M","MG1"))
  #label text
  tb2.label=paste0(tb2.tmp2$phase,"(",
                   round(100*tb2.tmp2$cellNumber/sum(tb2.tmp2$cellNumber),0),"%",
                   ")" 
  );
  tb2.tmp2$label=tb2.label
  head(tb2.tmp2)
  pie(tb2.tmp2$cellNumber, labels=tb2.label,
      col=col6,
      clockwise=T,init.angle=90,
      radius =0.6,
      main=paste0('', rownames(tb2)[i],'') )
  #
  #must be factor! for color or fill!
  library(ggrepel)
  g.tmp=ggplot(tb2.tmp2, mapping=aes(x="", y=cellNumber,fill=phase))+
    geom_bar(stat="identity",width=0.5)+
    coord_polar(theta = 'y', direction = -1)+
    scale_fill_manual(values=col6)+
    theme_blank+
    theme(
      axis.title.x = element_blank(), #no x caption
      axis.title.y = element_blank(),
      #
      axis.text.x = element_blank(), #no x axis text
      axis.text.y = element_blank(),
      axis.ticks = element_blank(), 
      #
      axis.line=element_blank()
    )+
    labs(title=paste0(rownames(tb2)[i],'') )+
    geom_text_repel(stat="identity",aes(x=1.3,y=cellNumber, label = label), size=4,
                    position=position_stack(vjust = 0.5))
  assign(paste0('gtmp',i),g.tmp )
}
dev.off()


#
# ggplot2: wjl=c(gtmp1,gtmp2,gtmp3,gtmp4)
CairoPDF( paste0(outPath,'02_cell_number_at_each_phase_coloredSe2_ggplot2.pdf'),height=6,width=6)
library(grid)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,2)) )
#
vplayout=function(x,y){
  viewport(layout.pos.row=x, layout.pos.col=y)
}
print(gtmp1,vp=vplayout(1,1))
print(gtmp2,vp=vplayout(1,2))
print(gtmp3,vp=vplayout(2,1))
print(gtmp4,vp=vplayout(2,2))
dev.off()
#

################
# end
