#2_5_6_motif_plot.R
# plot motif distribution

setwd("/home/wangjl/data/apa/191111Figure/f2/validation/");
#
motifs=c('AATAAA', 'ATTAAA','TATAAA','AGTAAA','AATACA','CATAAA','AATATA','GATAAA');
regions=c("PA","UTR3","extended3UTR", "exon","intron",  "intergenic","TSS", "Promoter")

#
options(stringsAsFactors=FALSE)
bedData=read.table('/data/jinwf/wangjl/apa/190705PAS/bed/pasPostions_Location_transcriptName-noChrM_noInnerPrime_PY-filterCountCell-motif.bed')
dim(bedData) #20222    12
colnames(bedData)=c('PASid','chr','pos', 'strand','count','region','transcript','down20', 'gene','up60','motif','distance')
bedData[1:5,]
bedData$distance=as.numeric(bedData$distance)
str(bedData)
#
#write.table(bedData,
# '/data/jinwf/wangjl/apa/190705PAS/bed/pasPostions_Location_transcriptName-noChrM_noInnerPrime_PY-filterCountCell-motif-Header.bed',quote = F)


#
# output bed for IGV 
igvBed=data.frame(
  chr=bedData$chr,
  start=bedData$pos,
  end=bedData$pos,
  
  pName=bedData$PASid,
  score=bedData$count,
  strand=bedData$strand
);
dim(igvBed) #20222    6
head(igvBed)
#
igvFile="/home/wangjl/data/apa/191111Figure/f2/validation/pas_bed/pasV5_20222.bed"
write.table(igvBed, igvFile, col.names = F, quote = F, row.names = F)
#




region_motif=table(bedData$region,bedData$motif);
region_motif
#
library(ggplot2)
drawStackBarPlot=function(df, posi2="stack"){
  print(dim(df))
  print(table(df$region,df$motif))
  #
  g=ggplot(df, aes(region, fill=factor(motif, levels = motifs )) )+
    geom_histogram(stat="count", position=posi2)+
    theme_classic()+
    theme(axis.text.x=element_text(angle=60, hjust=1,color='black')
    )+#,size=10,color="grey50"
    labs(x=NULL, y="Percent of pA site", title="Motif at each region")+
    scale_x_discrete(limits=regions)
  return(g)
}
p1=drawStackBarPlot(bedData)+
  scale_fill_hue("Motif")+
  labs(title="Motif at each region(including non motif)",y="Count")
p1
#
p2=drawStackBarPlot(bedData, posi2='fill')+
  scale_fill_hue("Motif")+
  ggtitle("Motif at each region(including non motif)")
p2




#
####### discard non notif rows;
bedData2=bedData[which(bedData$motif!="."),]
dim(bedData2) #[1] 13261    12

#
library(RColorBrewer)
#display.brewer.all() 
myColors=brewer.pal(9,"Set1")[1:9] 
n=8

#

#
p3=drawStackBarPlot(bedData2)+
  scale_x_discrete(limits=regions[-7])+
  #scale_fill_hue("Motif")+
  labs(title="Motif at each region",y="Count")+
  scale_fill_manual("Motif",values=myColors )
p3
#
p4=drawStackBarPlot(bedData2, posi='fill')+
  scale_x_discrete(limits=regions[-7])+
  #scale_fill_hue("Motif")+
  scale_fill_manual("Motif",values=myColors )
p4
##
#plot to file
library(gridExtra)
CairoPDF(file="2_5_6_motif_region.pdf", width=8, height=8)
grid.arrange(p1,p2,p3,p4, nrow=2)
dev.off()




###############
# test
tmpDF=bedData2[which(bedData2$motif=='AATAAA'),]
t1=density(-tmpDF$distance)
plot(t1)
#
plot(t1$x,t1$y*nrow(tmpDF), type="l", col=myColors[1], 
     xlab="Distance to cleavage site(nt)",ylab="Motif count")
# batch
##############
CairoPDF(file="2_4_6_motif_distanceToClevageSite.pdf", width=4,height=4);
# distance
for(i in 1:length(motifs) ){
  motif=motifs[i]
  print( paste(i, motif) )
  #
  tmpDF=bedData2[which(bedData2$motif==motif),]
  dim(tmpDF) # 8845   12
  #head(tmpDF)
  #
  t1=density(-tmpDF$distance)
  if(i==1){
    plot(t1$x,t1$y*nrow(tmpDF), type="l", col=myColors[i], 
         xlab="Distance to cleavage site(nt)",ylab="Motif count")
  }else{
    points(t1$x,t1$y*nrow(tmpDF), type="l", col=myColors[i])
  }
}
abline(v=-20,lty=2, col="grey")
legend(-63,0.1*8000,legend=motifs, box.col="white", col=myColors[1:8], lty=1, 
       cex=0.8,
       inset=-0.35, xpd=TRUE)
#
# count ratio
for(i in 1:length(motifs) ){
  motif=motifs[i]
  print( paste(i, motif) )
  #
  tmpDF=bedData2[which(bedData2$motif==motif),]
  dim(tmpDF) # 8845   12
  #head(tmpDF)
  #
  t1=density(-tmpDF$distance)
  if(i==1){
    plot(t1$x,t1$y, type="l", col=myColors[i], 
         xlab="Distance to cleavage site(nt)",ylab="Motif frequency")
  }else{
    points(t1$x,t1$y, type="l", col=myColors[i])
  }
}
abline(v=-20,lty=2, col="grey")
legend(-63,0.1,legend=motifs, box.col="white", col=myColors[1:8], lty=1, 
       cex=0.8,
       inset=-0.35, xpd=TRUE)
#
dev.off()
#




################# More ditails
bedData3=bedData2
head(bedData3)
bedData3$region=factor(bedData3$region, levels = regions)
bedData3$motif=factor(bedData3$motif, levels = motifs)
str(bedData3)
#
g_mr=ggplot(bedData3, aes(-distance))+
  facet_grid(motif~region)+
  geom_density()+
  labs(x="Distance to cleavage site(nt)")+
  theme_bw();g_mr

##############
#