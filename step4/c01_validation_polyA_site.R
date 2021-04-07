# validation polyA info V5(20222 sites)
# v0.3 
# v0.4 2020.12.9

# settings
setwd('/home/wangjl/data/apa/191111Figure/f2/validation/')
getwd()

library(Cairo)

# get data
apaInfo0=read.table("/data/jinwf/wangjl/apa/190705PAS/bed/pasPostions_Location_transcriptName-noChrM_noInnerPrime_PY-filterCountCell.bed"); 
dim(apaInfo0) #[1] 20222     10
apaInfo0[1:10, ]
#                   V1      V2      V3              V4    V5 V6   V7         V8                    V9 V10
#chr10:320403:-  chr10  320403  320403  chr10:320403:- 10945  - UTR3  DIP2C-201 TGAAAATAGCAGTTTCTTAAT   0
#chr10:855485:-  chr10  855485  855485  chr10:855485:-  3583  -   PA LARP4B-201 TCAATTTTTCTTTTGATTTTT   0
#chr10:1085971:- chr10 1085971 1085971 chr10:1085971:- 10905  - UTR3   IDI1-201 AAAGTTAAATGTTTGTTAAAT   0


#delete non-annotated site; Also can retain intergenic regions.
apaInfo0$V8=as.character(apaInfo0$V8)
#apaInfo=apaInfo0[which(apaInfo0$V8!="."),]
apaInfo=apaInfo0
dim(apaInfo) #18662    10
head(apaInfo)



#######################
# pic1: pie plot of regions
tb=table(apaInfo$V7)
tb
df=data.frame(
  region=names(tb),
  number=as.numeric(tb)
)
df$region=as.character(df$region)
df
#order
df=df[order(-df$number),]
#add percentage
df$pct= round( df$number/sum(df$number)*100, 2)
head(df)
str(df)
#add factor
df$region=factor(df$region, levels = df$region)
head(df)
#

# pie plot with annotation
library(ggrepel)
library(ggplot2)

#lable
myLabel=paste0(df$number,"\n",df$pct,'%')
myLabel[5:length(myLabel)]=""


# colors
library(RColorBrewer)
colDark2=brewer.pal(8,"Dark2")[1:8]
#barplot(rep(1,8),col= colDark2  )
#
n=nrow(df); #n=15
colDark2_more= colorRampPalette(colors = colDark2)( n )
barplot(rep(1, n ),col= colDark2_more)

#########
# draw method 1: simple
df2=df
df2$region2=paste0(df$region, '(', df$number, ')')
df2$region2=factor(df2$region2, levels = df2$region2)
#
df2$desc=paste0(df$region,'(',df$number, ', ', df$pct, '%)')
df2$desc=factor(df2$desc, levels = df2$desc)
#
head(df2)
g=ggplot(df2, mapping=aes(x="", y=number,fill=desc))+
  geom_bar( stat="identity", width = 1)+
  coord_polar("y", start=0, direction = -1)+
  theme_void()+
  #scale_fill_brewer('',palette="Dark2")+
  scale_fill_manual('', values=colDark2_more)+
  theme(
    legend.margin=margin(l = -1, unit='line'),
    
    legend.key.height=unit(1,"line"),
    legend.key.width=unit(0.5,"line"),
    
    plot.title=element_text(size=14, face="bold",
                            hjust = 0.5)
  )+
  labs(title="Genomic location of Poly(A) sites")+
  geom_text(stat="identity",aes(x=1.1, y=number, label=myLabel),
                  position=position_stack(vjust = 0.5));g
#
CairoPDF(file="01-piePlot_polyA_location.pdf",width=5.5,height=5)
g
dev.off()



#########
# draw method 2: legend blow pie

library(ggplot2)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold",
                            vjust=1,
                            hjust = 0.5)
  )
ggPie=function(df2){
  #df=df2
  g=ggplot(df2, mapping=aes(x="", y=count,fill=type))+
    geom_bar(stat="identity",width=0.5)+
    #coord_polar("y", start=0)+ 
    coord_polar(theta = 'y', start=0, direction = -1)+ 
    scale_fill_manual('', values=colDark2_more)+
    #scale_fill_brewer('',palette="Dark2")+
    blank_theme+
    theme(legend.position="bottom",
          
          legend.margin=margin(t = -2, unit='line'),
          
          legend.spacing.x = unit(2, 'pt'),
          legend.spacing.y = unit(2,"pt"),
          
          #plot.margin=unit(c(1,0,1,0),"lines"),
          
          
          legend.box = "vertical", 
          
          #legend.key.size = unit(0.2, "cm"), 
          legend.key.height=unit(0.7,"line"), 
          legend.key.width=unit(0.3,"line"), 
          
          #标签right距离，left距离
          legend.text = element_text(margin = margin(r = 10, l=2,t = -3, unit = "pt")), 
          #legend.text = element_text(margin = margin(t = -5,unit='pt')), #

          legend.background = element_blank() )+
    guides(fill = guide_legend(ncol = 3, 
                               inset=-0.1,
                               byrow=F))+ 
    labs(title="Genomic location of Poly(A) sites", 
         #subtitle ="subtitles here, pie from count data.",
         x="",y="")+
    geom_text(stat="identity",aes(x=1.1,y=count, label = tag), size=4, 
                    position=position_stack(vjust = 0.5));g
  #geom_text(stat="identity",aes(y=x, label = scales::percent(x/100)), 
  #          size=4, position=position_stack(vjust = 0.5))
  return(g)
}
CairoPDF(file="01-piePlot_2_polyA_location.pdf",width=4,height=4)
ggPie(data.frame(count=df$number, type=df$region, tag=myLabel))
dev.off()
##













#######################
# pic2: barplot, how many polyA site per gene?
head(apaInfo)

#by transcripts //Not use, use gene instead, see following part.
############################################
# annovar?  NO NO, use Python script
############################################
#tb2 = table(apaInfo$an_gene)
tb2 = table(apaInfo$V8)
head(tb2)


#
polyA_per_gene=data.frame(
  gene= as.character( names(tb2) ),
  number=as.numeric( tb2 )
)
#order
polyA_per_gene=polyA_per_gene[order(polyA_per_gene$number),]
head(polyA_per_gene)
tail(polyA_per_gene)
# remove intergenic region
polyA_per_gene=polyA_per_gene[which(polyA_per_gene$gene!='.'),]
#
tb3=table(polyA_per_gene$number)
head(tb3)
df3=data.frame(
  apaPerGene=as.character( names(tb3) ),
  number=as.numeric( tb3 )
)
df3$apaPerGene=as.character(df3$apaPerGene)
df3=df3[order(-df3$number),]
head(df3)
str(df3)
#
df4=df3[1:9,]
df4[10,]=c("M", sum(df3$number[10: nrow(df3)])  )
df4


CairoPDF(file=paste0("02-barPlot_pA-per-gene-PY_gene.pdf"),width=4,height=4.2)
posX=barplot(as.numeric( df4$number ),
             #ylim=c(0, 5000), 
             col="#DD7E00",border =NA,
             main='By transcript', #'poly(A) sites per gene',
             mgp=c(2,0.5,0),
             #xlab="poly(A) sites per gene",
             ylab="Gene number")
text(x=posX,y=-400,
     srt = 0,
     cex = 0.8,
     
     xpd = TRUE,
     labels =df4$apaPerGene)

text(x=posX, y=as.numeric( df4$number ),
     srt = 90,
     adj = c(0,0.5), 
     cex = 0.8,
     xpd = TRUE, 
     labels =df4$number)
text(x=6,y=-1300, labels=c("poly(A) sites per gene"), xpd = TRUE)

############
df5=df4[2:10,]
df5
posX=barplot(as.numeric( df5$number ),
             #ylim=c(0, 5000), 
             col="#DD7E00",border =NA,
             main='By transcript', #'poly(A) sites per gene',
             mgp=c(2,0.5,0),
             #xlab="poly(A) sites per gene",
             ylab="Gene number")
text(x=posX,y=-100, 
     srt = 0, 
     cex = 0.8,
     
     xpd = TRUE, 
     labels =df5$apaPerGene)

text(x=posX, y=as.numeric( df5$number ), 
     srt = 90, 
     adj = c(0,0.5), 
     cex = 0.8,
     xpd = TRUE, 
     labels =df5$number)
text(x=6,y=-250, labels=c("poly(A) sites per gene"), xpd = TRUE)

dev.off()
###################




############################################
# by gene
############################################
dim(apaInfo) #20222    10
head(apaInfo)

getGene=function(arr){
  rs=c();
  for(transcript in arr){
    #empty
    if(transcript=="."){
      rs=c(rs,"")
      next;
    }
    #get gene name from transcript;
    gene=strsplit(transcript, split = '-')[[1]][1]
    rs=c(rs,gene)
  }
  return(rs)
}
apaInfo$gene=getGene(as.character( apaInfo$V8 ) )
dim(apaInfo) #18662    11
head(apaInfo)
tail(apaInfo)

#####
tb2 = table(apaInfo$gene)
head(tb2)


#
polyA_per_gene=data.frame(
  gene= as.character( names(tb2) ),
  number=as.numeric( tb2 )
)
#order
polyA_per_gene=polyA_per_gene[order(polyA_per_gene$number),]
head(polyA_per_gene)
tail(polyA_per_gene)
# remove internenic region
polyA_per_gene=polyA_per_gene[which(polyA_per_gene$gene!=''),]
#
tb3=table(polyA_per_gene$number)
head(tb3)
df3=data.frame(
  apaPerGene=as.character( names(tb3) ),
  number=as.numeric( tb3 )
)
df3$apaPerGene=as.character(df3$apaPerGene)
df3=df3[order(-df3$number),]
head(df3)
str(df3)
df3$apaPerGene=as.numeric(df3$apaPerGene)
write.table(df3, file="apaPerGene_geneNumber.txt")
#
sum(df3$apaPerGene*df3$number)/sum(df3$number) #[1] 1.973979
#
df4=df3[1:9,]
df4[10,]=c("M", sum(df3$number[10: nrow(df3)])  )
df4

CairoPDF(file=paste0("03-barPlot_pA-per-gene-PY_gene-ByGene.pdf"),width=4,height=4.2)
posX=barplot(as.numeric( df4$number ),
             #ylim=c(0, 5000), 
             col="#DD7E00",border =NA,
             main='By gene', #'poly(A) sites per gene',
             mgp=c(2,0.5,0),
             #xlab="poly(A) sites per gene",
             ylab="Gene number")
text(x=posX,y=-200,
     srt = 0, 
     cex = 0.8,
     
     xpd = TRUE, 
     labels =df4$apaPerGene)

text(x=posX, y=as.numeric( df4$number ), 
     srt = 90, 
     adj = c(0,0.5), 
     cex = 0.8,
     xpd = TRUE, 
     labels =df4$number)
text(x=6,y=-700, labels=c("poly(A) sites per gene"), xpd = TRUE)

############
df5=df4[2:10,]
df5
posX=barplot(as.numeric( df5$number ),
             #ylim=c(0, 5000), 
             col="#DD7E00",border =NA,
             main='By gene', #'poly(A) sites per gene',
             mgp=c(2,0.5,0),
             #xlab="poly(A) sites per gene",
             ylab="Gene number")
text(x=posX,y=-100, 
     srt = 0, 
     cex = 0.8,
     #adj = 0.5, 
     xpd = TRUE, 
     labels =df5$apaPerGene)

text(x=posX, y=as.numeric( df5$number ), 
     srt = 90, 
     adj = c(0,0.5),
     cex = 0.8,
     xpd = TRUE,
     labels =df5$number)
text(x=6,y=-300, labels=c("poly(A) sites per gene"), xpd = TRUE)

dev.off()
