# filter by reads and cell

setwd('/data/jinwf/wangjl/apa/191111Figure/f2/filterByReadCell/')
getwd()

library(gridExtra)

##############
# loading data

# polyA site: 
pasSite=read.table('/data/jinwf/wangjl/apa/190705PAS/bed/pasPostions_Location_transcriptName-noChrM_noInnerPrime_PY.bed',header = F)
row.names(pasSite)=pasSite$V4
pasSite[1:4,]
dim(pasSite) #458419     10

# count matrix：#742289
apaMatrix=read.table('/data/jinwf/wangjl/apa/190705PAS/bed/freq2/all225_matrix_APACountsV2.txt',
                     header = T, row.names = 1)
dim(apaMatrix) #742288    225
colnames(apaMatrix)=gsub('_','',colnames(apaMatrix))
apaMatrix[1:4,1:6]

# get the no internal priming apa matrix
apaMatrix2=apaMatrix[row.names(pasSite),]
dim(apaMatrix2) #[1] 458419    225
apaMatrix2[1:4,1:6]





############################
# QC-0: by counts per site
counts_per_site=apply(apaMatrix2,1,sum)
length(counts_per_site) #458419
head(counts_per_site)
#filter out 0 lines
table(counts_per_site==0) #none
library(Cairo)
CairoPDF(file="01_counts_per_site.pdf", width=5,height=4)
hist( log2(counts_per_site), n=100,
      xlab="log2(counts per polyA site)", ylab="polyA site number", main='Total polyA site: 285,291' )
abline(v=4, lty=2, col='red') #2**4=16
dev.off()
#

############################
#QC-1: by cell nubmer per site
expressed_cell_number=apply(apaMatrix2>0, 1 ,sum)
length(expressed_cell_number) #285291
head(expressed_cell_number)
max(expressed_cell_number)
min(expressed_cell_number)
#
CairoPDF(file="01_expressed_cell_number.pdf", width=5,height=4)
hist( log2(expressed_cell_number), n=50,
      xlab="log2(Cell number per polyA site)", ylab="Freq", 
      main='Total polyA site: 285,291', ylim=c(0, 10000) )
abline(v=log2(15), lty=2, col='red') #2**4=16




# Cell number distribution after filtering by counts/site.
table(counts_per_site>15)
# FALSE   TRUE 
# 258033 200386 
dim(apaMatrix2)
apaMatrix2[1:4,1:5]
apaMatrix3=apaMatrix2[counts_per_site>15,]
dim(apaMatrix3) #[1] 200386    225

cellNumbers=log2(apply(apaMatrix3>0, 1 ,sum))

hist( cellNumbers, n=50,
      xlab="log2(expressed cell number)", main='Total polyA site(>15 per site): 218137' )
hist( cellNumbers[which(cellNumbers>log2(225*0.1))], n=50,
      xlab="log2(expressed cell number)", main='Total polyA site(>15 per site): 218137' )
dev.off()

#
length(cellNumbers) #200386
length(cellNumbers[which(cellNumbers>log2(225*0.1))]) #20222 少了一个数量级


############################
# QC-1: filter begin

# function
library(RColorBrewer)
filterByCount_cell=function(countNum, cellNum){
  keep.f01=apply(apaMatrix2>0,1,sum)>=cellNum;
  keep.f02=apply(apaMatrix2,1,sum)>=countNum;
  
  keep.f0= (keep.f01*keep.f02 == 1)
  tb0=table(keep.f0)
  print(tb0) #true means the site would keep according to the rule, false means discarded sites.
  #
  mt=apaMatrix2[keep.f0, ]; #print( dim(mt) ) #10904   225
  
  # pie plot
  tb=table(  pasSite[row.names(mt),]$V7)
  print(tb)
  df=data.frame(
    region=names(tb),
    number=as.numeric(tb)
  )
  df$pct=round(df$number/sum(df$number)*100,2); 
  #order
  df=df[order(-df$number),]
  
  #label on Pie
  myLabel=paste0(df$number,"\n",df$pct,'%')
  myLabel[4:nrow(df)]=""
  
  # colors
  colDark2=brewer.pal(8,"Dark2")[1:8] #Dark2: 8 colors
  #barplot(rep(1,8),col= colDark2  )
  #
  n=nrow(df); #n=15
  colDark2_more= colorRampPalette(colors = colDark2)( n )
  #barplot(rep(1, n ),col= colDark2_more) #
  
  df2=df
  df2$region2=paste0(df$region, '(', df$number, ', ', df$pct, '%)')
  #df2$region2=paste0(df$region, '(', df$pct, '%)')
  df2$region2=factor(df2$region2, levels = df2$region2)
  g=ggplot(df2, mapping=aes(x="", y=number,fill=region2))+
    geom_bar( stat="identity", width = 1)+
    coord_polar("y", start=0, direction = -1)+
    theme_void()+
    #scale_fill_brewer('',palette="Dark2")+
    scale_fill_manual('', values=colDark2_more)+
    theme(
      legend.margin=margin(l = -1.5, unit='line'),
      
      legend.key.height=unit(1,"line"),
      legend.key.width=unit(0.3,"line"),
      
      plot.title=element_text(size=14, face="bold",
                              hjust = 1)
    )+
    #labs(title="Genomic location of PolyA_DB3")+ #Poly(A) sites of
    labs(title= paste0("Genomic location\n", "count>=", countNum,
                       ", cell>=",cellNum, '\nPolyA site Number=',tb0[[2]]))+
    geom_text(stat="identity",aes(x=1.1, y=number,label=myLabel),
              position=position_stack(vjust = 0.5));
  return(g)
  # return the data frame according to the rule. 
  #return (mt)
}


# final >15 reads per site(>=16); and >22.5 cells(>=23)
g=filterByCount_cell(16,23)

CairoPDF(file="B02_filter_by_count_cell-final.pdf", width=5,height=4)
g
dev.off();
## end


# save to csv
save2csv=function(countNum, cellNum){
  keep.f01=apply(apaMatrix2>0,1,sum)>=cellNum;
  keep.f02=apply(apaMatrix2,1,sum)>=countNum;
  
  keep.f0= (keep.f01*keep.f02 == 1)
  tb0=table(keep.f0)
  print(tb0) #true means the site would keep according to the rule, false means discarded sites.
  #
  mt=apaMatrix2[keep.f0, ];
  print(dim(mt))
  #mt[1:4,1:4]
  return(mt);
}
  
mt=save2csv(16,23) #[1] 20222   225
dim(mt) #[1] 20222   225
mt[1:4,1:4]
write.csv(mt, file="/data/jinwf/wangjl/apa/190705PAS/bed/apa_Matrix_afterFilterByCounts_cell_20222.csv")



######
# check apa
head(pasSite)
dim(pasSite) #[1] 458419     10
#
dim(mt) #20222   225
#
pasSite2=pasSite[rownames(mt),]
dim(pasSite2) #20222    10
pasSite2[1:10,]
#
table( is.na(pasSite2$V1) )
#
# save bed;
write.table(pasSite2,quote = F,
            file='/data/jinwf/wangjl/apa/190705PAS/bed/pasPostions_Location_transcriptName-noChrM_noInnerPrime_PY-filterCountCell.bed')
#



#



#






########## testing below
g1=filterByCount_cell(10,10)
g2=filterByCount_cell(20,10)
g3=filterByCount_cell(30,10)


CairoPDF(file="B02_filter_by_count_cell.pdf", width=12,height=4)
grid.arrange(
g1,
g2,
g3,nrow=1)
dev.off()


#more test
p2=filterByCount_cell(2,2)
p3=filterByCount_cell(3,3)
p4=filterByCount_cell(4,4)
p5=filterByCount_cell(5,5)
p6=filterByCount_cell(6,6)
p7=filterByCount_cell(7,7)

p8=filterByCount_cell(8,8)
p9=filterByCount_cell(9,9)

########
CairoPDF(file="B02_2_filter_by_count_cell-moreTest.pdf", width=12,height=12)
grid.arrange(p2, p3, p4,
             p5, p6, p7,
             p8,p9,g1,
             nrow=3)
dev.off()
