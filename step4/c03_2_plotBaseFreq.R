#03_2_plotBaseFreq
#v2.1
# get para from outside
#myArgs<-commandArgs(TRUE)
#cid=myArgs[1] #"all";

bedDirName="/home/wangjl/data/apa/191111Figure/f2/validation/pas_bed"
setwd(bedDirName)
getwd()

#span=200;
#keyword='200_'

drawCurve=function(bedDirName, keyword, span){
  fname=paste0(bedDirName, '/',span,'_', keyword,".fa.freq");
  
  ###########
  #load data
  ###########
  aa=read.table(fname,row.names = 1)
  #dim(aa);  aa[1:3,1:3]
  bb=as.data.frame(t(aa))
  #
  A1=bb$A
  T1=bb$U
  G1=bb$G
  C1=bb$C
  
  ###########
  #draw
  ###########
  x=-100:100
  icolor = colorRampPalette(c("#CC0000","#008000","#FFB300","#0000CC"))(4)
  #IGV color: A(009600)green; T(FF0000)red; G(D17105)orange; C(0000FF)blue
  #MEME motif colors: A(CC0000)red;  T(008000)green; G(FFB300)orange; C(0000CC)blue
  
  plot(x,A1,type='l',col=icolor[1],
       main=paste0("base freq(", sub('_','',keyword),')'),
       xlab="Relative to pA(nt)", ylab="Nucleotide frequency",lwd=0.8,
       xlim=c(-span/2,span/2), ylim=c(0,1),
       
       xaxt="n")
  axis(1,labels=seq(-100,100,20),at=seq(-100,100,20),las=2)
  
  lines(x,T1,type='l',col=icolor[2])
  lines(x,G1,type='l',col=icolor[3])
  lines(x,C1,type='l',col=icolor[4])
  #abline(v=0,lty=2, col="#eeeeee")
  #40,1
  legend('topright', box.lwd=0.25,
         lty=1, col=icolor[1:4],lwd=1,
         text.col =icolor[1:4], legend=c('A','U','G','C'))

  print(pdfpath);
}


library(Cairo)
pdfpath=paste0(bedDirName, "/../03_base_freq_200_ntUD.pdf");
CairoPDF(file=pdfpath,width=3.5*3,height=4*3)
par(mfrow=c(3,3))

# begin
drawCurve(bedDirName, keyword='all', span=200) #200

regions=c('PA', 'UTR3', 'exon', 'intergenic', 'intron', 'extended3UTR', 'Promoter', 'TSS')
for(region in regions){
  print(region)
  drawCurve(bedDirName, keyword=region, span=200) 
}
dev.off()





#debug
drawCurve2=function(bedDirName, keyword, span){
  fname=paste0(bedDirName, '/', keyword,".fa.freq");
  
  ###########
  #load data
  ###########
  aa=read.table(fname,row.names = 1)
  #dim(aa);  aa[1:3,1:3]
  bb=as.data.frame(t(aa[,1:200]))
  #
  A1=bb$A
  T1=bb$U
  G1=bb$G
  C1=bb$C
  
  ###########
  #draw
  ###########
  library(Cairo)
  pdfpath=paste0(bedDirName, "/../03.test_base_freq_",keyword,".",span,"ntUD.pdf");
  CairoPDF(file=pdfpath,width=3.5,height=4)
  x=-100:99
  icolor = colorRampPalette(c("#CC0000","#008000","#FFB300","#0000CC"))(4)
  #IGV color: A(009600)green; T(FF0000)red; G(D17105)orange; C(0000FF)blue
  #MEME motif配色： A(CC0000)red;  T(008000)green; G(FFB300)orange; C(0000CC)blue
  
  plot(x,A1,type='l',col=icolor[1],
       main=paste0("base freq(", sub('_','',keyword),' nt)'),
       xlab="Relative to pA(nt)", ylab="Nucleotide frequency",lwd=0.8,
       xlim=c(-span/2,span/2), ylim=c(0,1),
       
       xaxt="n")
  axis(1,labels=seq(-50,50,2),at=seq(-50,50,2),las=2)
  
  lines(x,T1,type='l',col=icolor[2])
  lines(x,G1,type='l',col=icolor[3])
  lines(x,C1,type='l',col=icolor[4])
  #abline(v=0,lty=2, col="#eeeeee")
  #40,1
  legend('topright', box.lwd=0.25,
         lty=1, col=icolor[1:4],lwd=1,
         text.col =icolor[1:4], legend=c('A','U','G','C'))
  dev.off()
  return(pdfpath);
}

# begin
drawCurve2(bedDirName, keyword='200_', span=20)
drawCurve2(bedDirName, keyword='200_plus', span=20) # 200 plus
drawCurve2(bedDirName, keyword='200_minus', span=20) # 200 minus
#drawCurve(bedDirName, keyword='60_pas', span=200)
