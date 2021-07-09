## load data after b01


## annotate p value in figure
## v2.0: if p>=0.01, then show p directly with 3 digits.
# if p<0.01, then using scientific mark, with 2 digits.
transPvalue2=function(p0){
  #p0=1.2e-52; p0
  if(p0<0.01){
    p0=formatC(p0, format="e", digits=2);#p0 #"1.2e-52"
    # split by e
    p0=strsplit(p0, "e")[[1]]; #p0 #[1] "1.2" "-52"
    # plot
    label = paste0( 'italic(P)~"="~',p0[1],'~"\u00d7"~10^',p0[2])
  }else{
    p0=round(p0,3); #p0 #保留2位
    label = paste0( 'italic(P)~"=" ~', p0 )
  }
  return(label)
}



###
plot_apa_RNA_by_gene2=function(gene, cid, x=16,y=10, xmin=0, xmax=100){
  #addAlpha("#FF0000FF", '22')
  gene.common=intersect(rownames(rnaM.logcpm), rownames(gDPAU))
  if(!(gene %in% gene.common)){
    stop( paste(gene, " not in both matrix.") )
  }
  rna=rnaM.logcpm[gene, cid]
  apa=gDPAU[gene, cid]
  # filter out NA
  rna2=as.numeric(rna[,which(!is.na(rna) & !(is.na(apa)))])
  apa2=as.numeric(apa[,which(!is.na(rna) & !(is.na(apa)))]*100)
  # get correlation and p
  if(length(rna2)<3){stop(paste('Error: less than 3 observation, can NOT get correlation.'))}; # at leat 3;
  
  #
  rs=cor.test(rna2, apa2, method = 'spearman')
  r=rs$estimate;#[1] -0.8427099
  r0=round(r,2);r0
  p0=rs$p.value;
  

  ggplot(data.frame(apa=apa2, rna=rna2), aes(apa, rna))+
    geom_point(size=0.3, color='black')+ #xlim(xmin, xmax)+
    labs(x='gDPAU', y='Gene expression(log2cpm)', title=gene)+theme_bw()+
    stat_smooth(method='lm')+
    annotate(geom="text", x, y, label=paste0("rho = ",r0),  )+
    annotate(geom="text", x, y-0.7, label=transPvalue2(p0), parse=T )
}


pdf("cor_BC/BC03_point_apa_RNA-3.pdf", width=3, height=3)
plot_apa_RNA_by_gene2('NUTF2', cid.BC, x=25)
plot_apa_RNA_by_gene2('CDKN2B', cid.BC, y=6, x=25)
# neg
plot_apa_RNA_by_gene2('BLOC1S5', cid.BC, y=1, x=25)
plot_apa_RNA_by_gene2('KRAS', cid.BC, x=26, y=2.5)
# plot_apa_RNA('BST2')
plot_apa_RNA_by_gene2('CDK6', cid.BC,x=26, y=3 )
# plot_apa_RNA('CCNB2',y=9)
plot_apa_RNA_by_gene2('MRPL43', cid.BC, y=1, x=25)
plot_apa_RNA_by_gene2('ISCU', cid.BC, y=1, x=25)
dev.off()
# 
