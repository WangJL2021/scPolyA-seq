# src function for f1_correlaton_gDPAU_RNA.R
# source("../src/base/tool.df.R")


##########################
# (1) get cor
##########################
get_apa_rna_cor=function(cid){
  print(paste('cell number=', length(cid)))
  #addAlpha("#FF0000FF", '22')
  gene.common=intersect(rownames(rnaM.logcpm), rownames(gDPAU))
  df.rna=rnaM.logcpm[gene.common, cid]
  df.apa=gDPAU[gene.common, cid]
  pd(df.rna) #4038  222
  pd(df.apa) #4038  222
  #
  df_cor=NULL;
  i=0
  #myColors=rainbow(1.2*nrow(df.rna))
  #plot(NULL, NULL, type='n', add=T, col='white', 
  #     xlim=c(0,15), ylim=c(0,150), xlab="RNA expression level(log2cpm)", ylab="gDPAU")
  
  start=as.numeric(Sys.time())
  
  for(gene in gene.common){
    i=i+1
    if(i%%500==0){
      print(paste(i, '; Elapse', round(as.numeric(Sys.time())-start,2), 'seconds'))
    }
    #if(i>10)break;
    rna=df.rna[gene,]
    apa=df.apa[gene,]
    # filter out NA
    rna2=as.numeric(rna[,which(!is.na(rna) & !(is.na(apa)))])
    apa2=as.numeric(apa[,which(!is.na(rna) & !(is.na(apa)))]*100)
    
    # get correlation and p
    if(length(rna2)<=2) next; # at leat 3;
    rs=cor.test(rna2, apa2, method = 'spearman')
    df.tmp=data.frame(
      gene=gene,
      n=length(rna2),
      cor=rs$estimate,
      p=rs$p.value
    )
    df_cor=rbind(df_cor, df.tmp)
    #lines(rna2, apa2, type='o', col=addAlpha(myColors[i], '33'), cex=0.5)
  }
  # adjust p and order;
  rownames(df_cor)=df_cor$gene;
  df_cor=df_cor[complete.cases(df_cor),];
  df_cor=df_cor[order(df_cor$p),];
  df_cor$adj.p=p.adjust(df_cor$p, method='fdr');
  #
  # plot
  g1=ggplot(df_cor, aes(cor, -log10(adj.p), color=n))+
    geom_point(size=0.5)+theme_bw()+
    scale_color_gradient2(low="navy", mid='lightblue', high='red', midpoint = 100)
  print(g1)
  
  return( df_cor )
} # With plot consuming 15min; No plot 56s;



##########################
# (2) save and plot cor;
##########################
save_plot_cor=function(df_cor, keyword="BC", output='cor_BC/'){
  #define significance
  print(dim(df_cor))
  # filter out too few cell items
  df_cor=df_cor[which(df_cor$n>=10),]
  #
  df_cor$sig='n.s.'
  df_cor$sig[which( (df_cor$cor>0.5) & df_cor$adj.p<0.05 )]='pos_high'
  df_cor$sig[which( (df_cor$cor< -0.5) & df_cor$adj.p<0.05 )]='neg_high'
  #df_cor$sig[which( (df_cor$cor> -0.3 & df_cor$cor< 0.3) & df_cor$adj.p<0.05 )]='low'
  tb=table(df_cor$sig);
  df_cor$sig=factor(df_cor$sig, levels=c('pos_high','n.s.', 'neg_high'));
  print(tb)
  # low      mid neg_high pos_high 
  #213     3544       69       46 
  
  #df_cor.high=df_cor[which( (df_cor$cor<0.3 & (df_cor$cor> -0.3)) & df_cor$adj.p<0.01 ),]; dim(df_cor.high)
  write.table(df_cor, paste0(output, keyword, '_01-volcano_Spearman_cor_gDPAU_RNA.df.txt') )
  
  
  #writeLines( rownames(df_cor[which(df_cor$sig=='low'),]), paste0(output, keyword, '_02-Spearman_cor_gDPAU_RNA.low.gene.txt'))
  writeLines( rownames(df_cor[which( df_cor$cor> -0.3 & df_cor$cor< 0.3 & df_cor$adj.p<0.05 ),]), paste0(output, keyword, '_02-Spearman_cor_gDPAU_RNA.low.gene.txt'))
  
  
  writeLines( rownames(df_cor[which(df_cor$sig=='pos_high'),]), paste0(output, keyword, '_02-Spearman_cor_gDPAU_RNA.pos_high.gene.txt'))
  writeLines( rownames(df_cor[which(df_cor$sig=='neg_high'),]), paste0(output, keyword, '_02-Spearman_cor_gDPAU_RNA.neg_high.gene.txt'))
  #
  g1=ggplot(df_cor, aes(cor, -log10(adj.p), color=sig))+
    geom_point(size=0.2)+theme_bw()+
    scale_color_manual('Correlation', values=c('red', 'grey80', 'blue'), limits=c('pos_high', 'n.s.', 'neg_high') )+
    #scale_color_manual('Correlation', values=c('grey80', 'grey80', 'purple', 'red'))+
    labs(x='Spearman correlation between gDPAU and RNA', title=keyword);
  library(ggrepel)
  df3=df_cor; dd_text=df3[which( ((df3$cor >0.8 | df3$cor < (-0.8)) & df3$adj.p<1e-10) | df3$adj.p<1e-18 ),];
  print('==> add text, point number: ')
  print( dim(dd_text) )
  g2=g1+geom_text_repel(data=dd_text, aes(cor, -log10(adj.p), label=rownames(dd_text)),
                        color="black",size=3,alpha=0.8)+ #xlim(-1,1)+
    guides(colour = guide_legend(override.aes = list(alpha = 1,size=2)))#放大图例的点
  #
  CairoPDF(file=paste0(output, keyword, "01-volcano_Spearman_cor_gDPAU_RNA.pdf"), width=4.5, height=4)
  print(g2)
  dev.off()
  #scale_color_gradient2(low="navy", mid='lightblue', high='red', midpoint = 120)
  # dim(df_cor.high) #90  5
  #
  print(paste('Saved: ', output, keyword))
}




##########################
# (3) by gene;
##########################
plot_apa_RNA_by_gene=function(gene, cid, x=16,y=10, xmin=0, xmax=100){
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
  #model.lm<-lm(formula = rna2 ~ apa2)
  #p0=format(summary(model.lm)$coefficients[2,4], digits = 2);p0
  
  #
  rs=cor.test(rna2, apa2, method = 'spearman')
  r=rs$estimate;#[1] -0.8427099
  r0=round(r,2);r0
  p0=format(rs$p.value, digits = 2);p0
  #p0=formatC(rs$p.value, format = "e", digits = 2)
  #plot(apa2, rna2, col='red', cex=0.5, xlim=c(0,100), # ylim=c(0,15), 
  #     xlab="apa", ylab="RNA expression level(log2cpm)")
  library(ggrepel)
  df_=data.frame( x=x, y=y, label=paste0("rho = ",r0,"\n p = ",p0) )
  g1=ggplot(data.frame(apa=apa2, rna=rna2), aes(apa, rna))+
    geom_point(size=0.3, color='black')+ #xlim(xmin, xmax)+
    labs(x='gDPAU', y='Gene expression(log2cpm)', title=gene)+theme_bw()+
    stat_smooth(method='lm')+
    #geom_text(aes(x,y,label='text'), family="Times", fontface="italic", size=4, lineheight=.8)
    geom_text(data=df_, aes(x, y, label=label), color="black",size=4, alpha=1 )
    #geom_text(data=NULL, aes(x, y, label=paste0("r = ",r0,"\n p = ",p0)), color="black",size=3,alpha=0.8 )
  print(g1)
}
#plot_apa_RNA_by_gene('NUTF2', cid.BC)



#
getRNAHighGene=function(){
  dt1=data.frame(
    gene=row.names(rnaM.logcpm),
    cv2=apply(rnaM.logcpm, 1, function(x){
      (sd(x)/mean(x) )**2
    }),
    mean=apply(rnaM.logcpm, 1, mean)
  )
  rownames( dt1[which(dt1$mean>=1),] )
}
