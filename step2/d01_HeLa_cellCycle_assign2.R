# for HeLa cycle assignment.

setwd('/data/jinwf/wangjl/apa/20200701Fig/f3/cell_cycle/')

####### load data
# read matrix
rnaM=read.csv('/data/jinwf/wangjl/apa/20200701Fig/f2/BC_HeLa/BC_HeLa.222cells.count.V4.csv',
              header = T,row.names = 1)
dim(rnaM) #18662   222
rnaM[1:8,1:10]

####
# load cell info
cellInfo=read.table('/data/jinwf/wangjl/apa/20200701Fig/f2/BC_HeLa/cellInfo.V5.txt',header = T,row.names = 1)
dim(cellInfo) #222   7
head(cellInfo)

# order
table(colnames(rnaM) == row.names(cellInfo)) #222 T
#cellInfo[colnames(rna.counts),]
#
table(cellInfo$cellType)
# BC_0        BC_1 HeLa_normal   HeLa_sync 
#  92          73          30          27 





#######################
# subset
getCountsMatrix=function(keyword){
  if( !(keyword %in% c('BC','HeLa')) ){
    stop("Error: keyword must in c('BC','HeLa')")
  }
  data3=rnaM[,row.names(cellInfo[which( substr(cellInfo$cellType,1,2)==substr(keyword,1,2)  ),])]
  dim(data3) #[1] 18662    165
  #remove all 0 rows
  keep=apply(data3>0,1,sum)>0
  table(keep) #
  #FALSE  TRUE 
  #422 18240     
  data3=data3[keep,]
  print( dim(data3) ) #422 18240
  data3
}
# get HeLa
rnaM.HeLa=getCountsMatrix('HeLa') #16989    57

######
#get normal HeLa
getCountsMatrix2=function(){
  data3=rnaM[,row.names(cellInfo[which( cellInfo$cellType=="HeLa_normal" ),])]
  dim(data3) #[1] 18662    165
  #remove all 0 rows
  keep=apply(data3>0,1,sum)>0
  table(keep) #
  #FALSE  TRUE 
  #422 18240     
  data3=data3[keep,]
  print( dim(data3) ) #422 18240
  data3
}
rnaM.HeLaNormal=getCountsMatrix2() #15638    30





##################
#
outputRoot="result/HeLa2_"
phaseGenesPath="/home/wangjl/data/apa/190530Mix/cell_cycle/"

#
data2=rnaM.HeLa
data_normal=rnaM.HeLaNormal
#
#helper: return cell phase names
getCellPhaseList=function(){
  c('G1S', 'S','G2M','M','MG1');
}

# read gene set of each phase, intersect with rna matrix
geneSets=list();
for(i in 1:length(getCellPhaseList() )){
  setName=getCellPhaseList()[i]
  print( paste(i, setName) )
  tmpGenes=readLines( paste0(phaseGenesPath,setName,'.txt') )
  print( length(tmpGenes) )
  #
  geneSets[[setName]]=intersect(tmpGenes, row.names(data_normal))
  print( length(geneSets[[setName]])) #95
}


#step1: exclude genes cor<0.3 with mean of the set
geneSets2=list();
for(i in 1:length(getCellPhaseList() )){
  #gene mean of each phase
  setName=getCellPhaseList()[i] #'G1S';
  print( paste(i, setName) )
  setMean=apply(data_normal[geneSets[[setName]], ], 2,mean)
  #head(setMean)
  #cor of gene and phase mean
  tmpGenes=c()
  for(g in geneSets[[setName]]){
    rs=cor(t(data_normal[g,]), setMean)
    if(rs>=0.3){
      tmpGenes=c(tmpGenes,g)
    }
  }
  #[1] "ACD 0.324866931495592"
  #[1] "ACYP1 0.0858240480040331"
  #[1] "ADAMTS1 -0.128818193750771"
  print( length(tmpGenes) )
  geneSets2[[setName]]=tmpGenes
}
geneSets2
# write to file by phase name
for(pname in names(geneSets2) ){
  print(pname)
  writeLines( geneSets2[[pname]], paste ( outputRoot , "CycleRelatedGene_RNA_",pname,".txt", sep = "" ) )
}





##########
# use all cell: data2
dim(data2)

# step2: depth norm; log2 norm;
getNormalizedCts <- function ( cts ) {
  #ctsPath
  #cts <- read.table ( ctsPath , header = T , as.is = T )
  apply ( cts , 2 , function ( x ) { log2 ( ( 10^6 ) * x / sum ( x ) + 1 ) })
}
normCts=getNormalizedCts(data2)
dim(normCts) #16989    57
normCts[1:10,1:5]
#apply(normCts,2,sum)




#step3: calculate 5 phase scores each cell(mean of phase genes)
# Tested. Passed.
assignSampleScore <- function ( phaseGenesList , normCts ) {
  scores <- lapply ( phaseGenesList , function ( pGenes ) {
    print(length(pGenes) )
    apply ( normCts , 2 , function ( x ) {
      mean ( x [ pGenes ] )
    } )
  } )
  do.call ( cbind , scores )
}
scoreMatrix <- assignSampleScore ( geneSets2 , normCts )
head(scoreMatrix)
#               G1S        S      G2M        M      MG1
#c12ROW10 4.875341 3.848990 3.811698 4.708646 5.232650
#c12ROW16 3.776021 3.511977 3.563992 4.507177 5.757131
write.table ( scoreMatrix , paste ( outputRoot , "PhaseScores.txt" , sep = "" ) )





#step4: z-norm(each phase, then each cell)
# Tested. Passed.
getNormalizedScores <- function ( scoreMatrix ) {
  norm1 <- apply ( scoreMatrix , 2 , scale )
  normScores <- t ( apply ( t ( norm1 ) , 2 , scale ) )
  rownames ( normScores ) <- rownames ( scoreMatrix )
  colnames ( normScores )  <- colnames ( scoreMatrix )
  normScores
}
normScores <- getNormalizedScores ( scoreMatrix )
head(normScores)
#                 G1S         S         G2M          M        MG1
#c12ROW10  1.6958768 -0.8047522 -0.4522628  0.07682433 -0.51568613
#c12ROW16 -0.1256429 -1.1564363 -0.4759331  0.21754814  1.54046407
write.table ( normScores , paste ( outputRoot , "PhaseNormScores.txt", sep = "" ) )

#plot cycle score for each cell
i=4;plot(normScores[i,],type='o',col=rainbow(i), main=rownames(normScores)[i])






#step5: assign phase for each cell
getReferenceProfiles <- function () {
  referenceProfiles <- list (
    "G1S" = c ( 1 , 0 , 0 , 0 , 0 ) ,
    "G1S.S" =  c ( 1 , 1 , 0 , 0 , 0 ) ,
    "S" =  c ( 0 , 1 , 0 , 0 , 0 ) ,
    "S.G2M" =  c ( 0 , 1 , 1 , 0 , 0 ) ,
    "G2M" =  c ( 0 , 0 , 1 , 0 , 0 ) ,
    "G2M.M" =  c ( 0 , 0 , 1 , 1 , 0 ) ,
    "M" =  c ( 0 , 0 , 0 , 1 , 0 ) ,
    "M.MG1" =  c ( 0 , 0 , 0 , 1 , 1 ) ,
    "MG1" =  c ( 0 , 0 , 0 , 0 , 1 ) ,
    "MG1.G1S" =  c ( 1 , 0 , 0 , 0 , 1 ) ,
    "all" =  c ( 1 , 1 , 1 , 1 , 1 ) )
  #referenceProfiles <- lapply ( referenceProfiles , function ( x ) { names ( x ) <- c ( "G1S", "S", "G2" , "G2M" , "MG1" ); x } )
  referenceProfiles <- lapply ( referenceProfiles , function ( x ) { names ( x ) <- c ( 'G1S', 'S','G2M','M','MG1' ); x } )
  
  do.call ( rbind , referenceProfiles )
}
# Tested. Passed.
assignRefCors <- function ( normScores ) {
  referenceProfiles <- getReferenceProfiles()
  t ( apply ( normScores , 1 , function ( sampleScores ) {
    apply ( referenceProfiles , 1 , function ( refProfile ) {
      cor ( sampleScores , refProfile ) } )
  } ) )
}
# Tested. Passed.
getPhases <- function ( ) {
  #phases <- c ( "G1S", "S", "G2" , "G2M" , "MG1" )
  phases <- c ( 'G1S', 'S','G2M','M','MG1' )
  names ( phases ) <- phases
  phases
}
# Tested. Passed.
assignPhase <- function ( refCors ) {
  phases <- getPhases ()
  apply ( refCors [ ,phases ] , 1  , function ( x ) {
    phases [ which.max ( x ) ]
  } )
}




###### Score cycle similarity
refCors <- assignRefCors ( normScores )
head(refCors)
assignedPhase <- assignPhase ( refCors )
assignedPhase
#
table(assignedPhase)
#assignedPhase
#G1S G2M   M MG1   S 
#13   9  10  11  14




getDFfromNamed=function(Namedxx){
  data.frame(
    id=attr(Namedxx,'names'),
    val=unname(Namedxx),
    row.names = 1
  )
}
rs=getDFfromNamed(assignedPhase)
head(rs)
# save tag to file
write.csv(rs, paste ( outputRoot , "cellCycle_phase.csv" , sep = "" ) )

#
#compaire to setp 5?
i=5;plot(normScores[i,],type='o',col=rainbow(i), main=rownames(normScores)[i]) # get the tag of the highest phase
#





#step6 ####### Order
# Tested. Passed.
orderSamples <- function ( refCors , assignedPhase) {
  phases <- getPhases()
  orderedSamples <- list()
  
  for ( phase in phases ) {
    phaseCor <- refCors [ assignedPhase == phase , ]
    
    phaseIndex <- which ( colnames ( phaseCor ) == phase )
    if ( phaseIndex == 1 ) { preceding = ncol ( phaseCor ) - 1 } else { preceding <- phaseIndex - 1 }
    if ( phaseIndex == ncol ( phaseCor ) - 1 ) { following = 1 } else { following <- phaseIndex + 1 }
    
    earlyIndex <- phaseCor [ , preceding ] > phaseCor [ , following ]
    earlyCor <- subset ( phaseCor , earlyIndex )
    earlySamples <- rownames ( earlyCor ) [ order ( earlyCor [ , preceding ] , decreasing = T ) ]
    
    lateCor <- subset ( phaseCor , ! earlyIndex )
    lateSamples <- rownames ( lateCor ) [ order ( lateCor [ , following ] , decreasing = F ) ]
    
    orderedSamples [[ phase ]] <- c ( earlySamples , lateSamples )
  }
  refCors [ do.call ( c , orderedSamples ) , ]
}

#
ordCor <- orderSamples ( refCors , assignedPhase )
write.table ( cbind ( ordCor , "assignedPhase" = assignedPhase [ rownames ( ordCor ) ] ) , 
              paste ( outputRoot , "PhaseRefCor.txt" , sep = "" ) )
#




#step7 ####### Plot
# Passed.
plotCycle <- function  ( phaseCorsMatrix ) {
  library("pheatmap")
  library("RColorBrewer")
  breaks <- seq ( -1 , 1 , length.out = 31 )
  heatColors <- rev (brewer.pal ( 9, 'RdBu'))
  heatColors <-colorRampPalette(heatColors)
  colorPallete <- heatColors((length ( breaks ) - 1 ))
  
  # create heatmap
  hm.parameters <- list(phaseCorsMatrix,
                        border=FALSE, 
                        color = colorPallete,
                        breaks = breaks,
                        cellwidth = NA, cellheight = NA, scale = "none",
                        treeheight_row = 50,
                        kmeans_k = NA,
                        show_rownames = T, show_colnames = F,
                        main = "HeLa 57 cells",
                        clustering_method = "average",
                        cluster_rows = FALSE, cluster_cols = FALSE,
                        clustering_distance_rows = "euclidean",
                        clustering_distance_cols = NA ,
                        legend = T , annotation_legend = F )
  
  do.call("pheatmap", hm.parameters )
}


phases <- getPhases()
CairoPDF(file=paste ( outputRoot , "PhasePlot.pdf" , sep = "") , width=5 , height=3 ) #units = "in" , res = 300 
plotCycle ( t ( ordCor [ , phases ] ) )
dev.off()

jpeg ( paste ( outputRoot , "PhasePlot.jpg" , sep = "" ) , 5 , 3 , units = "in" , res = 300 )
plotCycle ( t ( ordCor [ , phases ] ) )
dev.off()





#############
# check
(function(){
  dt1=read.csv('result/HeLa2_cellCycle_phase.csv', header = T, row.names = 1)
  dt2=cellInfo
  dt2=dt2[rownames(dt1), ]
  #
  print(head(dt2))
  #
  print(table(dt1$val))
  tb=table(dt2$cellType, dt1$val);print(tb)
  #            G1S G2M  M MG1  S
  #HeLa_normal  10   3  3  11  3
  #HeLa_sync     3   6  7   0 11
  chisq.test(tb) # p-value = 0.0002153
})()
## looks like reasonable
