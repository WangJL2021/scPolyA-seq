gene='ZBTB38'
gene='CCDC43'
gene="MTF2"
gene="RBMS1"
gene="GRB10"
gene="SLC39A11"
gene="HIST3H2A"
gene="RPS29"
gene="IFT20"
gene="NPM1P39"
gene="UBFD1"
#
#




tmp=DPAU[gene,!is.na(DPAU[gene,])];tmp
#
tmp[,intersect(names(tmp), cid.normal)]
tmp[,intersect(names(tmp), cid.sync)]
deltaDPAU_df[gene,]


###
