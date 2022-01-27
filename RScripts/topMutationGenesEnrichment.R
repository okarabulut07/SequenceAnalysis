library(clusterProfiler)
library(GSEABase)
library(org.Hs.eg.db)
filename<- "c7.all.v7.1.entrez.gmt"
gmtfile <- system.file(filename)
c6 <- read.gmt(gmtfile)
yourEntrezIdList<- c(9906,199990,54432,57705,64318,3048,26610,7439,1937,7178,23091,9445,115761,3981,83986,51005,752014,9235,7284,54849,
50488,23135,5605,9466,199786,94039,199746,23645,2357,6233,130617,754,100423062,51386,6161,9045,23317,285513,6189,1611,6206,7402,
55146,8972,6132)
ImmunSigEnrich <- enricher(yourEntrezIdList, TERM2GENE=c6, pvalueCutoff = 0.01)
ImmunSigEnrich <- setReadable(ImmunSigEnrich, OrgDb = org.Hs.eg.db, keyType =
"ENTREZID")
write.csv(ImmunSigEnrich,"MyImmunePathwayRelatedGenes.csv")
goEnrich<-enrichGO(gene= yourEntrezIdList,OrgDb= org.Hs.eg.db, ont=
"ALL",pAdjustMethod=
"BH",pvalueCutoff = 0.01,readable= TRUE)
write.csv(goEnrich,"MyGORelatedGenes.csv")
keggEnrich<-enrichKEGG(gene= yourEntrezIdList,organism= "hsa",pAdjustMethod="BH",
pvalueCutoff = 0.01)
write.csv(keggEnrich,"MyKEGGRelatedGenes.csv")