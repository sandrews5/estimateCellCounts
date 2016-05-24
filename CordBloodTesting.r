
library(minfi)
library(genefilter)
library(quadprog)

source("/dcs01/feinberglab/roadmap/CordCellSorting/CombinedPanel/estimateCellCountsRevised.r")
setwd("/dcs01/feinberglab/roadmap/CordCellSorting/CombinedPanel")
load("RGset.Canada.rda")
load("RGset.Hopkins.rda")
load("RGset.Norway.rda")
load("drop.Canada.rda")
load("drop.Hopkins.rda")
load("drop.Norway.rda")

#compositeCellType argument should be one of:
#"Hopkins"
#"Norway"
#"Canada"
#"Hopkins&Norway"
#"Hopkins&Canada"
#"Canada&Norway"
#"Hopkins&Canada&Norway"

EARLI<-get(load(file="/amber3/feinbergLab/personal/kbakulsk/EARLI/450k-round2/RGset-All-Samples-EARLI-both-rounds.rda"))
EARLI.cord<-EARLI[,EARLI$Phenotype=="blood.delivery.cord" | (EARLI$Phenotype=="Pregnancy/Sibling" & EARLI$Status=="Delivery Visit")]
EARLI.pd<-pData(EARLI.cord)
EARLI.cord<-EARLI.cord[,1:10]

test<-estimateCellCountsRevised(EARLI.cord, compositeCellType="Hopkins", cellTypes=c("CD8T","CD4T", "NK","Bcell","Mono","Gran", "nRBC"))
load("/amber3/feinbergLab/personal/kbakulsk/Cord/Projected-Cell-Types-Cord-Ref-EARLI-Fcn-20160218.rda")
sum(test-EARLI.test.counts[1:10,]) #GOOD

test<-estimateCellCountsRevised(EARLI.cord, compositeCellType="Hopkins&Norway", cellTypes=c("CD8T","CD4T", "NK","Bcell","Mono","Gran", "nRBC"))
sum(test-EARLI.test.counts[1:10,]) #GOOD

test<-estimateCellCountsRevised(EARLI.cord, compositeCellType="Hopkins&Canada&Norway", cellTypes=c("CD8T","CD4T", "NK","Bcell","Mono","Gran", "nRBC"))
sum(test-EARLI.test.counts[1:10,]) #GOOD
