
library(minfi)
setwd("/dcs01/feinberglab/roadmap/CordCellSorting/CombinedPanel")
load("/dcs01/feinberglab/roadmap/CordCellSorting/PackageCreation/deGoede_nRBC_paper/deGoede_RGset.rda") #RGset.C
load("/dcs01/feinberglab/roadmap/CordCellSorting/PackageCreation/deGoede_nRBC_paper/Hopkins_RGset.rda") #RGset.H
load("/dcs01/feinberglab/roadmap/CordCellSorting/PackageCreation/deGoede_nRBC_paper/Norway_RGset.rda") #Rgset.N

# > dim(RGset.C)
# Features  Samples
  # 622399       62
# > dim(RGset.H)
# Features  Samples
  # 622399      104
# > dim(RGset.N)
# Features  Samples
  # 622399       77

canada<-pData(RGset.C)
hopkins<-pData(RGset.H)
norway<-pData(RGset.N)
# > colnames(canada)
# [1] "Sample"   "Rownames" "Cell"     "Method"   "Indiv"    "Dataset"
# > colnames(hopkins)
 # [1] "X"             "Plate_ID"      "Sample.Well"   "Hyb.date"
 # [5] "Image.date"    "Sex"           "Age"           "Individual.ID"
 # [9] "Array"         "Slide"         "Basename"      "predictedSex"
# [13] "CellType"      "Dataset"
# > colnames(norway)
# [1] "Sample_Name"  "SampleID"     "CellTypeLong" "CellType"     "Sex"
# [6] "Array"        "Slide"        "Basename"     "Dataset"

table(canada$Cell)
table(hopkins$CellType)
table(norway$CellType)
# > table(canada$Cell)

      # Bcell    CD4Tcell    CD8Tcell Granulocyte    Monocyte      NKcell
          # 7           7           6           7          12           6
       # nRBC       Tcell
         # 12           5
# > table(hopkins$CellType)

     # Bcell       CD4T       CD8T       Gran       Mono         NK       nRBC
        # 15         15         14         12         15         14          4
# WholeBlood
        # 15
# > table(norway$CellType)

# Bcell  CD4T  CD8T  Gran  Mono    NK   WBC
   # 11    11    11    11    11    11    11

canada$CellType<-as.character(canada$Cell)
canada$CellType<-ifelse(canada$CellType=="CD4Tcell","CD4T",canada$CellType)
canada$CellType<-ifelse(canada$CellType=="CD8Tcell","CD8T",canada$CellType)
canada$CellType<-ifelse(canada$CellType=="Granulocyte","Gran",canada$CellType)
canada$CellType<-ifelse(canada$CellType=="Monocyte","Mono",canada$CellType)
canada$CellType<-ifelse(canada$CellType=="NKcell","NK",canada$CellType)
canada$CellType<-factor(canada$CellType)
norway$CellType<-ifelse(norway$CellType=="WBC","WholeBlood",norway$CellType)
norway$CellType<-factor(norway$CellType)

table(canada$CellType)
table(hopkins$CellType)
table(norway$CellType)
class(canada$CellType)
class(hopkins$CellType)
class(norway$CellType)
levels(canada$CellType)
levels(hopkins$CellType)
levels(norway$CellType)

#Make sure the levels are the same for each cell type. First the content:
levels(canada$CellType)<-c(levels(canada$CellType),"WholeBlood")
levels(hopkins$CellType)<-c(levels(hopkins$CellType),"Tcell")
levels(norway$CellType)<-c(levels(norway$CellType),"nRBC","Tcell")
#Now the order:
canada$CellType<-factor(canada$CellType,levels=levels(hopkins$CellType))
norway$CellType<-factor(norway$CellType,levels=levels(hopkins$CellType))

#Reassign the pData
pData(RGset.C)<-canada
pData(RGset.H)<-hopkins
pData(RGset.N)<-norway

# #make sure combining works
trial<-combine(RGset.C,RGset.H) #GOOD
pData(trial)$CellType
trial<-combine(RGset.H,RGset.N) #FIX IT
Error in combine(pDataX, pDataY) : data.frames contain conflicting data:
        non-conforming colname(s): Sex, Array, Slide, Basename
norway$Sex<-factor(norway$Sex,levels=c("F","M"))
hopkins$Array<-as.character(hopkins$Array)
norway$Slide<-as.numeric(norway$Slide)
hopkins$Basename<-as.character(hopkins$Basename)
		
trial<-combine(RGset.C,RGset.N)	#GOOD
	
#Re-assign all the pdata
pData(RGset.C)<-canada
pData(RGset.H)<-hopkins
pData(RGset.N)<-norway
trial<-combine(RGset.H,RGset.N) #GOOD!
	
	
#make the naming easier for the function
RGset.Canada<-RGset.C
RGset.Hopkins<-RGset.H
RGset.Norway<-RGset.N

save(RGset.Canada,file="RGset.Canada.rda")
save(RGset.Hopkins,file="RGset.Hopkins.rda")
save(RGset.Norway,file="RGset.Norway.rda")

#Not able to get back to Kelly's # of probes removed. She started here
load("/dcs01/feinberglab/roadmap/CordCellSorting/RGset-Cord-All.rda")
	#let's see if we can find the differences
	#okay well could just be the one kelly used had more samples (controls, etc.)
	#Let's just go with the ones we determine in here

#We have to clean all of these and produced "bad probes" objects for each
#Cross reactive:
load("/legacy/amber2/scratch/claddaco/SEEDmethylation/Gaphunting/cross.probes.info.rda")
crossrx<-as.character(cross.probes.info$TargetID)
#Sex chromosomes
temp<-preprocessRaw(RGset.C)
temp<-mapToGenome(temp)
chrnames<-as.character(seqnames(temp))
sexchr<-rownames(temp)[which(chrnames=="chrX"|chrnames=="chrY")]
#SNP mapping probes
anno<-getAnnotation(temp)
cpg<-rownames(anno)[which(!is.na(anno$CpG_rs))]
sbe<-rownames(anno)[which(!is.na(anno$SBE_rs))]
dropprobes<-unique(c(crossrx,sexchr,cpg,sbe))

##################################
#Detection p-value
canada.detp<-detectionP(RGset.C)
	failed.canada<-canada.detp>0.01
	samplesfailed.canada<-which(colMeans(failed.canada)>0.01) #0 samples failed via detection p-value
	probesfailed.canada<-which(rowMeans(failed.canada)>0.1) #718 probes failed via detection p-value
	drop.Canada<-unique(c(dropprobes,names(probesfailed.canada))) #55715

hopkins.detp<-detectionP(RGset.H)
	failed.hopkins<-hopkins.detp>0.01
	samplesfailed.hopkins<-which(colMeans(failed.hopkins)>0.01) #0 samples failed via detection p-value
	probesfailed.hopkins<-which(rowMeans(failed.hopkins)>0.1) #674 probes failed via detection p-value
	drop.Hopkins<-unique(c(dropprobes,names(probesfailed.hopkins))) #55753
	
norway.detp<-detectionP(RGset.N)
	failed.norway<-norway.detp>0.01
	samplesfailed.norway<-which(colMeans(failed.norway)>0.01) #0 samples failed via detection p-value
	probesfailed.norway<-which(rowMeans(failed.norway)>0.1) #876 probes failed via detection p-value
	drop.Norway<-unique(c(dropprobes,names(probesfailed.norway))) #55776

save(drop.Canada,file="drop.Canada.rda")
save(drop.Hopkins,file="drop.Hopkins.rda")
save(drop.Norway,file="drop.Norway.rda")
	
#################################################################################################
#################################################################################################
#Based on all this below it doesn't seem like somethings gets screwed up because there are other
#cell types besides the ones we are interested in. This makes sense of course...
mSet<-preprocessQuantile(trial)
p <- getBeta(mSet)
pd <- as.data.frame(pData(mSet))
cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran")
pd$CellType <- factor(pd$CellType, levels = cellTypes)
library(matrixStats)
library(genefilter)
ffComp <- rowFtests(p, pd$CellType)
    prof <- sapply(splitit(pd$CellType), function(i) rowMeans(p[,i]))
