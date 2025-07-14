
####################################################################################
####################################################################################
### (1) CIRCOS-Plots including CN_Gains, CN_Losses and UPDs:
####################################################################################
####################################################################################

#try(library(knitr))
#try(opts_chunk$set(fig.pos = ""))
library(circlize)


### Fix COLUMN-NAMES:
#colnames(segments) <- c("chr",)



###### ###### ###### ###### ###### ######
### INPUT-dataframes -- Including PERIS, AECC & Stage-II (in total, N=262, including "Close_mucosa" samples):
###### ###### ###### ###### ###### ######

###### ####### ######
### INPUT dataframes: 
###### ####### ######
setwd("~/Desktop/Analysis/Prediction_Tumor-stage_Integrative_2022/Segments-CNAs_Cohorts")

###### ###### ######
###### 1) Cohort PERIS (N=92):
peris <- read.table("PERIS_1_Segments_FACETS_All-samples_5May2022.txt",
                    header=TRUE,stringsAsFactors=FALSE)
colnames(peris)[c(1,5)] <- c("chrom","cnlr.median")
dim(peris)
#[1] 9457   16


###### ###### ######
###### 2) Cohort AECC (N=86):
aecc <- read.table("AECC_1_Segments_FACETS_All-samples_17Mar2022.txt",
                   header=TRUE,stringsAsFactors=FALSE)
colnames(aecc)[c(1,5,9,10,13)] <- c("chrom","cnlr.median","tcn.em","lcn.em","Lesion")
dim(aecc)
#[1] 6421   13


###### ###### ######
###### 3) Cohort Stage II (N=84):
#st2 <- read.table("StageII_Nexus_segments_ReSegmented_Clínic-cohort_84samples_17June2020.txt",
#                    header=TRUE,stringsAsFactors=FALSE)
#st2$Lesion <- "cancer"
#colnames(st2)[c(2:4,7)] <- c("chrom","start","end","cnlr.median")
#dim(st2)
#[1] 7827    7


### ### ### ### ### ###
### Nexus-segments of the Stage-II cohort [N=84]:

### PATIENTS included in the initial Nexus-analysis:
setwd("~/Desktop/Analysis/FIS_Stage-II_2019/SNP-arrays/CNA-Load_Stage-II/Recalculate_CNA-Load_June2020/Outputs/15June2020/Nexus_segments")
patients <- read.table("NAMES_Nexus-segments-folders_June2020.txt",
                       header=TRUE,sep="\t",stringsAsFactors=FALSE)
dim(patients)
#[1] 92  1

### LOOP:
all_files <- data.frame()

#
for(i in 1:nrow(patients))
{
  folder <- patients[i,"sample"]
  
  ### Log2Ratios:
  logr_route <- paste("~/Desktop/Analysis/FIS_Stage-II_2019/SNP-arrays/ASCAT_CCF_StageII/Desktop_files/Nexus_segments_all/",folder,"/segments.txt",
                      sep="")
  logr_file <- read.table(logr_route,header=TRUE,sep="\t",stringsAsFactors=FALSE)
  colnames(logr_file)[4] <- "mean_LogR"
  logr_file$sample <- folder
  
  ### BAFs:
  baf_route <- paste("~/Desktop/Analysis/FIS_Stage-II_2019/SNP-arrays/ASCAT_CCF_StageII/Desktop_files/Nexus_segments_all/",folder,"/snpsegments.txt",
                     sep="")
  baf_file <- read.table(baf_route,header=TRUE,sep="\t",stringsAsFactors=FALSE)
  
  ### JOIN the BAF-column to the previous LogR-file:
  logr_file$mean_BAF <- baf_file[,4]
  
  ### JOIN all segment-files:
  all_files <- rbind(all_files,logr_file)
  
  print(paste(i,folder,sep=" "))
}

dim(all_files)
#[1] 71572     6
levels(factor(all_files$sample)) %>% length
#[1] 92


### MASTER-Table of the Stage-II cohort (N=84):
setwd("~/Desktop/Analysis/FIS_Stage-II_2019/Master_Table_All-Clinico-Genomic-markers")
master <- read.table("Master-Table_Stage-II_84-patients_Clinico-Genomic-Mutational_28Oct-2020_GOOD.txt",
                     header=TRUE,sep="\t",stringsAsFactors=FALSE)
dim(master)
#[1] 84 40


### PATIENTS to be included in the final Stage-II analysis (N=84):
patients$sample <- gsub("_","",patients$sample)
patients$sample[patients$sample=="T09"] <- "T9"

#
all_files$sample <- gsub("_","",all_files$sample)
all_files$sample[all_files$sample=="T09"] <- "T9"

#
all_files2 <- all_files[all_files$sample %in% master$sample,]
levels(factor(all_files2$sample)) %>% length
#[1] 84

#
#setwd("~/Desktop/Analysis/FIS_Stage-II_2019/SNP-arrays/ASCAT_CCF_StageII/Desktop_files/Nexus_segments_all")
setwd("~/Desktop/Analysis/Prediction_Tumor-stage_Integrative_2022/Segments-CNAs_Cohorts")

#write.table(all_files2,"StageII_Nexus_segments_Clínic-cohort_84samples_25June2022.txt",
#            col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE)


### Parameter "length", and ELIMINATION of segments w/ Length <2000000:
all_files2$length <- all_files2$End - all_files2$Start +1

#
st2 <- all_files2[all_files2$length >= 2000000,]
dim(st2)  # Filtered-out segments w/ length <2000000
#[1] 22431     7


### COLNAMES of "st2":
colnames(st2) <- c("chrom","start","end","cnlr.median","sample","mean_BAF","length")



###### ####### ######
### RE-CENTRALIZE Log2Ratio values to the MEDIAN of each tumor/lesion:
###### ####### ######
#
#segments <- peris
segments <- aecc
#segments <- st2

#
for(i in 1:nrow(segments)){
  
  sample <- segments$sample[i]
  median <- median(segments$cnlr.median[segments$sample==sample])
  
  segments$log2_rec[i] <- segments$cnlr.median[i] - median
}

# All log2_rec MEDIANS should be near 0:
median(segments$log2_rec[segments$sample=="P2"])
#[1] 0

# 
#peris <- segments
aecc <- segments
#st2 <- segments



###### ###### ###### ###### ###### ######
### Create INPUT-Dataframes for circlize (equal to DNAcopy-derived segments, 
### such as "AECC_AE_CNL_agi_31_05_2019.txt"):
###### ###### ###### ###### ###### ######

###### ###### ######
### 0- How many CNAs of each-type?
###### ###### ######

### CN_Gains:
nrow(segments[segments$log2_rec >= 0.18,])
#[1] 944
head(segments[segments$log2_rec >= 0.18,c(1,5:6,13:14,18)])
sum(segments[segments$log2_rec >= 0.18,"length"],na.rm=TRUE) / sum(segments$length)
#[1] 0.09143614


### CN_Losses:
nrow(segments[segments$log2_rec <= -0.18,])
#[1] 842
head(segments[segments$log2_rec <= -0.18,c(1,5:6,13:14,18)])
sum(segments[segments$log2_rec <= -0.18,"length"],na.rm=TRUE) / sum(segments$length)
#[1] 0.08254297


### UPDs/UPPs:
nrow(segments[segments$log2_rec <= 0.18 & segments$log2_rec >= -0.18 & segments$tcn.em==2 & segments$lcn.em==0,])
#[1] 1684
head(segments[segments$log2_rec <= 0.18 & segments$log2_rec >= -0.18 & segments$tcn.em==2 & segments$lcn.em==0,c(1,5:6,13:14,18)])
sum(segments[segments$log2_rec <= 0.18 & segments$log2_rec >= -0.18 & segments$tcn.em==2 & segments$lcn.em==0,"length"],na.rm=TRUE) / sum(segments$length)
#[1] 0.03731358

# Range of BAFs for regions with UPD/UPP:
range(segments[segments$log2_rec <= 0.18 & segments$log2_rec >= -0.18 & segments$tcn.em==2 & segments$lcn.em==0,"mafR"],na.rm=TRUE)
#[1] -0.01149579  3.09888131
mean(segments[segments$log2_rec <= 0.18 & segments$log2_rec >= -0.18 & segments$tcn.em==2 & segments$lcn.em==0,"mafR"],na.rm=TRUE)
#[1] 0.580911


### [Else] -- Diploid regions:
nrow(segments[segments$log2_rec <= 0.18 & segments$log2_rec >= -0.18 & segments$lcn.em!=0,])
#[1] 4407
head(segments[segments$log2_rec <= 0.18 & segments$log2_rec >= -0.18 & segments$lcn.em!=0,c(1,5:6,13:14,18)])
sum(segments[segments$log2_rec <= 0.18 & segments$log2_rec >= -0.18 & segments$lcn.em!=0,"length"],na.rm=TRUE) / sum(segments$length)
#[1] 0.6358174

# Range of BAFs for Diploid-regions:
mean(segments[segments$log2_rec <= 0.18 & segments$log2_rec >= -0.18 & segments$lcn.em!=0,"mafR"],na.rm=TRUE)
#[1] 0.05508227



###### ###### ###### ######
### 1- Define in "segments" each segment-CNA (using re-centralized Log2R):
###### ###### ###### ######

### 1.a) For FACETS-derived segments:
#segments <- peris
segments <- aecc

#
for(i in 1:nrow(segments)){
  
  ### CRITERIA_1 (when the lcn.em=minor_CN is NA):
  if( is.na(segments$lcn.em[i]) | is.na(segments$tcn.em[i]) )
  {
    ## CN_Gain regions:
    if( segments$log2_rec[i] >= 0.18 )
    { segments$CNA_status[i] <- "CN_Gain" }
    
    ## CN_Loss regions:
    else if( segments$log2_rec[i] <= -0.18 )
    { segments$CNA_status[i] <- "CN_Loss" }
    
    ## UPD/UPP regions:
    else if(segments$log2_rec[i] <= 0.18 & segments$log2_rec[i] >= -0.18 & segments$mafR[i] >=1 )
    { segments$CNA_status[i] <- "UPD" }
    
    ## Diploid regions:
    else { segments$CNA_status[i] <- "Diploid" }
  }
  
  
  ### CRITERIA_2 (when the lcn.em=minor_CN is NOT NA):
  else {
    ## CN_Gain regions:
    if( segments$log2_rec[i] >= 0.18 )
    { segments$CNA_status[i] <- "CN_Gain" }
    
    ## CN_Loss regions:
    else if( segments$log2_rec[i] <= -0.18 )
    { segments$CNA_status[i] <- "CN_Loss" }
    
    ## UPD/UPP regions:
    else if(segments$log2_rec[i] <= 0.18 & segments$log2_rec[i] >= -0.18 & segments$tcn.em[i]==2 & segments$lcn.em[i]==0 )
    { segments$CNA_status[i] <- "UPD" }
    
    ## Diploid regions:
    else { segments$CNA_status[i] <- "Diploid" }
  }
}

#
table(segments$CNA_status) ## PERIS
#CN_Gain CN_Loss Diploid     UPD 
#   1846    1558    5710     343 
table(segments$CNA_status) ## AECC 
#CN_Gain CN_Loss Diploid     UPD 
#   1596    1213    3531      81 

#
#peris <- segments
aecc <- segments


### 1.b) For NEXUS-derived segments (Stage II cohort):
segments <- st2

#
for(i in 1:nrow(segments)){
  
  ## CN_Gain regions:
  if( segments$log2_rec[i] >= 0.18 )
  { segments$CNA_status[i] <- "CN_Gain" }
  
  ## CN_Loss regions:
  else if( segments$log2_rec[i] <= -0.18 )
  { segments$CNA_status[i] <- "CN_Loss" }
  
  ## UPD/UPP regions:
  else if(segments$log2_rec[i] <= 0.18 & segments$log2_rec[i] >= -0.18 & segments$mean_BAF[i] <=0.3 )
  { segments$CNA_status[i] <- "UPD" }
  
  ## Diploid regions:
  else { segments$CNA_status[i] <- "Diploid" }
}

#
table(segments$CNA_status)
#CN_Gain CN_Loss Diploid     UPD 
#   2866    2934   15089    1542 

#
st2 <- segments



###### ####### ######
### 2- UNIFY the 3 previous dataframes, so that we can JOIN all samples for the posterior CIRCOS-plots:
peris2 <- peris[,c("chrom","start","end","log2_rec","mafR","cf.em","tcn.em","lcn.em","sample","Lesion","CNA_status")]

aecc2 <- aecc[,c("chrom","start","end","log2_rec","mafR","ccf","tcn.em","lcn.em","sample","Lesion","CNA_status")]
colnames(aecc2)[6] <- "cf.em"
aecc2$chrom <- gsub("chr","",aecc2$chrom)

st22 <- st2[,c("chrom","start","end","log2_rec","mean_BAF","sample","CNA_status")]
colnames(st22)[5] <- "mafR"
st22$Lesion <- "cancer"
st22$cf.em <- NA
st22$tcn.em <- NA
st22$lcn.em <- NA
st22 <- st22[,c("chrom","start","end","log2_rec","mafR","cf.em","tcn.em","lcn.em","sample","Lesion","CNA_status")]
st22$chrom <- gsub("chr","",st22$chrom)


### All COHORT-SEGMENTS:
segments.all <- rbind(peris2,aecc2)
segments.all <- rbind(segments.all,st22)
dim(segments.all)
#[1] 38309    11

# Segments of CANCER Lesions:
segments.aes <- segments.all[segments.all$Lesion %in% c("cancer","Tumor"),]
segments.aes$Lesion <- "cancer"
dim(segments.aes)
#[1] 28655    11

# Segments of ADENOMA Lesions:
segments.ads <- segments.all[segments.all$Lesion %in% c("adenoma","Polyp"),]
segments.ads$Lesion <- "adenoma"
dim(segments.ads)
#[1] 6499   11


### By assembling the 3 cohorts (PERIS+AECC+Stage-II), we accumulate N=146 CANCER-Lesions:
length(table(segments.aes$sample))
#[1] 146
### By assembling the 3 cohorts (PERIS+AECC+Stage-II), we accumulate N=73 ADENOMA-Lesions:
length(table(segments.ads$sample))
#[1] 73



###### ###### ###### ######
### 3- Import POSITIONS for circlize-PLOTS:
###### ###### ###### ######
pos.circ <- read.table("~/Desktop/Analysis/Prediction_Tumor-stage_Integrative_2022/Circos_CNAs/Sorted_bedtable_hg19.txt",
                       header=TRUE,stringsAsFactors=FALSE)
dim(pos.circ)
#[1] 28650     5

## Positions for Cancers (N=146):
pos.aes <- data.frame(matrix(nrow=nrow(pos.circ),ncol=5+146))
colnames(pos.aes) <- c(colnames(pos.circ),names(table(segments.aes$sample))) 
dim(pos.aes)
#[1] 28650    151
pos.aes$chrom <- pos.circ$chrom
pos.aes$start <- pos.circ$start
pos.aes$end <- pos.circ$end
pos.aes$genename <- pos.circ$genename
pos.aes$n <- pos.circ$n

# For CN_Gains:
pos.aes.cng <- pos.aes
# For CN_Losses:
pos.aes.cnl <- pos.aes
# For UPD's:
pos.aes.upd <- pos.aes


## Positions for Adenomas (N=73):
pos.ads <- data.frame(matrix(nrow=nrow(pos.circ),ncol=5+73))
colnames(pos.ads) <- c(colnames(pos.circ),names(table(segments.ads$sample))) 
dim(pos.ads)
#[1] 28650    78
pos.ads$chrom <- pos.circ$chrom
pos.ads$start <- pos.circ$start
pos.ads$end <- pos.circ$end
pos.ads$genename <- pos.circ$genename
pos.ads$n <- pos.circ$n

# For CN_Gains:
pos.ads.cng <- pos.ads
# For CN_Losses:
pos.ads.cnl <- pos.ads
# For UPD's:
pos.ads.upd <- pos.ads



###### ###### ###### ######
### 4- FILL the previous 2x3 dataframes (i.e., "aes.cng", "ads.cng", "aes.cnl", etc.) with 1's & 0's in
### each position, according whether they present or NOT each CNA/UPD-type:
###### ###### ###### ######
#
segments <- rbind(segments.aes,segments.ads)
dim(segments)
#[1] 35154    11

#
for(i in 1:nrow(segments)){
  
  chrom <- segments$chrom[i]
  start <- segments$start[i]
  end <- segments$end[i]
  sample.id <- segments$sample[i]
  
  # Matrix "pos.aes.cng" -- For CN_Gains in Cancers (N=146):
  if(segments$CNA_status[i]=="CN_Gain" & segments$Lesion[i]=="cancer") {
    pos.aes.cng[pos.aes.cng$chrom==chrom & pos.aes.cng$start >= start & pos.aes.cng$end <= end,
                colnames(pos.aes.cng)==sample.id] <- 1 }
  else { pos.aes.cng[pos.aes.cng$chrom==chrom & pos.aes.cng$start >= start & pos.aes.cng$end <= end,
                     colnames(pos.aes.cng)==sample.id] <- 0 }
  
  # Matrix "pos.aes.cnl" -- For CN_Losses in Cancers (N=146):
  if(segments$CNA_status[i]=="CN_Loss" & segments$Lesion[i]=="cancer") {
    pos.aes.cnl[pos.aes.cnl$chrom==chrom & pos.aes.cnl$start >= start & pos.aes.cnl$end <= end,
                colnames(pos.aes.cnl)==sample.id] <- 1 }
  else { pos.aes.cnl[pos.aes.cnl$chrom==chrom & pos.aes.cnl$start >= start & pos.aes.cnl$end <= end,
                     colnames(pos.aes.cnl)==sample.id] <- 0 }
  
  # Matrix "pos.aes.upd" -- For UPD's in Cancers (N=146):
  if(segments$CNA_status[i]=="UPD" & segments$Lesion[i]=="cancer") {
    pos.aes.upd[pos.aes.upd$chrom==chrom & pos.aes.upd$start >= start & pos.aes.upd$end <= end,
                colnames(pos.aes.upd)==sample.id] <- 1 }
  else { pos.aes.upd[pos.aes.upd$chrom==chrom & pos.aes.upd$start >= start & pos.aes.upd$end <= end,
                     colnames(pos.aes.upd)==sample.id] <- 0 }
  
  
  # Matrix "pos.ads.cng" -- For CN_Gains in Adenomas (N=73):
  if(segments$CNA_status[i]=="CN_Gain" & segments$Lesion[i]=="adenoma") {
    pos.ads.cng[pos.ads.cng$chrom==chrom & pos.ads.cng$start >= start & pos.ads.cng$end <= end,
                colnames(pos.ads.cng)==sample.id] <- 1 }
  else { pos.ads.cng[pos.ads.cng$chrom==chrom & pos.ads.cng$start >= start & pos.ads.cng$end <= end,
                     colnames(pos.ads.cng)==sample.id] <- 0 }
  
  # Matrix "pos.ads.cnl" -- For CN_Losses in Adenomas (N=73):
  if(segments$CNA_status[i]=="CN_Loss" & segments$Lesion[i]=="adenoma") {
    pos.ads.cnl[pos.ads.cnl$chrom==chrom & pos.ads.cnl$start >= start & pos.ads.cnl$end <= end,
                colnames(pos.ads.cnl)==sample.id] <- 1 }
  else { pos.ads.cnl[pos.ads.cnl$chrom==chrom & pos.ads.cnl$start >= start & pos.ads.cnl$end <= end,
                     colnames(pos.ads.cnl)==sample.id] <- 0 }
  
  # Matrix "pos.ads.upd" -- For UPD's in Adenomas (N=73):
  if(segments$CNA_status[i]=="UPD" & segments$Lesion[i]=="adenoma") {
    pos.ads.upd[pos.ads.upd$chrom==chrom & pos.ads.upd$start >= start & pos.ads.upd$end <= end,
                colnames(pos.ads.upd)==sample.id] <- 1 }
  else { pos.ads.upd[pos.ads.upd$chrom==chrom & pos.ads.upd$start >= start & pos.ads.upd$end <= end,
                     colnames(pos.ads.upd)==sample.id] <- 0 }
  
  print(c(i,sample.id))
}


### FILL all the blank/NA-gaps with 0, EXCEPT if the rest of the chromosome is GAINED:
pos.aes.cng[is.na(pos.aes.cng)] <- 0
pos.aes.cnl[is.na(pos.aes.cnl)] <- 0
pos.aes.upd[is.na(pos.aes.upd)] <- 0

pos.ads.cng[is.na(pos.ads.cng)] <- 0
pos.ads.cnl[is.na(pos.ads.cnl)] <- 0
pos.ads.upd[is.na(pos.ads.upd)] <- 0


### CHECKPOINTS for the next-command LOOP:
setwd("~/Desktop/Analysis/Prediction_Tumor-stage_Integrative_2022/Circos_CNAs")

#write.table(pos.aes.cng,"Cancers_Gains_WES_June2022.txt",col.names=T,row.names=F,quote=F,sep="\t")
#write.table(pos.aes.cnl,"Cancers_Losses_WES_June2022.txt",col.names=T,row.names=F,quote=F,sep="\t")
#write.table(pos.aes.upd,"Cancers_UPDs_WES_June2022.txt",col.names=T,row.names=F,quote=F,sep="\t")
#write.table(pos.ads.cng,"Adenomas_Gains_WES_June2022.txt",col.names=T,row.names=F,quote=F,sep="\t")
#write.table(pos.ads.cnl,"Adenomas_Losses_WES_June2022.txt",col.names=T,row.names=F,quote=F,sep="\t")
#write.table(pos.ads.upd,"Adenomas_UPDs_WES_June2022.txt",col.names=T,row.names=F,quote=F,sep="\t")



###### ###### ###### ######
### 4.5- ONLY SELECT those Samples being: Adenomas (in case of 
### df's "pos.ads") or Stage I-II cancers (in case of df's "pos.aes"):
###### ###### ###### ######

### CHECKPOINTS:
pos.aes.cng0 <- pos.aes.cng
pos.aes.cnl0 <- pos.aes.cnl
pos.aes.upd0 <- pos.aes.upd
pos.ads.cng0 <- pos.ads.cng
pos.ads.cnl0 <- pos.ads.cnl
pos.ads.upd0 <- pos.ads.upd

### Dataframes with the FINAL included Patients: "polyps.all", "stage1.all", "stage2.all"
cancer.names <- c(stage1.all$Sample,stage2.all$Sample)
polyp.names <- c(polyps.all$Sample)

# Cancers:
pos.aes.cng2 <- cbind(pos.aes.cng[,1:5],
                      pos.aes.cng[,colnames(pos.aes.cng) %in% cancer.names])
pos.aes.cnl2 <- cbind(pos.aes.cnl[,1:5],
                      pos.aes.cnl[,colnames(pos.aes.cnl) %in% cancer.names])
pos.aes.upd2 <- cbind(pos.aes.upd[,1:5],
                      pos.aes.upd[,colnames(pos.aes.upd) %in% cancer.names])
dim(pos.aes.cng2)
#[1] 28650   137
dim(pos.aes.cnl2)
#[1] 28650   137
dim(pos.aes.upd2)
#[1] 28650   137

# Polyps:
pos.ads.cng2 <- cbind(pos.ads.cng[,1:5],
                      pos.ads.cng[,colnames(pos.ads.cng) %in% polyp.names])
pos.ads.cnl2 <- cbind(pos.ads.cnl[,1:5],
                      pos.ads.cnl[,colnames(pos.ads.cnl) %in% polyp.names])
pos.ads.upd2 <- cbind(pos.ads.upd[,1:5],
                      pos.ads.upd[,colnames(pos.ads.upd) %in% polyp.names])
dim(pos.ads.cng2)
#[1] 28650    78
dim(pos.ads.cnl2)
#[1] 28650    78
dim(pos.ads.upd2)
#[1] 28650    78


### CHECKPOINTS for the next-command LOOP:
setwd("~/Desktop/Analysis/Prediction_Tumor-stage_Integrative_2022/Circos_CNAs")

#write.table(pos.aes.cng2,"Cancers_Gains_WES_June2022_v2.txt",col.names=T,row.names=F,quote=F,sep="\t")
#write.table(pos.aes.cnl2,"Cancers_Losses_WES_June2022_v2.txt",col.names=T,row.names=F,quote=F,sep="\t")
#write.table(pos.aes.upd2,"Cancers_UPDs_WES_June2022_v2.txt",col.names=T,row.names=F,quote=F,sep="\t")
#write.table(pos.ads.cng2,"Adenomas_Gains_WES_June2022_v2.txt",col.names=T,row.names=F,quote=F,sep="\t")
#write.table(pos.ads.cnl2,"Adenomas_Losses_WES_June2022_v2.txt",col.names=T,row.names=F,quote=F,sep="\t")
#write.table(pos.ads.upd2,"Adenomas_UPDs_WES_June2022_v2.txt",col.names=T,row.names=F,quote=F,sep="\t")



###### ###### ###### ######
### 5- Create the CIRCOS-PLOTS using "circlize":
###### ###### ###### ######

### Define Lesion-CLASSES:
cancerSet <- c("Cancers","Adenomas")
cancerName <- c("Cancers \n stages I-II","Adenomas")


### LOOP:
for (k in 1:length(cancerSet)) {
  
  cancer <- cancerSet[k]
  
  
  ###### ###### ######
  ### (5.1) Define a CHROMOSOME-matrix from the internal-data in "circlize" package:
  d = read.table(file = paste0(system.file(package = "circlize"), "/extdata/cytoBand.txt"),
                 colClasses = c("character", "numeric", "numeric", "character", "character"))
  
  d <- d[which(d[,1]%in%paste("chr",1:22,sep="")),]
  chromosome = unique(d[[1]])
  chromosome.ind = gsub("chr", "", chromosome)
  chromosome.num = grep("^\\d+$", chromosome.ind, value = TRUE)
  chromosome.letter = chromosome.ind[!grepl("^\\d+$", chromosome.ind)]
  chromosome.num = sort(as.numeric(chromosome.num))
  chromosome.letter = sort(chromosome.letter)
  chromosome.num = paste("chr", chromosome.num, sep = "")
  chromosome.letter = paste("chr", chromosome.letter, sep = "")
  
  chromosome = c(chromosome.num, chromosome.letter)
  chromosome <- chromosome[1:22]
  
  xlim = matrix(nrow = 0, ncol = 2)
  for(chr in chromosome) {
    d2 = d[d[[1]] == chr, ]
    xlim = rbind(xlim, c(min(d2[[2]]), max(d2[[3]])))
  }
  
  
  ###### ###### ######
  ### (5.2) Define the IDEOGRAM (= circular-PLOT):
  setwd("~/Desktop/Analysis/Prediction_Tumor-stage_Integrative_2022/Circos_CNAs")
  
  png(paste("Circos_",cancer,"_June2022_v2.png",sep=""), 1800, 1800, res=200);
  
  circos.par(points.overflow.warning = FALSE)
  par(lwd = 0.5)
  circos.par(cell.padding = c(0, 0, 0, 0))
  circos.par("start.degree" = 90)
  circos.initialize(factors = factor(chromosome, levels = chromosome), xlim = xlim)
  #circos.initializeWithIdeogram(species = "hg19")
  
  
  ###### ###### ######
  ### (5.3) Define REGIONS + COLORS within the IDEOGRAM:
  circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.08)
  for(chr in chromosome) {
    # data in current `chr`
    d2 = d[d[[1]] == chr, ]
    n = nrow(d2)
    
    # Assign colors to Cytobands within the Ideogram:
    col = rep("#FFFFFF", n)
    col[d2[[5]] == "gpos100"] = rgb(0, 0, 0, maxColorValue = 255)
    col[d2[[5]] == "gpos"]    = rgb(0, 0, 0, maxColorValue = 255)
    col[d2[[5]] == "gpos75"]  = rgb(130, 130, 130, maxColorValue = 255)
    col[d2[[5]] == "gpos66"]  = rgb(160, 160, 160, maxColorValue = 255)
    col[d2[[5]] == "gpos50"]  = rgb(200, 200, 200, maxColorValue = 255)
    col[d2[[5]] == "gpos33"]  = rgb(210, 210, 210, maxColorValue = 255)
    col[d2[[5]] == "gpos25"]  = rgb(200, 200, 200, maxColorValue = 255)
    col[d2[[5]] == "gvar"]    = rgb(220, 220, 220, maxColorValue = 255)
    col[d2[[5]] == "gneg"]    = rgb(255, 255, 255, maxColorValue = 255)
    col[d2[[5]] == "acen"]    = rgb(217, 47, 39, maxColorValue = 255)
    col[d2[[5]] == "stalk"]   = rgb(100, 127, 164, maxColorValue = 255)
    
    # Rectangles for different locus:
    for(i in seq_len(n)) {
      circos.rect(d2[i, 2], 0, d2[i, 3], 0.4, sector.index = chr,
                  col = col[i], border = NA)
    }
    
    # Rectangle that covers the whole chromosome:
    circos.rect(d2[1, 2], 0, d2[n, 3], 0.4, sector.index = chr, border = "black")
    
    chr.xlim = get.cell.meta.data("xlim", sector.index = chr)
    
    # Chromosome names, only the number part or the letter part:
    circos.text(mean(chr.xlim), 0.8, labels = gsub("chr", "", chr),
                sector.index = chr, cex = 1.0)
  }
  
  
  ###### ###### ######
  ### (5.4) UPD / cnn-LOH:
  CNV <- read.csv(paste("~/Desktop/Analysis/Prediction_Tumor-stage_Integrative_2022/Circos_CNAs/",cancer,"_UPDs_WES_June2022_v2.txt",sep=""), header=T, sep='\t')
  #CNV <- pos.aes.upd
  CNV <- CNV[which(CNV[,"chrom"]<=22),]
  pos <- CNV$start  #(CNV$end + CNV$start)/2
  chrom <- CNV$chrom
  perc <- rowMeans( as.matrix(CNV[,6:(ncol(CNV))]) )
  circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.17)    # He canviat el ylim=(c(0,0.8) a ylim=c(0,1)
  for(chr in chromosome) {
    # data in current `chr`
    d2 = d[d[[1]] == chr, ]
    n = nrow(d2)
    
    pos1  <-  c(0, pos[chrom==gsub("chr", "", chr)][1], pos[chrom==gsub("chr", "", chr)],
                seq(d2[n, 3], 0, length.out=100) )
    
    perc1 <- c(0, 0, perc[chrom==gsub("chr", "", chr)], rep(0,100))
    
    circos.lines(c(0,d2[n, 3]), c(0.1,0.1), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.2,0.2), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.3,0.3), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.4,0.4), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.5,0.5), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.6,0.6), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.7,0.7), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.8,0.8), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.9,0.9), sector.index = chr, col="grey", lty=2 )
    
    circos.polygon(pos1, perc1, sector.index = chr,
                   col=rgb(0/256, 100/256, 0/256, 1), border=NA)   #He canviat el 0.6 per 1!!!
    
    circos.rect(d2[1, 2], 0, d2[n, 3], 1, sector.index = chr, border = "black")
    
  }
  
  
  ###### ###### ######
  ### (5.5) LOSSES:
  CNV <- read.csv(paste("~/Desktop/Analysis/Prediction_Tumor-stage_Integrative_2022/Circos_CNAs/",cancer,"_Losses_WES_June2022_v2.txt",sep=""), header=T, sep='\t')
  #CNV <- pos.aes.cnl
  CNV <- CNV[which(CNV[,"chrom"]<=22),]
  pos <- CNV$start  #(CNV$end + CNV$start)/2
  chrom <- CNV$chrom
  perc<- rowMeans( as.matrix(CNV[,6:(ncol(CNV))]) )
  circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.17)    # He canviat el ylim=(c(0,0.8) a ylim=c(0,1)
  for(chr in chromosome) {
    # data in current `chr`
    d2 = d[d[[1]] == chr, ]
    n = nrow(d2)
    
    pos1  <-  c(0, pos[chrom==gsub("chr", "", chr)][1], pos[chrom==gsub("chr", "", chr)],
                seq(d2[n, 3], 0, length.out=100) )
    
    perc1 <- c(0, 0, perc[chrom==gsub("chr", "", chr)], rep(0,100))
    
    circos.lines(c(0,d2[n, 3]), c(0.1,0.1), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.2,0.2), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.3,0.3), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.4,0.4), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.5,0.5), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.6,0.6), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.7,0.7), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.8,0.8), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.9,0.9), sector.index = chr, col="grey", lty=2 )
    
    circos.polygon(pos1, perc1, sector.index = chr,
                   col=rgb(0/256, 0/256, 205/256, 1), border=NA)   #He canviat el 0.6 per 1!!!
    
    circos.rect(d2[1, 2], 0, d2[n, 3], 1, sector.index = chr, border = "black")
    
  }
  
  
  ###### ###### ######
  ### (5.6) GAINS:
  CNV <- read.csv(paste("~/Desktop/Analysis/Prediction_Tumor-stage_Integrative_2022/Circos_CNAs/",cancer,"_Gains_WES_June2022_v2.txt",sep=""), header=T, sep='\t')
  #CNV <- pos.aes.cng
  CNV <- CNV[which(CNV[,"chrom"]<=22),]
  pos <- CNV$start # (CNV$end + CNV$start)/2
  chrom <- CNV$chrom
  perc<- rowMeans( as.matrix(CNV[,6:(ncol(CNV))]) )
  circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.17)      # He canviat el ylim=(c(0,0.8) a ylim=c(0,1)
  for(chr in chromosome) {
    # data in current `chr`
    d2 = d[d[[1]] == chr, ]
    n = nrow(d2)
    
    pos1  <-  c(0, pos[chrom==gsub("chr", "", chr)][1], pos[chrom==gsub("chr", "", chr)],
                seq(d2[n, 3], 0, length.out=100) )
    
    perc1 <- c(0, 0, perc[chrom==gsub("chr", "", chr)], rep(0,100))
    
    circos.lines(c(0,d2[n, 3]), c(0.1,0.1), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.2,0.2), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.3,0.3), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.4,0.4), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.5,0.5), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.6,0.6), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.7,0.7), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.8,0.8), sector.index = chr, col="grey", lty=2 )
    circos.lines(c(0,d2[n, 3]), c(0.9,0.9), sector.index = chr, col="grey", lty=2 )
    
    circos.polygon(pos1, perc1, sector.index = chr,
                   col=rgb(205/256, 0/256, 0/256, 1), border=NA)   #He canviat el 0.6 per 1!!!
    
    circos.rect(d2[1, 2], 0, d2[n, 3], 1, sector.index = chr, border = "black")
    
  }
  
  text(0,0.0,cancerName[k],cex=1.8)
  #text(0,-0.02,paste(nTSG[k]," known TSG",sep=""),cex=1.2)
  #text(0,-0.057,paste("pval=",round(pvalg, 5),sep=""),cex=1.2)
  
  dev.off()
  
  #circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA)
  #circos.text(88888888, 0.2, sector.index = "chr6", labels = "site", adj = c(0.5, 1))
  #circos.lines(c(88888888, 88888888), c(0.3, 1), sector.index = "chr6", straight = TRUE)
  circos.clear()
  
}


