##############################################
#Load data, separate tissues
##############################################

library(minfi)
library(sva)

setwd('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/')

#load noob data
load("/dcl01/NDEpi/data/Projects/InProgress/kbakulsk/kbakulsk/EARLI/450k-round2/Noob-beta-All-Samples-EARLI-both-rounds.rda")

#load pd
load("/dcl01/NDEpi/data/MasterCohortData/EARLI/EARLIBiologicalSamples/450k/EARLI_450k_all1155samples_pd_with_ancestry.rda")
load("/dcl01/NDEpi/data/MasterCohortData/EARLI/EARLIBiologicalSamples/450k/earli.450k.pd.N=1225.20161117.rda")

#drop poor detection probes
#load('/dcl01/NDEpi/data/Projects/InProgress/kbakulsk/kbakulsk/EARLI/450k-round2/Failed-Probes-By-DetP-01.rda')
load("/dcl01/NDEpi/data/Projects/InProgress/taung/dnam_asd/dnam_asd_proj/Failed-Probes-By-DetP-01.rda")
dim(probe.fail) #635 failed probes
noob.beta <- noob.beta[!(rownames(noob.beta)%in%rownames(probe.fail)),]
rm(probe.fail)

#drop cross reactive
cross <- read.csv("/dcl01/NDEpi/data/Projects/InProgress/jdou/48639-non-specific-probes-Illumina450k.csv", stringsAsFactors=FALSE)
noob.beta <- noob.beta[!rownames(noob.beta)%in%cross$TargetID,]
rm(cross)

#save beta file with dropped probes
save(noob.beta, file='noob_beta_drop-poor-det_drop-cross.rda')
load('noob_beta_drop-poor-det_drop-cross.rda')

#get the cord blood
pd.cord <- pd[pd$Tissue=='Cord.Blood',]
pd.cord <- merge(pd.cord,pd.ancest[,c('Tissue.Cat','predictedSex')],by='row.names',sort=FALSE)
rownames(pd.cord) <- pd.cord$Row.names
#drop those with discordant predicted and observed sex
discordant.sex <- pd.cord[!is.na(pd.cord$Sex) & pd.cord$Sex != pd.cord$predictedSex,]
#save(discordant.sex,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/sex_mismatched_samples.rda')
#load('/dcl01/NDEpi/data/Projects/InProgress/jdou/sex_mismatched_samples.rda')
pd.cord <- pd.cord[!rownames(pd.cord)%in%rownames(discordant.sex),]
#cell type
load("/dcl01/NDEpi/data/MasterCohortData/EARLI/EARLIBiologicalSamples/450k/EARLI_cord_cell_type_7_estimates_cord_reference.rda")
counts <- data.frame(counts)
pd.cord <- merge(pd.cord[,-1],counts,by='row.names',sort=F,all.x=T)
rownames(pd.cord) <- pd.cord$Row.names

#get the maternal blood
pd.maternal <- pd[pd$Tissue=='blood' & pd$DCC_Subject_Type=='Biological Mother',]
pd.earlypreg <- pd.maternal[pd.maternal$DCC_Visit=='Enrollment Clinic Visit',]
pd.latepreg <- pd.maternal[startsWith(pd.maternal$DCC_Visit,'Home Visit'),]
#cell type
load("/dcl01/NDEpi/data/Projects/InProgress/kbakulsk/kbakulsk/EARLI/450k-round2/Blood-cell-type-estimates.rda")
pd.earlypreg <- merge(pd.earlypreg,cell,by='row.names',sort=F,all.X=T)
pd.latepreg <- merge(pd.latepreg,cell,by='row.names',sort=F,all.X=T)
rownames(pd.earlypreg) <- pd.earlypreg$Row.names
rownames(pd.latepreg) <- pd.latepreg$Row.names
#T2 and T3 overlap, choose T3
fams <- data.frame(table(pd.latepreg$DCC_Family_ID))
fams <- fams[fams$Freq==2,]
remove <- which(pd.latepreg$DCC_Family_ID %in% fams[,1] & pd.latepreg$DCC_Visit=='Home Visit (2nd Trimster)')
pd.latepreg <- pd.latepreg[-remove,]
rm(remove)

#placenta
pd.placenta <- pd[pd$Tissue=='placenta',]
pd.placenta$predictedSex <- pd.ancest[rownames(pd.placenta),'predictedSex']
pd.plc <- pd.placenta[pd.placenta$Phenotype=='placenta.child',]
pd.plm <- pd.placenta[pd.placenta$Phenotype=='placenta.mom',]
#remove multiple births
multi <- pd.plc[duplicated(pd.plc$DCC_Family_ID),"DCC_Family_ID"]
pd.plc <- pd.plc[!pd.plc$DCC_Family_ID %in% multi,]
pd.plm <- pd.plm[!pd.plm$DCC_Family_ID %in% multi,]


#36 month SRS score
SRS36 <- read.csv("/dcl01/NDEpi/data/MasterCohortData/EARLI/EARLIBehavioralAssessments/Social Responsiveness Scale (SRS)/srs36.csv")
reducedSRS36 <- SRS36[,c('Subject_ID', 'totalrawscore', 'overall_t_score')]
table(SRS36$Subject_ID %in% pd.cord$Subject.ID)
# FALSE  TRUE 
# 50    96 

#12 month AOSI score
AOSI12 <- read.csv("/dcl01/NDEpi/data/MasterCohortData/EARLI/EARLIBehavioralAssessments/AOSI/12moAOSI_01282015.csv")
reducedAOSI12 <- AOSI12[,c('Subject_ID','AOSI_AOSI_TotalScore_Sibling..Postpartum..12.months.')]
table(AOSI12$Subject_ID %in% pd.cord$Subject.ID)
# FALSE  TRUE 
# 66   149 

#ASD Diangosis
ASD <- read.csv("/dcl01/NDEpi/data/MasterCohortData/EARLI/EARLIBehavioralAssessments/36MonthOutcomes/outcome36mos_12_16.csv")
reducedASD <- ASD[,c('Family_ID','ASD_DSM5_Sib_36mos','BSRC_group_Sib_36mos')]
table(pd.earlypreg$DCC_Family_ID %in% ASD$Family_ID)
# FALSE  TRUE
#   43   158
table(pd.latepreg$DCC_Family_ID %in% ASD$Family_ID)
# FALSE  TRUE 
# 19    84 
table(pd.plc$DCC_Family_ID %in% ASD$Family_ID)
# FALSE  TRUE 
# 19   109 
table(pd.plm$DCC_Family_ID %in% ASD$Family_ID)
# FALSE  TRUE 
# 20   109

#merge SRS and AOSI into cord blood phenotype file
pd.cord.bm <- merge(pd.cord,reducedSRS36, by.x='Subject.ID',by.y='Subject_ID',all.x=TRUE,sort=FALSE)
pd.cord.bm <- merge(pd.cord.bm,reducedAOSI12, by.x='Subject.ID',by.y='Subject_ID',all.x=TRUE,sort=FALSE)
rownames(pd.cord.bm) <- pd.cord.bm$Row.names

#merge ASD diagnosis into maternal blood phenotype file
pd.earlypreg.asd <- merge(pd.earlypreg,reducedASD, by.x='DCC_Family_ID',by.y='Family_ID',all.x=T,sort=F )
pd.latepreg.asd <- merge(pd.latepreg,reducedASD, by.x='DCC_Family_ID',by.y='Family_ID',all.x=T,sort=F )
#remove multiple births
pd.earlypreg.asd <- pd.earlypreg.asd[!(duplicated(pd.earlypreg.asd$DCC_Family_ID) | duplicated(pd.earlypreg.asd$DCC_Family_ID,fromLast=T)),]
rownames(pd.earlypreg.asd) <- pd.earlypreg.asd$Row.names
pd.latepreg.asd <- pd.latepreg.asd[!(duplicated(pd.latepreg.asd$DCC_Family_ID) | duplicated(pd.latepreg.asd$DCC_Family_ID,fromLast=T)),]
rownames(pd.latepreg.asd) <- pd.latepreg.asd$Row.names
#baby gender
babygen <- pd.ancest[pd.ancest$Participant=='Sib',c('Individual.ID','FamilyID','predictedSex')]
babygen <- babygen[!duplicated(babygen$Individual.ID),c(2,3)]
pd.earlypreg.asd <- merge(pd.earlypreg.asd,babygen,by.x='DCC_Family_ID',by.y='FamilyID',all.x=T,sort=F)
rownames(pd.earlypreg.asd) <- pd.earlypreg.asd$Row.names
pd.latepreg.asd <- merge(pd.latepreg.asd,babygen,by.x='DCC_Family_ID',by.y='FamilyID',all.x=T,sort=F)
rownames(pd.latepreg.asd) <- pd.latepreg.asd$Row.names
pd.earlypreg.asd <- rename(pd.earlypreg.asd, replace=c('ASD_DSM5_Sib_36mos'='asd5','BSRC_group_Sib_36mos'='td'))
pd.latepreg.asd <- rename(pd.latepreg.asd, replace=c('ASD_DSM5_Sib_36mos'='asd5','BSRC_group_Sib_36mos'='td'))

#merge ASD into placenta phenotype file
pd.plc$Row.names <- rownames(pd.plc)
pd.plm$Row.names <- rownames(pd.plm)
ASD.nomult <- reducedASD[!(duplicated(reducedASD$Family_ID)|duplicated(reducedASD$Family_ID,fromLast=T)),]
plChild <- merge(pd.plc,ASD.nomult,by.x='DCC_Family_ID',by.y='Family_ID',sort=F,all.x=T)
plMom <- merge(pd.plm,ASD.nomult,by.x='DCC_Family_ID',by.y='Family_ID',sort=F,all.x=T)
pd.plc <- plChild
pd.plm <- plMom
pd.plc <- rename(pd.plc,replace=c('ASD_DSM5_Sib_36mos'='asd5','BSRC_group_Sib_36mos'='td'))
pd.plm <- rename(pd.plm,replace=c('ASD_DSM5_Sib_36mos'='asd5','BSRC_group_Sib_36mos'='td'))
rm(ASD.nomult,plChild,plMom)

### files for cord blood and ASD
rASD <- ASD[,c('Subject_ID','ASD_DSM5_Sib_36mos','BSRC_group_Sib_36mos')]
#maternal age
interview <- read.csv("/dcl01/NDEpi/data/MasterCohortData/EARLI/EARLISelfReportedInterviewData/Data Extraction Codebook 20150203_export_20150203_1422987480_wide_20150204_1423070501.csv")
interview.fam <- interview[interview$Family_ID %in% pd.cord$DCC_Family_ID,]
questions <- c('m1_age_q2','m1_fatherage1_q17','m1_education_q25','m1_income_q11')
interview.q <- interview.fam[,c('Family_ID',questions)]
interview.q <- plyr::rename(interview.q,replace=c('m1_age_q2'='Maternal.Age','m1_fatherage1_q17'='Paternal.Age','m1_education_q25'='mom.edu','m1_income_q11'='income'))
interview.q$Maternal.Age <- as.numeric(as.character(interview.q$Maternal.Age))
interview.q$Paternal.Age <- as.numeric(as.character(interview.q$Paternal.Age))
# interview <- read.csv("/dcl01/NDEpi/data/MasterCohortData/EARLI/EARLISelfReportedInterviewData/Data Extraction Codebook 20150203_export_20150203_1422987480.csv", header=T)
# age <- interview[interview$Field_Name=='m1_age_q2',]
# age <- age[,c('Family_ID','Value')]
# names(age)[2] <- 'Maternal.Age'
# age$Maternal.Age <- as.numeric(as.character(age$Maternal.Age))

pd.cord.asd <- merge(pd.cord,rASD,by.x='Subject.ID',by.y='Subject_ID',all.x=T,all.y=F)
pd.cord.asd <- plyr::rename(pd.cord.asd,replace=c('ASD_DSM5_Sib_36mos'='asd5','BSRC_group_Sib_36mos'='td'))
rownames(pd.cord.asd) <- pd.cord.asd$Row.names
pd.cord.asd <- merge(pd.cord.asd,interview.q,by.x='DCC_Family_ID',by.y='Family_ID',all.x=T,all.y=F)
rownames(pd.cord.asd) <- pd.cord.asd$Row.names

#cut for sva
saveRDS(pd.cord.asd,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/pd.cord.rds')
noob.cord <- noob.beta[,rownames(pd.cord)]
saveRDS(noob.cord,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/noob.cord.rds')

#legacy cut, known covariates
pd.cord.asd <- pd.cord.asd[!is.na(pd.cord.asd$asd5) & !is.na(pd.cord.asd$PC1),]
save(pd.cord.asd,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/pd_cord-asd.rda')
noob.beta.cordasd <- noob.beta[,rownames(pd.cord.asd)]
save(noob.beta.cordasd,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/noob_beta_cord-ASD.rda')


### files for cord blood and 36 month SRS
pd.cord.srs <- pd.cord.bm[!is.na(pd.cord.bm$totalrawscore),]
quantile(pd.cord.srs$totalrawscore,c(0,0.25,0.5,0.75,1.0))
# 0%    25%    50%    75%   100% 
# 7.00  17.75  25.50  40.50 133.00
pd.cord.srs$totalraw_Quart <- cut(pd.cord.srs$totalrawscore,c(7,17.75,25.50,40.50,133.0),include.lowest=T)
pd.cord.srs$totalraw_Quart4 <- ifelse(pd.cord.srs$totalraw_Quart=='(40.5,133]',1,0)
pd.cord.srs$log_totalraw <- log(pd.cord.srs$totalrawscore)
noob.beta.cordsrs <- noob.beta[,rownames(pd.cord.srs)]
save(noob.beta.cordsrs,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_36MonthSRS/noob_beta_cord-SRS-20170811.rda')
save(pd.cord.srs,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_36MonthSRS/pd_cord-SRS.rda')

#files for cord blood and 12 month AOSI
pd.cord.aosi <- pd.cord.bm[!is.na(pd.cord.bm$AOSI_AOSI_TotalScore_Sibling..Postpartum..12.months.),]
pd.cord.aosi$AOSI_12Month <- pd.cord.aosi$AOSI_AOSI_TotalScore_Sibling..Postpartum..12.months. #for the sake of readability
pd.cord.aosi <- pd.cord.aosi[!is.na(pd.cord.aosi$PC1),] #drop missing PC data
#categorical quartiles
quantile(pd.cord.aosi$AOSI_12Month,c(0,0.25,0.5,0.75,1.0))
# 0%  25%  50%  75% 100% 
# 0    2    4    8   19 
pd.cord.aosi$AOSI_Quart <- cut(pd.cord.aosi$AOSI_12Month,c(0,2,4,8,19),include.lowest=T)
pd.cord.aosi$AOSI_Quart4 <- ifelse(pd.cord.aosi$AOSI_Quart=='(8,19]',1,0)
pd.cord.aosi$log_aosi <- log(pd.cord.aosi$AOSI_12Month+1)
noob.beta.cordaosi <- noob.beta[,rownames(pd.cord.aosi)]
save(noob.beta.cordaosi,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_12MonthAOSI/noob_beta_cord-AOSI-20170811.rda')
save(pd.cord.aosi,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_12MonthAOSI/pd_cord-AOSI.rda')

### files for maternal early preg blood and ASD diagnosis

#get maternal age
interview <- read.csv("/dcl01/NDEpi/data/MasterCohortData/EARLI/EARLISelfReportedInterviewData/Data Extraction Codebook 20150203_export_20150203_1422987480.csv", header=T)
age <- interview[interview$Field_Name=='m1_age_q2',]
age <- age[,c('Subject_ID','Value')]
names(age)[2] <- 'Maternal.Age'
age$Maternal.Age <- as.numeric(as.character(age$Maternal.Age))
pd.earlypreg.asd <- merge(pd.earlypreg.asd,age,by.x='Subject.ID',by.y='Subject_ID',all.x=T,sort=F)
rownames(pd.earlypreg.asd) <- pd.earlypreg.asd$Row.names
pd.earlypreg.asd$Maternal.Age <- as.numeric(as.character(pd.earlypreg.asd$Maternal.Age))

#cut for sva
saveRDS(pd.earlypreg.asd,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/pd.epreg.rds')
noob.epreg <- noob.beta[,rownames(pd.earlypreg.asd)]
saveRDS(noob.epreg,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/noob.epreg.rds')

#legacy cut
pd.earlypreg.asd <- pd.earlypreg.asd[!is.na(pd.earlypreg.asd$PC1),] #drop missing PC data
noob.beta.epreg <- noob.beta[,rownames(pd.earlypreg.asd)]
save(pd.earlypreg.asd,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/pd_epreg-ASD.rda')
save(noob.beta.epreg,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/noob_beta_momearly-ASD-20170811.rda')

### files for maternal late preg blood and ASD diagnosis
#get maternal age
pd.latepreg.asd <- merge(pd.latepreg.asd,age,by.x='Subject.ID',by.y='Subject_ID',all.x=T,sort=F)
rownames(pd.latepreg.asd) <- pd.latepreg.asd$Row.names
pd.latepreg.asd$Maternal.Age <- as.numeric(as.character(pd.latepreg.asd$Maternal.Age))

#cut for sva
saveRDS(pd.latepreg.asd,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/pd.lpreg.rds')
noob.lpreg <- noob.beta[,rownames(pd.latepreg.asd)]
saveRDS(noob.lpreg,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/noob.lpreg.rds')


#legacy cut
pd.latepreg.asd <- pd.latepreg.asd[!is.na(pd.latepreg.asd$asd5),]
pd.latepreg.asd <- pd.latepreg.asd[!is.na(pd.latepreg.asd$PC1),] #drop missing PC data
noob.beta.lpreg <- noob.beta[,rownames(pd.latepreg.asd)]
save(pd.latepreg.asd,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/pd_lpreg-ASD.rda')
save(noob.beta.lpreg,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/noob_beta_momlate-ASD-20170811.rda')

### files for placenta child and ASD
pd.plc.edit <- pd.plc
pd.plm.edit <- pd.plm
#get maternal age
interview <- read.csv("/dcl01/NDEpi/data/MasterCohortData/EARLI/EARLISelfReportedInterviewData/Data Extraction Codebook 20150203_export_20150203_1422987480.csv", header=T)
age <- interview[interview$Field_Name=='m1_age_q2',]
age <- age[,c('Subject_ID','Value')]
names(age)[2] <- 'Maternal.Age'
age$Maternal.Age <- as.numeric(as.character(age$Maternal.Age))
pd.plc.edit <- merge(pd.plc.edit,age[,c('Family_ID','Maternal.Age')],by.x='DCC_Family_ID',by.y='Family_ID',all.x=T,sort=F)
pd.plm.edit <- merge(pd.plm.edit,age[,c('Subject_ID','Maternal.Age')],by.x='Subject.ID',by.y='Subject_ID',all.x=T,sort=F)

#cut for sva
pd.plc <- pd.plc.edit; rownames(pd.plc) <- pd.plc$Row.names
pd.plm <- pd.plm.edit; rownames(pd.plm) <- pd.plm$Row.names
noob.plc <- noob.beta[,rownames(pd.plc)]
noob.plm <- noob.beta[,rownames(pd.plm)]
saveRDS(pd.plm,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/pd.plm.rds")
saveRDS(noob.plm,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/noob.plm.rds")
saveRDS(pd.plc,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/pd.plc.rds")
saveRDS(noob.plc,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/noob.plc.rds")



##############################################
#Table 1
##############################################
library(doBy)
library(dplyr)

pd.cord <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/pd.cord.rds")
pd <- pd.cord[!is.na(pd.cord$td) & pd.cord$td != "",]
pd$td <- factor(pd$td,levels=c('TD','Non-TD','ASD'))

table(pd$DCC_Locale)
prop.table(table(pd$DCC_Locale))
tbl <- table(pd$DCC_Locale,pd$td)
tbl
prop.table(tbl,2)
chisq.test(tbl)


table(pd$predictedSex)
prop.table(table(pd$predictedSex))
tbl <- table(pd$predictedSex,pd$td)
tbl
prop.table(tbl,2)
chisq.test(tbl)

mean(pd$Maternal.Age)
sd(pd$Maternal.Age)
summaryBy(Maternal.Age~td,data=pd)
summaryBy(Maternal.Age~td,data=pd,FUN=sd)
summary(aov(Maternal.Age~td,data=pd))

pd$age.q4 <- ntile(pd$Maternal.Age,4)
table(pd$age.q4)
prop.table(table(pd$age.q4))
tbl <- table(pd$age.q4,pd$td)
tbl
prop.table(tbl,2)
chisq.test(tbl)

mean(pd$Paternal.Age,na.rm=T)
sd(pd$Paternal.Age,na.rm=T)
summaryBy(Paternal.Age~td,data=pd,na.rm=T)
summaryBy(Paternal.Age~td,data=pd,FUN=sd,na.rm=T)
summary(aov(Paternal.Age~td,data=pd))

pd$agef.q4 <- ntile(pd$Paternal.Age,4)
table(pd$agef.q4,exclude=NULL)
prop.table(table(pd$agef.q4,exclude=NULL))
tbl <- table(pd$agef.q4,pd$td,exclude=NULL)
tbl
prop.table(tbl,2)
chisq.test(table(pd$agef.q4,pd$td))

table(pd$mom.edu,exclude=NULL)
pd$mom.edu4 <- ifelse(pd$mom.edu %in% c(4,6,9),'LE_high',
                      ifelse(pd$mom.edu %in% c(10,11,12,13),'some_asso',
                             ifelse(pd$mom.edu == 14, 'bach', 
                                    ifelse(pd$mom.edu %in% c(15,16,17),'grad',NA))))
pd$mom.edu4 <- factor(pd$mom.edu4,levels=c('LE_high','some_asso','bach','grad'))
table(pd$mom.edu4,exclude=NULL)
prop.table(table(pd$mom.edu4,exclude=NULL))
tbl <- table(pd$mom.edu4,pd$td,exclude=NULL)
tbl
prop.table(tbl,2)
chisq.test(table(pd$mom.edu4,pd$td))
fisher.test(table(pd$mom.edu4,pd$td))

pd$popid.questionnaire <- factor(pd$popid.questionnaire, levels=c('EARLI Eur','EARLI Afr','EARLI Latin','EARLI Mix'))
table(pd$popid.questionnaire,exclude=NULL)
prop.table(table(pd$popid.questionnaire,exclude=NULL))
tbl <- table(pd$popid.questionnaire,pd$td,exclude=NULL)
tbl
prop.table(tbl,2)

table(pd$income)
pd$inc3 <- ifelse(pd$income %in% c(1:4),'lt50k',
                  ifelse(pd$income %in% c(5,6),'ge50klt100k',
                         ifelse(pd$income %in% c(7,8),'gt100k',NA)))
pd$inc3 <- factor(pd$inc3,levels=c('lt50k','ge50klt100k','gt100k'))
table(pd$inc3,exclude=NULL)
prop.table(table(pd$inc3,exclude=NULL))
tbl <- table(pd$inc3,pd$td,exclude=NULL)
tbl
prop.table(tbl,2)
chisq.test(table(pd$inc3,pd$td))

table(pd$Round,exclude=NULL)
prop.table(table(pd$Round,exclude=NULL))
tbl <- table(pd$Round,pd$td)
tbl
prop.table(tbl,2)
chisq.test(tbl)

summary(pd$Gran); sd(pd$Gran)
summaryBy(Gran~td,data=pd); summaryBy(Gran~td,data=pd,FUN=sd)
summary(aov(Gran~td,data=pd))

summary(pd$CD8T); sd(pd$CD8T)
summaryBy(CD8T~td,data=pd); summaryBy(CD8T~td,data=pd,FUN=sd)
summary(aov(CD8T~td,data=pd))

summary(pd$CD4T); sd(pd$CD4T)
summaryBy(CD4T~td,data=pd); summaryBy(CD4T~td,data=pd,FUN=sd)
summary(aov(CD4T~td,data=pd))

summary(pd$NK*100); sd(pd$NK*100)
summaryBy(NK*100~td,data=pd); summaryBy(NK*100~td,data=pd,FUN=sd)
summary(aov(NK~td,data=pd))

summary(pd$Bcell*100); sd(pd$Bcell*100)
summaryBy(Bcell*100~td,data=pd); summaryBy(Bcell*100~td,data=pd,FUN=sd)
summary(aov(Bcell~td,data=pd))

summary(pd$Mono*100); sd(pd$Mono*100)
summaryBy(Mono*100~td,data=pd); summaryBy(Mono*100~td,data=pd,FUN=sd)
summary(aov(Mono~td,data=pd))

summary(pd$nRBC*100); sd(pd$nRBC*100)
summaryBy(nRBC*100~td,data=pd); summaryBy(nRBC*100~td,data=pd,FUN=sd)
summary(aov(nRBC~td,data=pd))


##############################################
#Surrogate Variables
##############################################
library(sva)
library(limma)
library(qqman)
library(corrplot)
source("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/ASD_sens_functions.R")


### cord
pd.cord <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/pd.cord.rds")
beta.cord <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/noob.cord.rds")

#asd vs td
FOLDER <- '/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/asd-v-td-sv/'
pd.cord.td <- pd.cord[!is.na(pd.cord$td),]
pd.cord.td <- pd.cord.td[pd.cord.td$td=='ASD' | pd.cord.td$td=='TD',]
pd.cord.td$td <- as.character(pd.cord.td$td)
pd.cord.td$td <- factor(pd.cord.td$td,levels=c('TD','ASD'))
table(pd.cord.td$td)
# TD ASD
# 47  29

pdsv.cord.td <- sva.fit(beta=beta.cord[,rownames(pd.cord.td)],pd=pd.cord.td,var='td')
saveRDS(pdsv.cord.td,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/pdsv.cord.td.rds")

pdsv.cord.td <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/pdsv.cord.td.rds")
bring.the.heatmap.sv(pdsv.cord.td,covar=c('Round','predictedSex','PC1','CD8T','CD4T','NK','Bcell','Mono','Gran','nRBC','Maternal.Age','asd5'),folder=FOLDER)

lambdas.cord <- sva.sequential(beta=beta.cord[,rownames(pdsv.cord.td)],pd=pdsv.cord.td,var='td',folder=FOLDER)
save(lambdas.cord,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/lambdas.cord.rda")
rm(pd.cord,pd.cord.td,pdsv.cord.td,lambdas.cord,beta.cord)

#ntd vs td
FOLDER <- '/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/ntd-v-td-sv/'
pd.cord.ntd <- pd.cord[!is.na(pd.cord$td),]
pd.cord.ntd <- pd.cord.ntd[pd.cord.ntd$td=='Non-TD' | pd.cord.ntd$td=='TD',]
pd.cord.ntd$td <- as.character(pd.cord.ntd$td)
pd.cord.ntd$td <- factor(pd.cord.ntd$td,levels=c('TD','Non-TD'))
table(pd.cord.ntd$td)
# TD Non-TD
# 47     57

pdsv.cord.ntd <- sva.fit(beta=beta.cord[,rownames(pd.cord.ntd)],pd=pd.cord.ntd,var='td')
saveRDS(pdsv.cord.ntd,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/pdsv.cord.ntd.rds")

pdsv.cord.td <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/pdsv.cord.ntd.rds")
bring.the.heatmap.sv(pdsv.cord.ntd,covar=c('Round','predictedSex','PC1','CD8T','CD4T','NK','Bcell','Mono','Gran','nRBC','Maternal.Age','td'),folder=FOLDER)

lambdas.cord <- sva.sequential(beta=beta.cord[,rownames(pdsv.cord.ntd)],pd=pdsv.cord.ntd,var='td',folder=FOLDER)
save(lambdas.cord,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/lambdas.cord.ntd.rda")

rm(pd.cord,pd.cord.td,pdsv.cord.td,lambdas.cord,beta.cord,pd.cord.ntd,pdsv.cord.ntd,lambdas.cord.ntd)


### early pregnancy
pd.epreg <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/pd.epreg.rds")
beta.epreg <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/noob.epreg.rds")

#asd vs td
FOLDER <- "/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/asd-v-td-sv/"
pd.epreg.td <- pd.epreg[!is.na(pd.epreg$td),]
pd.epreg.td <- pd.epreg.td[pd.epreg.td$td=='ASD' | pd.epreg.td$td=='TD',]
pd.epreg.td$td <- as.character(pd.epreg.td$td)
pd.epreg.td$td <- factor(pd.epreg.td$td,levels=c('TD','ASD'))
table(pd.epreg.td$td)
# TD ASD
# 55  28

pdsv.epreg.td <- sva.fit(beta=beta.epreg[,rownames(pd.epreg.td)],pd=pd.epreg.td,var='td')
saveRDS(pdsv.epreg.td,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/pdsv.epreg.td.rds")

pdsv.epreg.td <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/pdsv.epreg.td.rds")
bring.the.heatmap.sv(pdsv.epreg.td,covar=c('Round','predictedSex','PC1','CD8T','CD4T','NK','Bcell','Mono','Gran','Maternal.Age','asd5','td'),folder=FOLDER)

lambdas.epreg <- sva.sequential(beta=beta.epreg[,rownames(pdsv.epreg.td)],pd=pdsv.epreg.td,var='td',folder='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/asd-v-td-sv/')
save(lambdas.epreg,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/lambdas.epreg.rda")
rm(pd.epreg,pd.epreg.td,pdsv.epreg.td,lambdas.epreg,beta.epreg)


#ntd vs td
FOLDER <- '/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/ntd-v-td-sv/'
pd.epreg.ntd <- pd.epreg[!is.na(pd.epreg$td),]
pd.epreg.ntd <- pd.epreg.ntd[pd.epreg.ntd$td=='Non-TD' | pd.epreg.ntd$td=='TD',]
pd.epreg.ntd$td <- as.character(pd.epreg.ntd$td)
pd.epreg.ntd$td <- factor(pd.epreg.ntd$td,levels=c('TD','Non-TD'))
table(pd.epreg.ntd$td)
# TD Non-TD
# 55     60

pdsv.epreg.ntd <- sva.fit(beta=beta.epreg[,rownames(pd.epreg.ntd)],pd=pd.epreg.ntd,var='td')
saveRDS(pdsv.epreg.ntd,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/pdsv.epreg.ntd.rds")

pdsv.epreg.td <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/pdsv.epreg.ntd.rds")
bring.the.heatmap.sv(pdsv.epreg.ntd,covar=c('Round','predictedSex','PC1','CD8T','CD4T','NK','Bcell','Mono','Gran','Maternal.Age','td'),folder=FOLDER)

lambdas.epreg <- sva.sequential(beta=beta.epreg[,rownames(pdsv.epreg.ntd)],pd=pdsv.epreg.ntd,var='td',folder=FOLDER)
save(lambdas.epreg,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/lambdas.epreg.ntd.rda")


### late pregnancy
pd.lpreg <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/pd.lpreg.rds")
beta.lpreg <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/noob.lpreg.rds")

#asd vs td
FOLDER <- '/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/asd-v-td-sv/'
pd.lpreg.td <- pd.lpreg[!is.na(pd.lpreg$td),]
pd.lpreg.td <- pd.lpreg.td[pd.lpreg.td$td=='ASD' | pd.lpreg.td$td=='TD',]
pd.lpreg.td$td <- as.character(pd.lpreg.td$td)
pd.lpreg.td$td <- factor(pd.lpreg.td$td,levels=c('TD','ASD'))
table(pd.lpreg.td$td)
# TD ASD
# 27  15

pdsv.lpreg.td <- sva.fit(beta=beta.lpreg[,rownames(pd.lpreg.td)],pd=pd.lpreg.td,var='td')
saveRDS(pdsv.lpreg.td,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/pdsv.lpreg.td.rds")

pdsv.lpreg.td <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/pdsv.lpreg.td.rds")
bring.the.heatmap.sv(pdsv.lpreg.td,covar=c('Round','predictedSex','PC1','CD8T','CD4T','NK','Bcell','Mono','Gran','Maternal.Age','asd5','td'),folder=FOLDER)

lambdas.lpreg <- sva.sequential(beta=beta.lpreg[,rownames(pdsv.lpreg.td)],pd=pdsv.lpreg.td,var='td',folder=FOLDER)
save(lambdas.lpreg,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/lambdas.lpreg.rda")
rm(pd.lpreg.td,pdsv.lpreg.td,lambdas.lpreg)

#ntd vs td
FOLDER <- '/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/ntd-v-td-sv/'
pd.lpreg.ntd <- pd.lpreg[!is.na(pd.lpreg$td),]
pd.lpreg.ntd <- pd.lpreg.ntd[pd.lpreg.ntd$td=='Non-TD' | pd.lpreg.ntd$td=='TD',]
pd.lpreg.ntd$td <- as.character(pd.lpreg.ntd$td)
pd.lpreg.ntd$td <- factor(pd.lpreg.ntd$td,levels=c('TD','Non-TD'))
table(pd.lpreg.ntd$td)
# TD Non-TD
# 27     36

pdsv.lpreg.ntd <- sva.fit(beta=beta.lpreg[,rownames(pd.lpreg.ntd)],pd=pd.lpreg.ntd,var='td')
saveRDS(pdsv.lpreg.ntd,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/pdsv.lpreg.ntd.rds")

pdsv.lpreg.td <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/pdsv.lpreg.ntd.rds")
bring.the.heatmap.sv(pdsv.lpreg.ntd,covar=c('Round','predictedSex','PC1','CD8T','CD4T','NK','Bcell','Mono','Gran','Maternal.Age','td'),folder=FOLDER)

lambdas.lpreg <- sva.sequential(beta=beta.lpreg[,rownames(pdsv.lpreg.ntd)],pd=pdsv.lpreg.ntd,var='td',folder=FOLDER)
save(lambdas.lpreg,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/lambdas.lpreg.ntd.rda")
rm(pd.lpreg,pd.lpreg.ntd,pdsv.lpreg.ntd,lambdas.lpreg,beta.lpreg)


### placenta fetal
pd.plc <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/pd.plc.rds")
beta.plc <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/noob.plc.rds")

#asd vs td
FOLDER <- '/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/asd-v-td-sv/'
pd.plc.td <- pd.plc[!is.na(pd.plc$td),]
pd.plc.td <- pd.plc.td[pd.plc.td$td=='ASD' | pd.plc.td$td=='TD',]
pd.plc.td$td <- as.character(pd.plc.td$td)
pd.plc.td$td <- factor(pd.plc.td$td,levels=c('TD','ASD'))
table(pd.plc.td$td)
# TD ASD
# 35  18

pdsv.plc.td <- sva.fit(beta=beta.plc[,rownames(pd.plc.td)],pd=pd.plc.td,var='td')
saveRDS(pdsv.plc.td,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/pdsv.plc.td.rds")

pdsv.plc.td <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/pdsv.plc.td.rds")
bring.the.heatmap.sv(pdsv.plc.td,covar=c('Hyb.date','predictedSex','PC1','Maternal.Age','asd5','td'),folder=FOLDER)

lambdas.plc <- sva.sequential(beta=beta.plc[,rownames(pdsv.plc.td)],pd=pdsv.plc.td,var='td',folder=FOLDER)
save(lambdas.plc, file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/lambdas.plc.rda")
rm(pd.plc.td,pdsv.plc.td,lambdas.plc)

#ntd vs td
FOLDER <- '/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/ntd-v-td-sv/'
pd.plc.ntd <- pd.plc[!is.na(pd.plc$td),]
pd.plc.ntd <- pd.plc.ntd[pd.plc.ntd$td=='Non-TD' | pd.plc.ntd$td=='TD',]
pd.plc.ntd$td <- as.character(pd.plc.ntd$td)
pd.plc.ntd$td <- factor(pd.plc.ntd$td,levels=c('TD','Non-TD'))
table(pd.plc.ntd$td)
# TD Non-TD
# 35     48

pdsv.plc.ntd <- sva.fit(beta=beta.plc[,rownames(pd.plc.ntd)],pd=pd.plc.ntd,var='td')
saveRDS(pdsv.plc.ntd,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/pdsv.plc.ntd.rds")

pdsv.plc.td <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/pdsv.plc.ntd.rds")
bring.the.heatmap.sv(pdsv.plc.ntd,covar=c('Hyb.date','predictedSex','PC1','Maternal.Age','td'),folder=FOLDER)

lambdas.plc <- sva.sequential(beta=beta.plc[,rownames(pdsv.plc.ntd)],pd=pdsv.plc.ntd,var='td',folder=FOLDER)
save(lambdas.plc,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/lambdas.plc.ntd.rda")
rm(pd.plc,pd.plc.ntd,pdsv.plc.ntd,lambdas.plc,beta.plc)


### placenta maternal
pd.plm <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/pd.plm.rds")
beta.plm <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/noob.plm.rds")

#asd vs td
FOLDER <- '/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/asd-v-td-sv/'
pd.plm.td <- pd.plm[!is.na(pd.plm$td),]
pd.plm.td <- pd.plm.td[pd.plm.td$td=='ASD' | pd.plm.td$td=='TD',]
pd.plm.td$td <- as.character(pd.plm.td$td)
pd.plm.td$td <- factor(pd.plm.td$td,levels=c('TD','ASD'))
table(pd.plm.td$td)
# TD ASD
# 38  15

pdsv.plm.td <- sva.fit(beta=beta.plm[,rownames(pd.plm.td)],pd=pd.plm.td,var='td')
saveRDS(pdsv.plm.td,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/pdsv.plm.td.rds")

pdsv.plm.td <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/pdsv.plm.td.rds")
bring.the.heatmap.sv(pdsv.plm.td,covar=c('Hyb.date','predictedSex','PC1','Maternal.Age','asd5','td'),folder=FOLDER)

lambdas.plm <- sva.sequential(beta=beta.plm[,rownames(pdsv.plm.td)],pd=pdsv.plm.td,var='td',folder=FOLDER)
save(lambdas.plm, file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/lambdas.plm.rda")
rm(pd.plm,pd.plm.td,pdsv.plm.td,lambdas.plm,beta.plm)

#ntd vs td
FOLDER <- '/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/ntd-v-td-sv/'
pd.plm.ntd <- pd.plm[!is.na(pd.plm$td),]
pd.plm.ntd <- pd.plm.ntd[pd.plm.ntd$td=='Non-TD' | pd.plm.ntd$td=='TD',]
pd.plm.ntd$td <- as.character(pd.plm.ntd$td)
pd.plm.ntd$td <- factor(pd.plm.ntd$td,levels=c('TD','Non-TD'))
table(pd.plm.ntd$td)
# TD Non-TD
# 38     48

pdsv.plm.ntd <- sva.fit(beta=beta.plm[,rownames(pd.plm.ntd)],pd=pd.plm.ntd,var='td')
saveRDS(pdsv.plm.ntd,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/pdsv.plm.ntd.rds")

pdsv.plm.td <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/pdsv.plm.ntd.rds")
bring.the.heatmap.sv(pdsv.plm.ntd,covar=c('Hyb.date','predictedSex','PC1','Maternal.Age','td'),folder=FOLDER)

lambdas.plm <- sva.sequential(beta=beta.plm[,rownames(pdsv.plm.ntd)],pd=pdsv.plm.ntd,var='td',folder=FOLDER)
save(lambdas.plm,file="/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/lambdas.plm.ntd.rda")
rm(pd.plm,pd.plm.ntd,pdsv.plm.ntd,lambdas.plm,beta.plm)



#############################################
#Models
#############################################
library(limma)
library(qqman)
source("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/ASD_sens_functions.R")

#############################################
#cord blood
setwd('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/')

#asd vs td
pd.cord <- readRDS('pdsv.cord.td.rds')
beta.cord <- readRDS('noob.cord.rds')

design <- model.matrix(~pd.cord$td + pd.cord$sv1 + pd.cord$sv2 + pd.cord$sv3)
fit <- lmFit(beta.cord[,rownames(pd.cord)],design)
fit <- eBayes(fit)
res.sv.cord <- topTable(fit,coef=2,nrow(beta.cord))
save(res.sv.cord,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/asd-v-td-sv/res.sv.cord.rda')
rm(pd.cord,beta.cord,res.sv.cord)

#ntd vs td
pd.cord <- readRDS('pdsv.cord.ntd.rds')
beta.cord <- readRDS('noob.cord.rds')

design <- model.matrix(~pd.cord$td + pd.cord$sv1 + pd.cord$sv2)
fit <- lmFit(beta.cord[,rownames(pd.cord)],design)
fit <- eBayes(fit)
res.sv.cord <- topTable(fit,coef=2,nrow(beta.cord))
save(res.sv.cord,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/ntd-v-td-sv/res.sv.cord.rda')
rm(pd.cord,beta.cord,res.sv.cord)


#############################################
#ASD diagnosis, early pregnancy
setwd('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/')

#asd vs td
pd.epreg <- readRDS('pdsv.epreg.td.rds')
beta.epreg <- readRDS('noob.epreg.rds')

design <- model.matrix(~pd.epreg$td + pd.epreg$sv1 + pd.epreg$sv2 + pd.epreg$sv3)
fit <- lmFit(beta.epreg[,rownames(pd.epreg)],design)
fit <- eBayes(fit)
res.sv.epreg <- topTable(fit,coef=2,nrow(beta.epreg))
save(res.sv.epreg,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/asd-v-td-sv/res.sv.epreg.rda')
rm(pd.epreg,beta.epreg,res.sv.epreg)

#ntd vs td
pd.epreg <- readRDS('pdsv.epreg.ntd.rds')
beta.epreg <- readRDS('noob.epreg.rds')

design <- model.matrix(~pd.epreg$td + pd.epreg$sv1 + pd.epreg$sv2 + pd.epreg$sv3 + pd.epreg$sv4 + pd.epreg$sv5 + pd.epreg$sv6)
fit <- lmFit(beta.epreg[,rownames(pd.epreg)],design)
fit <- eBayes(fit)
res.sv.epreg <- topTable(fit,coef=2,nrow(beta.epreg))
save(res.sv.epreg,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/ntd-v-td-sv/res.sv.epreg.rda')
rm(pd.epreg,beta.epreg,res.sv.epreg)

#############################################
#ASD diagnosis, late pregnancy
setwd('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/')

#asd vs td
pd.lpreg <- readRDS('pdsv.lpreg.td.rds')
beta.lpreg <- readRDS('noob.lpreg.rds')

design <- model.matrix(~pd.lpreg$td + pd.lpreg$sv1 + pd.lpreg$sv2 + pd.lpreg$sv3)
fit <- lmFit(beta.lpreg[,rownames(pd.lpreg)],design)
fit <- eBayes(fit)
res.sv.lpreg <- topTable(fit,coef=2,nrow(beta.lpreg))
save(res.sv.lpreg,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/asd-v-td-sv/res.sv.lpreg.rda')
rm(pd.lpreg,beta.lpreg,res.sv.lpreg)

#ntd vs td
pd.lpreg <- readRDS('pdsv.lpreg.ntd.rds')
beta.lpreg <- readRDS('noob.lpreg.rds')

design <- model.matrix(~pd.lpreg$td + pd.lpreg$sv1 + pd.lpreg$sv2 + pd.lpreg$sv3 + pd.lpreg$sv4 + pd.lpreg$sv5)
fit <- lmFit(beta.lpreg[,rownames(pd.lpreg)],design)
fit <- eBayes(fit)
res.sv.lpreg <- topTable(fit,coef=2,nrow(beta.lpreg))
save(res.sv.lpreg,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/ntd-v-td-sv/res.sv.lpreg.rda')
rm(pd.lpreg,beta.lpreg,res.sv.lpreg)

###############################################
#placenta child asd
setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/")

#asd vs td
pd.plc <- readRDS('pdsv.plc.td.rds')
beta.plc <- readRDS('noob.plc.rds')

design <- model.matrix(~pd.plc$td + pd.plc$sv1 + pd.plc$sv2 + pd.plc$sv3)
fit <- lmFit(beta.plc[,rownames(pd.plc)],design)
fit <- eBayes(fit)
res.sv.plc <- topTable(fit,coef=2,nrow(beta.plc))
save(res.sv.plc,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/asd-v-td-sv/res.sv.plc.rda')
rm(pd.plc,beta.plc,res.sv.plc)

#ntd vs td
pd.plc <- readRDS('pdsv.plc.ntd.rds')
beta.plc <- readRDS('noob.plc.rds')

design <- model.matrix(~pd.plc$td + pd.plc$sv1 + pd.plc$sv2 + pd.plc$sv3 + pd.plc$sv4 + pd.plc$sv5 + pd.plc$sv6 + pd.plc$sv7 + pd.plc$sv8 + pd.plc$sv9 + pd.plc$sv10 + pd.plc$sv11 + pd.plc$sv12)
fit <- lmFit(beta.plc[,rownames(pd.plc)],design)
fit <- eBayes(fit)
res.sv.plc <- topTable(fit,coef=2,nrow(beta.plc))
save(res.sv.plc,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/ntd-v-td-sv/res.sv.plc.rda')
rm(pd.plc,beta.plc,res.sv.plc)

###############################################
#placenta mom asd
setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/")

#asd vs td
pd.plm <- readRDS('pdsv.plm.td.rds')
beta.plm <- readRDS('noob.plm.rds')

design <- model.matrix(~pd.plm$td + pd.plm$sv1 + pd.plm$sv2)
fit <- lmFit(beta.plm[,rownames(pd.plm)],design)
fit <- eBayes(fit)
res.sv.plm <- topTable(fit,coef=2,nrow(beta.plm))
save(res.sv.plm,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/asd-v-td-sv/res.sv.plm.rda')
rm(pd.plm,beta.plm,res.sv.plm)

#ntd vs td
pd.plm <- readRDS('pdsv.plm.ntd.rds')
beta.plm <- readRDS('noob.plm.rds')

design <- model.matrix(~pd.plm$td + pd.plm$sv1 + pd.plm$sv2 + pd.plm$sv3 + pd.plm$sv4 + pd.plm$sv5)
fit <- lmFit(beta.plm[,rownames(pd.plm)],design)
fit <- eBayes(fit)
res.sv.plm <- topTable(fit,coef=2,nrow(beta.plm))
save(res.sv.plm,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/ntd-v-td-sv/res.sv.plm.rda')
rm(pd.plm,beta.plm,res.sv.plm)


###############
# qq plots

load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/asd-v-td-sv/res.sv.cord.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/asd-v-td-sv/res.sv.epreg.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/asd-v-td-sv/res.sv.lpreg.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/asd-v-td-sv/res.sv.plc.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/asd-v-td-sv/res.sv.plm.rda")

setwd('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/')

plot.qq <- function(ss.hits,lab,file){
  png(file, type='cairo')
  qq(ss.hits$P.Value,main=paste0(lab," Lambda = ",round(qchisq(median(ss.hits$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1),2)))
  dev.off()
  
}

plot.qq(res.sv.cord,'Cord','/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/qq-sv-cord.png')
plot.qq(res.sv.epreg,'Mom Blood Early Preg','/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/qq-sv-epreg.png')
plot.qq(res.sv.lpreg,'Mom Blood Late Preg','/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/qq-sv-lpreg.png')
plot.qq(res.sv.plc,'Placenta Fetal','/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/qq-sv-plc.png')
plot.qq(res.sv.plm,'Placenta Maternal','/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/qq-sv-plm.png')




#############################################
#no sv model fits
#############################################

#cord blood
setwd('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/')
load('noob_beta_cord-ASD.rda')
load('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/pd.cord.td.rda')
rownames(pd.cord.asd) <- gsub('.*/','',pd.cord.asd$Basename)

library(sas7bdat)
covars <- read.sas7bdat("/dcl01/NDEpi/data/Projects/InProgress/jdou/MARBLES/e_m_covars_v4.sas7bdat")
table(pd.cord.asd$Subject.ID %in% covars$coi_id)
cov <- covars[match(pd.cord.asd$Subject.ID, covars$coi_id),]

### mat edu
pd.cord.asd$mom_edu <- cov$MomEdu
pd.cord.asd$mom_edu_bn <- ifelse(pd.cord.asd$mom_edu<4, "No Degree", "College Degree")

### maternal age
pd.cord.asd$mom_age <- cov$MomAgeYr

### gest age
pd.cord.asd$gest_age <- cov$GA_deliv_wks

pd.cord.asd <- pd.cord.asd[!is.na(pd.cord.asd$mom_edu_bn),]
design <- model.matrix(~ pd.cord.asd$asd5 + pd.cord.asd$predictedSex + pd.cord.asd$Round + pd.cord.asd$mom_age + pd.cord.asd$gest_age + pd.cord.asd$mom_edu_bn + pd.cord.asd$PC1 + pd.cord.asd$PC2 + pd.cord.asd$Gran + pd.cord.asd$nRBC)

fit <- lmFit(noob.beta.cordasd[,rownames(pd.cord.asd)],design)
fit <- eBayes(fit)
res.nsv.cord <- topTable(fit,coef=2,nrow(noob.beta.cordasd))
save(res.nsv.cord, file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/asd-v-td-sv/res.nsv.cord.rda')


#ASD diagnosis, early pregnancy
setwd('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/')

### old model fit
load('noob_beta_momearly-ASD-20170811.rda')
load('pd_epreg-ASD.rda')
sv.vars0 <- as.formula("~ pd.earlypreg.asd$ASD_DSM5_Sib_36mos + factor(pd.earlypreg.asd$Round) + pd.earlypreg.asd$PC1 + pd.earlypreg.asd$Gran + pd.earlypreg.asd$Maternal.Age")

design <- model.matrix(sv.vars0)
fit <- lmFit(noob.beta.epreg,design)
fit <- eBayes(fit)
res.nsv.epreg <- topTable(fit,coef=2,nrow(noob.beta.epreg))
save(res.nsv.epreg, file='Results_MomBloodEarlyPreg-ASD_NoSV-20170824.rda')

#######global methylation
library(metafor)

#cpg annotation
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
an450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
island_enh <- an450k[rownames(res.nsv.epreg),c('Relation_to_Island', 'Enhancer')]
res.nsv.epreg$island <- island_enh[,1]
res.nsv.epreg$enhancer <- island_enh[,2]
save(res.nsv.epreg,file='Results_MomBloodEarlyPreg-ASD_NoSV-20170824.rda')
rm(island_enh)

load('Results_MomBloodEarlyPreg-ASD_4SV.rda')
global <- island.relation(res.sv.epreg)

pdf("GlobalForestPlot.pdf")
par(mar=c(4,4,0,2))
forest(global[,2], ci.lb=global$lower, ci.ub=global$upper, 
       annotate=FALSE, slab=global$feature ,
       xlab="Difference in mean percent methylation")
par(font=2)
text(0, 10, "Maternal Blood Early Pregnancy, Mothers of ASD5 Cases vs Controls")
dev.off()

#############################################
#ASD diagnosis, late pregnancy
setwd('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/')
load('noob_beta_momlate-ASD-20170811.rda')
load('pd_lpreg-ASD.rda')


sv.vars0 <- as.formula("~ pd.latepreg.asd$ASD_DSM5_Sib_36mos + factor(pd.latepreg.asd$Round) + pd.latepreg.asd$PC1 + pd.latepreg.asd$Gran + pd.latepreg.asd$Maternal.Age")

design <- model.matrix(sv.vars0)
fit <- lmFit(noob.beta.lpreg,design)
fit <- eBayes(fit)
res.nsv.lpreg <- topTable(fit,coef=2,nrow(noob.beta.lpreg))
lambda <- qchisq(median(res.nsv.lpreg$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
save(res.nsv.lpreg, file='Results_MomBloodLatePreg-ASD_NoSV-20170824.rda')

#######global methylation
library(metafor)

#cpg annotation
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
an450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
island_enh <- an450k[rownames(res.nsv.lpreg),c('Relation_to_Island', 'Enhancer')]
res.nsv.lpreg$island <- island_enh[,1]
res.nsv.lpreg$enhancer <- island_enh[,2]
save(res.nsv.lpreg,file='Results_MomBloodLatePreg-ASD_NoSV-20170824.rda')
rm(island_enh)

load('Results_MomBloodLatePreg-ASD_3SV.rda')
global <- island.relation(res.sv.lpreg)

pdf("GlobalForestPlot.pdf")
par(mar=c(4,4,0,2))
forest(global[,2], ci.lb=global$lower, ci.ub=global$upper,
       annotate=FALSE, slab=global$feature, 
       xlab="Difference in mean percent methylation")
par(font=2)
text(0.07, 10, "Maternal Blood Late Pregnancy, Mothers of ASD5 Cases vs Controls")
dev.off()

###############################################
### placenta child asd
setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/")
load("pd_placentachild-ASD.rda")
load("noob_beta_placentachild.rda")

names(pd.plc)[27] <- "gest.age"
pd.plc <- pd.plc[pd.plc$gest.age > 1,]
pd.plc <- pd.plc[!is.na(pd.plc$PC1),]
noob.beta.plc <- noob.beta.plc[,rownames(pd.plc)]

vars <- as.formula("~ pd.plc$ASD_DSM5_Sib_36mos + pd.plc$PC1 + factor(pd.plc$predictedSex) + as.numeric(pd.plc$Maternal.Age)")
vars <- as.formula("~ pd.plc$ASD_DSM5_Sib_36mos + pd.plc$PC1 + pd.plc$PC2 + pd.plc$PC3 + pd.plc$PC4 + pd.plc$PC5 + pd.plc$PC6 + pd.plc$PC7 + factor(pd.plc$predictedSex) + factor(pd.plc$Hyb.date) + as.numeric(pd.plc$Maternal.Age)")
design <- model.matrix(vars)
fit <- lmFit(noob.beta.plc,design)
fit <- eBayes(fit)
#saveRDS(fit,file='plc.fit.rda')
res.nsv.plc <- topTable(fit,coef=2,nrow(noob.beta.plc))
lambda <- qchisq(median(res.nsv.plc$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
#0.5914
png('QQ-PlacentaChild-asd+pc1+sex+maternalage+gestage.png')
qq(res.nsv.plc$P.Value,  xlab='Expected -log10(p)', ylab='Observed -log10(p)')
dev.off()
save(res.nsv.plc, file='Results_PlacentaChild-ASD_NoSV-20170825.rda')

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
an450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
island_enh <- an450k[rownames(res.nsv.plc),c('Relation_to_Island', 'Enhancer')]
res.nsv.plc$island <- island_enh[,1]
res.nsv.plc$enhancer <- island_enh[,2]
save(res.nsv.plc,file='Results_PlacentaChild-ASD_NoSV-20170825.rda')

###############################################
### placenta mom asd
setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/")
load("pd_placentamom-ASD.rda")
load("noob_beta_placentamom.rda")

names(pd.plm)[27] <- "gest.age"
pd.plm <- pd.plm[pd.plm$gest.age > 1,]
pd.plm <- pd.plm[!is.na(pd.plm$PC1),]
noob.beta.plm <- noob.beta.plm[,rownames(pd.plm)]

vars <- as.formula("~ pd.plm$ASD_DSM5_Sib_36mos + pd.plm$PC1 + factor(pd.plm$predictedSex) + as.numeric(pd.plm$Maternal.Age) + pd.plm$gest.age")
vars <- as.formula("~ pd.plc$ASD_DSM5_Sib_36mos + pd.plc$PC1 + as.numeric(pd.plc$Maternal.Age) + pd.plc$gest.age")
design <- model.matrix(vars)
fit <- lmFit(noob.beta.plm,design)
fit <- eBayes(fit)
saveRDS(fit,file='plm.fit.rda')
res.nsv.plm <- topTable(fit,coef=2,nrow(noob.beta.plm))
lambda <- qchisq(median(res.nsv.plm$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
#1.2178
png('QQ-PlacentaMom-asd+pc1+sex+maternalage+gestage.png')
qq(res.nsv.plm$P.Value, main=paste('lambda = ',round(lambda,2),sep=''))
dev.off()
save(res.nsv.plm, file='Results_PlacentaMom-ASD_NoSV-20170825.rda')

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
an450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
island_enh <- an450k[rownames(res.nsv.plm),c('Relation_to_Island', 'Enhancer')]
res.nsv.plm$island <- island_enh[,1]
res.nsv.plm$enhancer <- island_enh[,2]
save(res.nsv.plm,file='Results_PlacentaMom-ASD_NoSV-20170825.rda')


### qq plots
library(qqman)

setwd('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/')

load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/Results_CordBlood_ASD_NoSV-20171026.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/Results_MomBloodEarlyPreg-ASD_NoSV-20170824.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/Results_MomBloodLatePreg-ASD_NoSV-20170824.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/Results_PlacentaChild-ASD_NoSV-20170825.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/Results_PlacentaMom-ASD_NoSV-20170825.rda")

plot.qq <- function(ss.hits,lab,file){
  png(file, type='cairo')
  qq(ss.hits$P.Value,main=paste0(lab," Lambda = ",round(qchisq(median(ss.hits$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1),2)))
  dev.off()
  
}

plot.qq(res.nsv.asd,'Cord','/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/qq-cord.png')
plot.qq(res.nsv.epreg,'Mom Blood Early Preg','/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/qq-epreg.png')
plot.qq(res.nsv.lpreg,'Mom Blood Late Preg','/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/qq-lpreg.png')
plot.qq(res.nsv.plc,'Placenta Fetal','/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/qq-plc.png')
plot.qq(res.nsv.plm,'Placenta Maternal','/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/qq-plm.png')

###############################################
#svs in cord asd model?
###############################################
mod <- read.csv("/dcl01/NDEpi/data/Projects/InProgress/taung/dnam_asd/dnam_asd_proj/modSv.csv"); mod <- mod[,3:16]
library(limma)
mod1 <- mod[,-c(9:14)]
fit <- lmFit(beta.n,mod)
fit <- eBayes(fit)
ss.hits <- topTable(fit,coef=2,nrow(beta.n))
fit <- lmFit(beta.n,mod1)
fit <- eBayes(fit)
ss.hits.nsv <- topTable(fit,coef=2,nrow(beta.n))
island_enh_mod <- an450k[rownames(ss.hits),c('Relation_to_Island', 'Enhancer', 'chr')]
ss.hits$island <- island_enh_mod[,1]
ss.hits$enhancer <- island_enh_mod[,2]
island_enh_mod <- an450k[rownames(ss.hits.nsv),c('Relation_to_Island', 'Enhancer', 'chr')]
ss.hits.nsv$island <- island_enh_mod[,1]
ss.hits.nsv$enhancer <- island_enh_mod[,2]
global.mod <- island.relation(ss.hits)
global.mod1 <- island.relation(ss.hits.nsv)

res.nsv.asd <- ss.hits.nsv
save(res.nsv.asd,file='Results_CordBlood_ASD_NoSV-20170824.rda')

ss.hits.nsvo <- ss.hits.nsv[rownames(ss.hits),]
aver = (ss.hits[ss.hits$island=='S_Shelf',]$logFC + ss.hits.nsvo[ss.hits.nsvo$island=='S_Shelf',]$logFC)/2
diff = ss.hits[ss.hits$island=='S_Shelf',]$logFC - ss.hits.nsvo[ss.hits.nsvo$island=='S_Shelf',]$logFC
smoothScatter(aver,diff)
smoothScatter(ss.hits$P.Value,ss.hits.nsvo$P.Value)

coefficients <- read.csv("/dcl01/NDEpi/data/Projects/InProgress/taung/dnam_asd/dnam_asd_proj/beta.adj.csv")
rownames(coefficients) <- coefficients$X
coefficients <- coefficients[,-1]
design <- read.csv("/dcl01/NDEpi/data/Projects/InProgress/taung/dnam_asd/dnam_asd_proj/modSv.csv")
design <- design[,3:16]
coefficients <- coefficients[,c(1,2,5,3,4,6,7,8,9,10,11,12,13,14)]
coefficients <- as.matrix(coefficients)
design <- as.matrix(design)
design <- t(design)
coef <- coefficients[,-2]
des <- design[-2,]
adj.dnaM <- coef%*%des 
design <- t(design)
beta.asd1 <- adj.dnaM[,which(design[,2]==1)]
beta.asd0 <- adj.dnaM[,which(design[,2]==0)]
beta.asd1 <- merge(beta.asd1,island_enh_cord,by='row.names',sort=F,all.x=T,all.y=F)
beta.asd0 <- merge(beta.asd0,island_enh_cord,by='row.names',sort=F,all.x=T,all.y=F)
beta.asd1$Rowmeans <- rowMeans(beta.asd1[,2:32])
beta.asd0$Rowmeans <- rowMeans(beta.asd0[,2:111])
beta.asd1 <- beta.asd1[beta.asd1$chr!='chrX'&beta.asd1$chr!='chrY',]
beta.asd0 <- beta.asd0[beta.asd0$chr!='chrX'&beta.asd0$chr!='chrY',]
global.matrixmethod <- island.relation.matrix(beta.asd1,beta.asd0)
global.matrixmethod

setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/")
nosv.vars <- as.formula("~ pd.cord.asd$asd5 +pd.cord.asd$PC1 + pd.cord.asd$predictedSex + factor(pd.cord.asd$Round) + pd.cord.asd$Maternal.Age + pd.cord.asd$Gran + pd.cord.asd$nRBC")
design <- model.matrix(nosv.vars)
fit <- lmFit(noob.beta.cordasd,design)
fit <- eBayes(fit)
res.nsv.asd <- topTable(fit,coef=2,nrow(noob.beta.cordasd))
save(res.nsv.asd,file='Results_CordBlood_ASD_NoSV-20171026.rda')

#svs in mom T3?
setwd('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/')
load('noob_beta_momlate-ASD.rda')
load('pd_lpreg-ASD.rda')
nosv.vars <- as.formula("~ pd.latepreg.asd$ASD_DSM5_Sib_36mos + factor(pd.latepreg.asd$Round) + pd.latepreg.asd$PC1 + pd.latepreg.asd$Gran + pd.latepreg.asd$Maternal.Age")
design <- model.matrix(nosv.vars)
fit <- lmFit(noob.beta.lpreg,design)
fit <- eBayes(fit)
res.nsv <- topTable(fit,coef=2,nrow(noob.beta.lpreg))
island_enh <- an450k[rownames(res.nsv),c('Relation_to_Island', 'Enhancer', 'chr')]
res.nsv$island <- island_enh[,1]
res.nsv$enhancer <- island_enh[,2]
global.momT3.nosv <- island.relation(res.nsv)
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/Results_MomBloodLatePreg-ASD_3SV.rda")
island_enh <- an450k[rownames(res.sv.lpreg),c('Relation_to_Island', 'Enhancer', 'chr')]
res.sv.lpreg$island <- island_enh[,1]
res.sv.lpreg$enhancer <- island_enh[,2]
global.momT3.sv <- island.relation(res.sv.lpreg)

res.sv.lpreg <- res.sv.lpreg[rownames(res.nsv),]
aver <- (res.nsv$logFC + res.sv.lpreg$logFC) / 2
diff <- (res.nsv$logFC - res.sv.lpreg$logFC)
smoothScatter(aver,diff,xlab='Average of Two Betas',ylab='No SV Beta - With SV Beta Coefficient')

#svs in mom T1?
setwd('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/')
load('noob_beta_momearly-ASD.rda')
load('pd_epreg-ASD.rda')
nosv.vars <- as.formula("~ pd.earlypreg.asd$ASD_DSM5_Sib_36mos + factor(pd.earlypreg.asd$Round) + pd.earlypreg.asd$PC1 + pd.earlypreg.asd$Gran + pd.earlypreg.asd$Maternal.Age")
design <- model.matrix(nosv.vars)
fit <- lmFit(noob.beta.epreg,design)
fit <- eBayes(fit)
res.nsv <- topTable(fit,coef=2,nrow(noob.beta.epreg))
island_enh <- an450k[rownames(res.nsv),c('Relation_to_Island', 'Enhancer', 'chr')]
res.nsv$island <- island_enh[,1]
res.nsv$enhancer <- island_enh[,2]
global.momT1.nosv <- island.relation(res.nsv)
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/Results_MomBloodEarlyPreg-ASD_4SV-20170811.rda")
island_enh <- an450k[rownames(res.sv.epreg),c('Relation_to_Island', 'Enhancer', 'chr')]
res.sv.epreg$island <- island_enh[,1]
res.sv.epreg$enhancer <- island_enh[,2]
global.momT1.sv <- island.relation(res.sv.lpreg)

#check cases and controls, unadjusted, stratified by sex
beta.asd1m <- beta.n[,rownames(pd[pd$asd5==1 & pd$babygender==1,])]
beta.asd1f <- beta.n[,rownames(pd[pd$asd5==1 & pd$babygender==0,])]
beta.asd0m <- beta.n[,rownames(pd[pd$asd5==0 & pd$babygender==1,])]
beta.asd0f <- beta.n[,rownames(pd[pd$asd5==0 & pd$babygender==0,])]
beta.asd1m <- merge(beta.asd1m,island_enh_cord,by='row.names',sort=F,all.x=T,all.y=F)
beta.asd1f <- merge(beta.asd1f,island_enh_cord,by='row.names',sort=F,all.x=T,all.y=F)
beta.asd0m <- merge(beta.asd0m,island_enh_cord,by='row.names',sort=F,all.x=T,all.y=F)
beta.asd0f <- merge(beta.asd0f,island_enh_cord,by='row.names',sort=F,all.x=T,all.y=F)
beta.asd1m$Rowmeans <- rowMeans(beta.asd1m[,2:24])
beta.asd1f$Rowmeans <- rowMeans(beta.asd1f[,2:9])
beta.asd0m$Rowmeans <- rowMeans(beta.asd0m[,2:53])
beta.asd0f$Rowmeans <- rowMeans(beta.asd0f[,2:59])
beta.asd1m <- beta.asd1m[beta.asd1m$chr!='chrX'&beta.asd1m$chr!='chrY',]
beta.asd1f <- beta.asd1f[beta.asd1f$chr!='chrX'&beta.asd1f$chr!='chrY',]
beta.asd0m <- beta.asd0m[beta.asd0m$chr!='chrX'&beta.asd0m$chr!='chrY',]
beta.asd0f <- beta.asd0f[beta.asd0f$chr!='chrX'&beta.asd0f$chr!='chrY',]
# beta.asd1m <- beta.asd1m[,c("Relation_to_Island","Enhancer","Rowmeans")]
# names(beta.asd1m) <- c("island","enhancer","logFC")
# beta.asd1f <- beta.asd1f[,c("Relation_to_Island","Enhancer","Rowmeans")]
# names(beta.asd1f) <- c("island","enhancer","logFC")
# beta.asd0m <- beta.asd0m[,c("Relation_to_Island","Enhancer","Rowmeans")]
# names(beta.asd0m) <- c("island","enhancer","logFC")
# beta.asd0f <- beta.asd0f[,c("Relation_to_Island","Enhancer","Rowmeans")]
# names(beta.asd0f) <- c("island","enhancer","logFC")
global.case <- island.relation.matrix(beta.asd1f,beta.asd1m)
global.con <- island.relation.matrix(beta.asd0f,beta.asd0m)
global <- rbind(global.case,global.con)
pdf("CordASD-GlobalForestPlot-Sex.pdf")
par(mar=c(4,4,0,2))
forest(global[,2], ci.lb=global$lower, ci.ub=global$upper,
       annotate=FALSE, slab=global$feature, 
       rows=c(0:7,10:17), ylim=c(0,22),
       xlab="Difference in mean percent methylation")
par(font=2)
text(0.1, 22, "Unadjusted Cord Blood Methylation Female - Males")
text(0, c(8,18), pos=4, c("Case", "Control"))
dev.off()


######################################################################
#beta vs beta and -logP vs -logP pairwise

setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/")

cord.asd <- read.csv("/dcl01/NDEpi/data/Projects/InProgress/taung/dnam_asd/dnam_asd_proj/4Kelly/ss.hits.042417.csv")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_36MonthSRS/Results_CordBlood-SRS_7SV.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_12MonthAOSI/Results_CordBlood-AOSI_10SV.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/Results_MomBloodEarlyPreg-ASD_4SV.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/Results_MomBloodLatePreg-ASD_3SV.rda")
rownames(cord.asd) <- cord.asd$X
cord.asd <- cord.asd[,-1]
cord.srs <- res.sv.srs
cord.aosi <- res.sv.aosi
momE.asd <- res.sv.epreg
momL.asd <- res.sv.lpreg
rm(res.sv.aosi,res.sv.srs,res.sv.lpreg,res.sv.epreg)

#beta vs beta
beta_cord.srs <- cord.srs[,'logFC']
beta_cord.asd <- cord.asd[rownames(cord.srs),'logFC']
beta_cord.aosi <- cord.aosi[rownames(cord.srs),'logFC']
beta_momE.asd <- momE.asd[rownames(cord.srs),'logFC']
beta_momL.asd <- momL.asd[rownames(cord.srs),'logFC']
betas <- cbind(beta_cord.asd,beta_cord.srs,beta_cord.aosi,beta_momE.asd,beta_momL.asd)
rownames(betas) <- rownames(cord.srs)

pdf("Beta vs Beta Plots.pdf")
upper.plots <- function(...){
  smoothScatter(...,colramp=colorRampPalette(c('white','slateblue3','slateblue4')), bandwidth=0.003,add=T)
}

lower.plots <- function(...){
  op <- par("usr")
  par(usr=c(0,1,0,1))
  text(0.5,0.5,paste("r=",round(cor(...,use="pairwise.complete.obs"),2),sep=""),cex=1.8)
  par(usr=op)
}
pairs(betas,lower.panel=lower.plots,upper.panel=upper.plots)
dev.off()

#logP vs logP
p_cord.srs <- -log(cord.srs[,'P.Value'],10)
p_cord.asd <- -log(cord.asd[rownames(cord.srs),'P.Value'],10)
p_cord.aosi <- -log(cord.aosi[rownames(cord.srs),'P.Value'],10)
p_momE.asd <- -log(momE.asd[rownames(cord.srs),'P.Value'],10)
p_momL.asd <- -log(momL.asd[rownames(cord.srs),'P.Value'],10)
logp <- cbind(p_cord.asd,p_cord.srs,p_cord.aosi,p_momE.asd,p_momL.asd)
rownames(logp) <- rownames(cord.srs)

pdf("logP vs logP Plots.pdf")
upper.plots <- function(...){
  smoothScatter(...,colramp=colorRampPalette(c('white','slateblue3','slateblue4')), bandwidth=0.003,add=T)
}

lower.plots <- function(...){
  op <- par("usr")
  par(usr=c(0,1,0,1))
  text(0.5,0.5,paste("r=",round(cor(...,use="pairwise.complete.obs"),2),sep=""),cex=1.8)
  par(usr=op)
}
pairs(logp,lower.panel=lower.plots,upper.panel=upper.plots)
dev.off()

#check with no sv regressions to see how different
setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/")

cord.asd <- read.csv("/dcl01/NDEpi/data/Projects/InProgress/taung/dnam_asd/dnam_asd_proj/4Kelly/ss.hits.042417.csv")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_36MonthSRS/Results_CordBlood-SRS_NoSV.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_12MonthAOSI/Results_CordBlood-AOSI_NoSV.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/Results_MomBloodEarlyPreg-ASD_NoSV.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/Results_MomBloodLatePreg-ASD_NoSV.rda")
rownames(cord.asd) <- cord.asd$X
cord.asd <- cord.asd[,-1]
cord.srs <- res.nsv.srs
cord.aosi <- res.sv0
momE.asd <- res.nsv.epreg
momL.asd <- res.nsv.lpreg
rm(res.sv0,res.nsv.srs,res.nsv.lpreg,res.nsv.epreg)

#beta vs beta
beta_cord.srs <- cord.srs[,'logFC']
beta_cord.asd <- cord.asd[rownames(cord.srs),'logFC']
beta_cord.aosi <- cord.aosi[rownames(cord.srs),'logFC']
beta_momE.asd <- momE.asd[rownames(cord.srs),'logFC']
beta_momL.asd <- momL.asd[rownames(cord.srs),'logFC']
betas <- cbind(beta_cord.asd,beta_cord.srs,beta_cord.aosi,beta_momE.asd,beta_momL.asd)
rownames(betas) <- rownames(cord.srs)

pairs(betas, panel = function(...){smoothScatter(...,colramp=colorRampPalette(c('white','slateblue3','slateblue4')), bandwidth=0.003,add=T)})

#logP vs logP
p_cord.srs <- -log(cord.srs[,'P.Value'],10)
p_cord.asd <- -log(cord.asd[rownames(cord.srs),'P.Value'],10)
p_cord.aosi <- -log(cord.aosi[rownames(cord.srs),'P.Value'],10)
p_momE.asd <- -log(momE.asd[rownames(cord.srs),'P.Value'],10)
p_momL.asd <- -log(momL.asd[rownames(cord.srs),'P.Value'],10)
logp <- cbind(p_cord.asd,p_cord.srs,p_cord.aosi,p_momE.asd,p_momL.asd)
rownames(logp) <- rownames(cord.srs)

pairs(logp, panel = function(...){smoothScatter(...,colramp=colorRampPalette(c('white','slateblue3','slateblue4')), bandwidth=0.003,add=T)})


#######################################################################################
#Compare sv to known covariate models
#######################################################################################
setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/")
source("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/ASD_sens_functions.R")

load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/Results_CordBlood_ASD_NoSV-20171026.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/Results_MomBloodEarlyPreg-ASD_NoSV-20170824.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/Results_MomBloodLatePreg-ASD_NoSV-20170824.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/Results_PlacentaChild-ASD_NoSV-20170825.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/Results_PlacentaMom-ASD_NoSV-20170825.rda")

load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/asd-v-td-sv/res.sv.cord.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/asd-v-td-sv/res.sv.epreg.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/asd-v-td-sv/res.sv.lpreg.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/asd-v-td-sv/res.sv.plc.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/asd-v-td-sv/res.sv.plm.rda")

betavsbeta(res.sv.cord,res.nsv.asd,'Cord: SV vs Known Covariate Model','/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/cord-betavsbeta.pdf')
betavsbeta(res.sv.epreg,res.nsv.epreg,'Mom Blood Early Preg: SV vs Known Covariate Model','/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/epreg-betavsbeta.pdf')
betavsbeta(res.sv.lpreg,res.nsv.lpreg,'Mom Blood Late Preg: SV vs Known Covariate Model','/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/lpreg-betavsbeta.pdf')
betavsbeta(res.sv.plc,res.nsv.plc,'Placenta Fetal: SV vs Known Covariate Model','/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/plc-betavsbeta.pdf')
betavsbeta(res.sv.plm,res.nsv.plm,'Placenta Mom: SV vs Known Covariate Model','/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/plm-betavsbeta.pdf')

pvp(res.sv.cord,res.nsv.asd,'Cord: SV vs Known Covariate Model','/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/cord-pvp.pdf')
pvp(res.sv.epreg,res.nsv.epreg,'Mom Blood Early Preg: SV vs Known Covariate Model','/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/epreg-pvp.pdf')
pvp(res.sv.lpreg,res.nsv.lpreg,'Mom Blood Late Preg: SV vs Known Covariate Model','/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/lpreg-pvp.pdf')
pvp(res.sv.plc,res.nsv.plc,'Placenta Fetal: SV vs Known Covariate Model','/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/plc-pvp.pdf')
pvp(res.sv.plm,res.nsv.plm,'Placenta Mom: SV vs Known Covariate Model','/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/plm-pvp.pdf')


#######################################################################################
#SFARI enrichment
#######################################################################################

setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/")
source("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/ASD_sens_functions.R")

# with old known covariates model
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/Results_CordBlood_ASD_NoSV-20171026.rda")
# # load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_36MonthSRS/Results_CordBlood-SRS_NoSV-20170824.rda")
# # load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_12MonthAOSI/Results_CordBlood-AOSI_NoSV-20170824.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/Results_MomBloodEarlyPreg-ASD_NoSV-20170824.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/Results_MomBloodLatePreg-ASD_NoSV-20170824.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/Results_PlacentaChild-ASD_NoSV-20170825.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/Results_PlacentaMom-ASD_NoSV-20170825.rda")

# sparse <- T
# enrich.asd <- sfari.enrich(res.nsv.asd, sparse.p=sparse)
# #enrich.srs <- sfari.enrich(res.nsv.srs, sparse.p=sparse, cut=enrich.asd$peak.p.cutoff)
# #enrich.aosi <- sfari.enrich(res.nsv.aosi, sparse.p=sparse, cut=enrich.asd$peak.p.cutoff)
# enrich.epreg <- sfari.enrich(res.nsv.epreg, sparse.p=sparse, cut=enrich.asd$peak.p.cutoff)
# enrich.lpreg <- sfari.enrich(res.nsv.lpreg, sparse.p=sparse, cut=enrich.asd$peak.p.cutoff)
# enrich.plc <- sfari.enrich(res.nsv.plc, sparse.p=sparse, cut=enrich.asd$peak.p.cutoff)
# enrich.plm <- sfari.enrich(res.nsv.plm, sparse.p=sparse, cut=enrich.asd$peak.p.cutoff)


### sv models - pick one set

#asd vs td
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/asd-v-td-sv/res.sv.cord.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/asd-v-td-sv/res.sv.epreg.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/asd-v-td-sv/res.sv.lpreg.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/asd-v-td-sv/res.sv.plc.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/asd-v-td-sv/res.sv.plm.rda")

# #ntd vs td set
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/ntd-v-td-sv/res.sv.cord.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/ntd-v-td-sv/res.sv.epreg.rda")


sparse <- T
enrich.epreg <- sfari.enrich(res.sv.epreg, sparse.p=sparse,cut=0.05)
enrich.asd <- sfari.enrich(res.sv.cord, sparse.p=sparse,cut=enrich.epreg$peak.p.cutoff)
enrich.lpreg <- sfari.enrich(res.sv.lpreg, sparse.p=sparse, cut=enrich.epreg$peak.p.cutoff)
enrich.plc <- sfari.enrich(res.sv.plc, sparse.p=sparse, cut=enrich.epreg$peak.p.cutoff)
enrich.plm <- sfari.enrich(res.sv.plm, sparse.p=sparse, cut=enrich.epreg$peak.p.cutoff)


#combine 5 plots for ASD models
library(RColorBrewer)
col <- brewer.pal(5, 'Set1')
pdf('sfari-enrichment-all-asd-models-sv.pdf')
par(mar=c(4.5,4.5,3.8,2))
plot(enrich.epreg$cutoffs,enrich.epreg$p.values,type='n',xlab='EWAS P-Value Cut-off',ylab=expression('-Log'[10]*'(P-Value)'),main='SFARI Enrichment',las=1)
axis(1, at=enrich.asd$cutoffs, labels=F, tck=0.04, lwd.ticks=1.2)
lines(enrich.asd$cutoffs,enrich.asd$p.values,col=col[1])
lines(enrich.epreg$cutoffs,enrich.epreg$p.values,col=col[2])
lines(enrich.lpreg$cutoffs,enrich.lpreg$p.values,col=col[3])
lines(enrich.plc$cutoffs,enrich.plc$p.values,col=col[4])
lines(enrich.plm$cutoffs,enrich.plm$p.values,col=col[5])
legend('topright',legend=c('Cord Blood','Maternal Blood Early Pregnancy','Maternal Blood Late Pregnancy','Fetal Placenta','Maternal Placenta'),
       col=col,lty=1, text.col=col, title='Tissue', title.col='black')
dev.off()

#overlap table
over.tab <- matrix(0,nrow=9,ncol=5)
rownames(over.tab) <- c('0.01N','0.01Expected','0.01OR','0.05N','0.05Expected','0.05OR','0.1N','0.1Expected','0.1OR')
colnames(over.tab) <- c('CordBlood','MaternalBloodEarlyPreg','MaternalBloodLatePreg','PlacentaFetal','PlacentaMaternal')

enrich.list <- list(enrich.asd,enrich.epreg,enrich.lpreg,enrich.plc,enrich.plm)
for(enr in 1:5){
  over.tab[1,enr] <- enrich.list[[enr]]$tables$'EWAS P 0.01'[4]
  over.tab[2,enr] <- chisq.test(enrich.list[[enr]]$tables$'EWAS P 0.01')$expected[4]
  over.tab[3,enr] <- enrich.list[[enr]]$OR[4][[1]]
  over.tab[4,enr] <- enrich.list[[enr]]$tables$'EWAS P 0.05'[4]
  over.tab[5,enr] <- chisq.test(enrich.list[[enr]]$tables$'EWAS P 0.05')$expected[4]
  over.tab[6,enr] <- enrich.list[[enr]]$OR[8][[1]]
  over.tab[7,enr] <- enrich.list[[enr]]$tables$'EWAS P 0.1'[4]
  over.tab[8,enr] <- chisq.test(enrich.list[[enr]]$tables$'EWAS P 0.1')$expected[4]
  over.tab[9,enr] <- enrich.list[[enr]]$OR[9][[1]]
}

write.csv(over.tab,file='sv-overlap.table.csv')

#which genes appear in enrichment for each group
sfari<-read.csv("/dcl01/NDEpi/data/Projects/InProgress/kbakulsk/SFARI-Gene_genes_export28-08-2017.csv")
sfari.na<-sfari[!(is.na(sfari$gene.symbol)),]
#sfari.gene.hunt <- data.frame(matrix(0,nrow=881,ncol=8))
sfari.gene.hunt <- data.frame(matrix(0,nrow=881,ncol=6))
#colnames(sfari.gene.hunt) <- c("Gene.Name","CordASD","CordSRS","CordAOSI","MomEPreg","MomLPreg","PlacentaC","PlacentaM")
colnames(sfari.gene.hunt) <- c("Gene.Name","CordASD","MomEPreg","MomLPreg","PlacentaC","PlacentaM")
sfari.gene.hunt[,1] <- as.character(sfari.na$gene.symbol)
sfari.gene.hunt$CordASD <- ifelse(sfari.gene.hunt$Gene.Name %in% enrich.asd$genes$symbols,1,0)
#sfari.gene.hunt$CordSRS <- ifelse(sfari.gene.hunt$Gene.Name %in% enrich.srs$genes$symbols,1,0)
#sfari.gene.hunt$CordAOSI <- ifelse(sfari.gene.hunt$Gene.Name %in% enrich.aosi$genes$symbols,1,0)
sfari.gene.hunt$MomEPreg <- ifelse(sfari.gene.hunt$Gene.Name %in% enrich.epreg$genes$symbols,1,0)
sfari.gene.hunt$MomLPreg <- ifelse(sfari.gene.hunt$Gene.Name %in% enrich.lpreg$genes$symbols,1,0)
sfari.gene.hunt$PlacentaC <- ifelse(sfari.gene.hunt$Gene.Name %in% enrich.plc$genes$symbols,1,0)
sfari.gene.hunt$PlacentaM <- ifelse(sfari.gene.hunt$Gene.Name %in% enrich.plm$genes$symbols,1,0)
#sfari.gene.hunt$freq <- sfari.gene.hunt$CordASD+sfari.gene.hunt$CordSRS+sfari.gene.hunt$CordAOSI+sfari.gene.hunt$MomEPreg+sfari.gene.hunt$MomLPreg+sfari.gene.hunt$PlacentaC+sfari.gene.hunt$PlacentaM
sfari.gene.hunt$freq <- sfari.gene.hunt$CordASD+sfari.gene.hunt$MomEPreg+sfari.gene.hunt$MomLPreg+sfari.gene.hunt$PlacentaC+sfari.gene.hunt$PlacentaM

save(sfari.gene.hunt,file='genes_sfari_enrich_in_models.rda')
rm(res.nsv.asd,res.nsv.srs,res.nsv.aosi,res.nsv.epreg,res.nsv.lpreg,res.nsv.plc,res.nsv.plm,enrich.asd,enrich.srs,enrich.aosi,enrich.epreg,enrich.lpreg,enrich.plc,enrich.plm,sparse)

load('genes_sfari_enrich_in_models.rda')
library(VennDiagram)
n1 = sum(sfari.gene.hunt$CordASD)
n2 = sum(sfari.gene.hunt$PlacentaC)
n3 = sum(sfari.gene.hunt$PlacentaM)
n4 = sum(sfari.gene.hunt$MomEPreg)
n5 = sum(sfari.gene.hunt$MomLPreg)
n12 = sum(sfari.gene.hunt$CordASD==1 & sfari.gene.hunt$PlacentaC==1)
n13 = sum(sfari.gene.hunt$CordASD==1 & sfari.gene.hunt$PlacentaM==1)
n14 = sum(sfari.gene.hunt$CordASD==1 & sfari.gene.hunt$MomEPreg==1)
n15 = sum(sfari.gene.hunt$CordASD==1 & sfari.gene.hunt$MomLPreg==1)
n23 = sum(sfari.gene.hunt$PlacentaC==1 & sfari.gene.hunt$PlacentaM==1)
n24 = sum(sfari.gene.hunt$PlacentaC==1 & sfari.gene.hunt$MomEPreg==1)
n25 = sum(sfari.gene.hunt$PlacentaC==1 & sfari.gene.hunt$MomLPreg==1)
n34 = sum(sfari.gene.hunt$PlacentaM==1 & sfari.gene.hunt$MomEPreg==1)
n35 = sum(sfari.gene.hunt$PlacentaM==1 & sfari.gene.hunt$MomLPreg==1)
n45 = sum(sfari.gene.hunt$MomEPreg==1 & sfari.gene.hunt$MomLPreg==1)
n123 = sum(sfari.gene.hunt$CordASD==1 & sfari.gene.hunt$PlacentaC==1 & sfari.gene.hunt$PlacentaM==1)
n124 = sum(sfari.gene.hunt$CordASD==1 & sfari.gene.hunt$PlacentaC==1 & sfari.gene.hunt$MomEPreg==1)
n125 = sum(sfari.gene.hunt$CordASD==1 & sfari.gene.hunt$PlacentaC==1 & sfari.gene.hunt$MomLPreg==1)
n134 = sum(sfari.gene.hunt$CordASD==1 & sfari.gene.hunt$PlacentaM==1 & sfari.gene.hunt$MomEPreg==1)
n135 = sum(sfari.gene.hunt$CordASD==1 & sfari.gene.hunt$PlacentaM==1 & sfari.gene.hunt$MomLPreg==1)
n145 = sum(sfari.gene.hunt$CordASD==1 & sfari.gene.hunt$MomEPreg==1 & sfari.gene.hunt$MomLPreg==1)
n234 = sum(sfari.gene.hunt$PlacentaC==1 & sfari.gene.hunt$PlacentaM==1 & sfari.gene.hunt$MomEPreg==1)
n235 = sum(sfari.gene.hunt$PlacentaC==1 & sfari.gene.hunt$PlacentaM==1 & sfari.gene.hunt$MomLPreg==1)
n245 = sum(sfari.gene.hunt$PlacentaC==1 & sfari.gene.hunt$MomEPreg==1 & sfari.gene.hunt$MomLPreg==1)
n345 = sum(sfari.gene.hunt$PlacentaM==1 & sfari.gene.hunt$MomEPreg==1 & sfari.gene.hunt$MomLPreg==1)
n1234 = sum(sfari.gene.hunt$CordASD==1 & sfari.gene.hunt$PlacentaC==1 & sfari.gene.hunt$PlacentaM==1 & sfari.gene.hunt$MomEPreg==1)
n1235 = sum(sfari.gene.hunt$CordASD==1 & sfari.gene.hunt$PlacentaC==1 & sfari.gene.hunt$PlacentaM==1 & sfari.gene.hunt$MomLPreg==1)
n1245 = sum(sfari.gene.hunt$CordASD==1 & sfari.gene.hunt$PlacentaC==1 & sfari.gene.hunt$MomEPreg==1 & sfari.gene.hunt$MomLPreg==1)
n1345 = sum(sfari.gene.hunt$CordASD==1 & sfari.gene.hunt$PlacentaM==1 & sfari.gene.hunt$MomEPreg==1 & sfari.gene.hunt$MomLPreg==1)
n2345 = sum(sfari.gene.hunt$PlacentaC==1 & sfari.gene.hunt$PlacentaM==1 & sfari.gene.hunt$MomEPreg==1 & sfari.gene.hunt$MomLPreg==1)
n12345 = sum(sfari.gene.hunt$CordASD==1 & sfari.gene.hunt$PlacentaC==1 & sfari.gene.hunt$PlacentaM==1 & sfari.gene.hunt$MomEPreg==1 & sfari.gene.hunt$MomLPreg==1)

library(RColorBrewer)
col <- brewer.pal(5, 'Set1')
pdf('sv-venn-diagram-5asd-models-v2.pdf')
draw.quintuple.venn(sum(sfari.gene.hunt$CordASD),sum(sfari.gene.hunt$PlacentaC),sum(sfari.gene.hunt$PlacentaM),sum(sfari.gene.hunt$MomEPreg),sum(sfari.gene.hunt$MomLPreg),
                    n12,n13,n14,n15,n23,n24,n25,n34,n35,n45,
                    n123,n124,n125,n134,n135,n145,n234,n235,n245,n345,
                    n1234,n1235,n1245,n1345,n2345,n12345,
                    category=c(paste0('Cord\nBlood (n=',n1,')'),paste0('Fetal\nPlacenta (n=',n2,')'),paste0('Maternal\nPlacenta (n=',n3,')'),paste0('Maternal Blood\nEarly Preg (n=',n4,')'),paste0('Maternal Blood\nLate Preg (n=',n5,')')),
                    fill=col,cat.col=col,mar=c(0.15,0.15,0.15,0.15),cat.dist=c(0.3,0.3,0.3,0.3,0.3),cat.fontfamily="sans",fontfamily="sans",cat.cex=0.9)
dev.off()


# overlap for non-sfari too
sfari<-read.csv("/dcl01/NDEpi/data/Projects/InProgress/kbakulsk/SFARI-Gene_genes_export28-08-2017.csv")
nsf.cord <- enrich.asd$genes.at.cut[!enrich.asd$genes.at.cut %in% sfari$gene.symbol]
nsf.epreg <- enrich.epreg$genes.at.cut[!enrich.epreg$genes.at.cut %in% sfari$gene.symbol]
nsf.lpreg <- enrich.lpreg$genes.at.cut[!enrich.lpreg$genes.at.cut %in% sfari$gene.symbol]
nsf.plc <- enrich.plc$genes.at.cut[!enrich.plc$genes.at.cut %in% sfari$gene.symbol]
nsf.plm <- enrich.plm$genes.at.cut[!enrich.plm$genes.at.cut %in% sfari$gene.symbol]

all <- unique(c(nsf.cord,nsf.epreg,nsf.lpreg,nsf.plc,nsf.plm))
counts <- data.frame(gene=all,cord=NA)
counts$cord <- ifelse(counts$gene %in% nsf.cord,1,0)
counts$epreg <- ifelse(counts$gene %in% nsf.epreg,1,0)
counts$lpreg <- ifelse(counts$gene %in% nsf.lpreg,1,0)
counts$plc <- ifelse(counts$gene %in% nsf.plc,1,0)
counts$plm <- ifelse(counts$gene %in% nsf.plm,1,0)
counts$total <- counts$cord + counts$epreg + counts$lpreg + counts$plc + counts$plm

venn.diagram(x=list(cord=nsf.cord,plm=nsf.plm,lpreg=nsf.lpreg,epreg=nsf.epreg,plc=nsf.plc),'venn-diagram-non-sfari.tiff',
             category=c(paste0('Cord Blood\n(n=',length(nsf.cord),')'),paste0('Maternal\nPlacenta\n(n=',length(nsf.plm),')'),paste0('Maternal Blood\nLate Preg (n=',length(nsf.lpreg),')'),paste0('Maternal Blood\nEarly Preg (n=',length(nsf.epreg),')'),paste0('Fetal Placenta\n(n=',length(nsf.plc),')')),
             mar=c(0.15,0.15,0.15,0.15),cat.dist=c(0.3,0.3,0.3,0.3,0.3),cex=0.8)


### permutation with random genes
library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(org.Hs.eg.db)
source("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/ASD_sens_functions.R")

sfari<-read.csv("/dcl01/NDEpi/data/Projects/InProgress/kbakulsk/SFARI-Gene_genes_export28-08-2017.csv")
ann <- missMethyl:::.getFlatAnnotation('450K')
ann.g <- unique(ann$symbol)

nperm <- 10
#set seeds for sake of reproducability of plot
seeds <- c(446,666,1453,1718,4353,325683,21460,24986,42,10)
random.permutations <- matrix(0,ncol=nperm,nrow=15)
rownames(random.permutations) <- as.character(c(0.0001, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04 , 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.99))
colnames(random.permutations) <- paste('p',1:nperm,sep='')

for(i in 1:nperm){
  set.seed(seeds[i])
  rando <- sample(ann.g,nrow(sfari))
  rand.cord <- gset.enrich(res.sv.cord,rando,missMethyl=F)
  random.permutations[,i] <- rand.cord$logp
}


library(RColorBrewer)
col <- brewer.pal(5, 'Set1')
pdf('permutation-enrichment.pdf')
par(mar=c(4.5,4.5,3.8,2))
plot(enrich.asd$cutoffs,enrich.asd$p.values,type='n',xlab='P-Value Cut-off in Single CpG Results to Select "Significant" Genes',
     ylab=expression('-Log'[10]*'(P-Value of SFARI Gene Set Enrichment)'),main='SFARI Enrichment',las=1)
axis(1, at=enrich.asd$cutoffs, labels=F, tck=0.04, lwd.ticks=1.2)
for(i in 1:ncol(random.permutations)){
  lines(enrich.asd$cutoffs,random.permutations[,i])
}
legend('topright',legend=c('Cord Blood','Maternal Blood Early Pregnancy','Maternal Blood Late Pregnancy','Fetal Placenta','Maternal Placenta','Random Gene Set'),
       col=c(col,'black'),lty=1, text.col=c(col,'black'), title='', title.col='black')
dev.off()


### enrichment with ntd
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/ntd-v-td-sv/res.sv.cord.rda")
enrich.ntd <- sfari.enrich(res.sv.cord, sparse.p=sparse,cut=enrich.epreg$peak.p.cutoff)

setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/")

library(RColorBrewer)
col <- brewer.pal(5, 'Set1')
pdf('ntd-enrichment.pdf')
par(mar=c(4.5,4.5,3.8,2))
plot(enrich.asd$cutoffs,enrich.asd$p.values,type='n',xlab='EWAS P-Value Cut-off',ylab=expression('-Log'[10]*'(P-Value)'),main='SFARI Gene Set Enrichment',las=1)
lines(enrich.asd$cutoffs,enrich.asd$p.values,col=col[1])
lines(enrich.ntd$cutoffs,enrich.ntd$p.values,col=col[2])
axis(1, at=enrich.asd$cutoffs, labels=F, tck=0.04, lwd.ticks=1.2)
legend('topright',legend=c('ASD vs TD','Non-TD vs TD'),
       col=c(col),lty=1, text.col=c(col))
dev.off()




#append all results to existing table
setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/")

####load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Results_CordBlood_ASD_NoSV-20170824.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/Results_CordBlood_ASD_NoSV-20171026.rda")
#load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_36MonthSRS/Results_CordBlood-SRS_NoSV-20170824.rda")
#load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_12MonthAOSI/Results_CordBlood-AOSI_NoSV-20170824.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/Results_MomBloodEarlyPreg-ASD_NoSV-20170824.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/Results_MomBloodLatePreg-ASD_NoSV-20170824.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/Results_PlacentaChild-ASD_NoSV-20170825.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/Results_PlacentaMom-ASD_NoSV-20170825.rda")

ss.hits.all <- res.nsv.asd[,]
# ss.hits.srs <- res.nsv.srs[rownames(ss.hits.all),c('logFC','P.Value')]
# names(ss.hits.srs) <- c('cord.srs.beta','cord.srs.p')
# identical(rownames(ss.hits.all),rownames(ss.hits.srs))
# ss.hits.all <- cbind(ss.hits.all,ss.hits.srs)
# rm(ss.hits.srs)
# ss.hits.aosi <- res.nsv.aosi[rownames(ss.hits.all),c('logFC','P.Value')]
# names(ss.hits.aosi) <- c('cord.aosi.beta','cord.aosi.p')
# identical(rownames(ss.hits.all),rownames(ss.hits.aosi))
# ss.hits.all <- cbind(ss.hits.all,ss.hits.aosi)
# rm(ss.hits.aosi)
ss.hits.epreg <- res.nsv.epreg[rownames(ss.hits.all),c('logFC','P.Value')]
names(ss.hits.epreg) <- c('matblood.earlypreg.asd.beta','matblood.earlypreg.asd.p')
identical(rownames(ss.hits.all),rownames(ss.hits.epreg))
ss.hits.all <- cbind(ss.hits.all,ss.hits.epreg)
rm(ss.hits.epreg)
ss.hits.lpreg <- res.nsv.lpreg[rownames(ss.hits.all),c('logFC','P.Value')]
names(ss.hits.lpreg) <- c('matblood.latepreg.asd.beta','matblood.latepreg.asd.p')
identical(rownames(ss.hits.all),rownames(ss.hits.lpreg))
ss.hits.all <- cbind(ss.hits.all,ss.hits.lpreg)
rm(ss.hits.lpreg)
ss.hits.plc <- res.nsv.plc[rownames(ss.hits.all),c('logFC','P.Value')]
names(ss.hits.plc) <- c('placenta.fetal.asd.beta','placenta.fetal.asd.p')
identical(rownames(ss.hits.all),rownames(ss.hits.plc))
ss.hits.all <- cbind(ss.hits.all,ss.hits.plc)
rm(ss.hits.plc)
ss.hits.plm <- res.nsv.plm[rownames(ss.hits.all),c('logFC','P.Value')]
names(ss.hits.plm) <- c('placenta.maternal.asd.beta','placenta.maternal.asd.p')
identical(rownames(ss.hits.all),rownames(ss.hits.plm))
ss.hits.all <- cbind(ss.hits.all,ss.hits.plm)
rm(ss.hits.plm)
#write.csv(ss.hits.all,file='ss.hits.20170830')

names(ss.hits.all)[1] <- "cord.asd.beta"
names(ss.hits.all)[4] <- "cord.asd.p"

#append chromosome and pos info
ss.hits.all <- read.csv('ss.hits.20170830')
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/450kannotation.rda")
anno <- an450k[rownames(an450k) %in% rownames(ss.hits.all),1:2]
anno <- anno[rownames(ss.hits.all),]

#annotated.results <- cbind(ss.hits.all[,1],anno,ss.hits.all[,2:ncol(ss.hits.all)])
annotated.results <- cbind(anno,ss.hits.all)

write.csv(annotated.results,file="ss.hits.all_annotated.csv")

annotated.results <- annotated.results[,-c(9:12)]
write.csv(annotated.results,file="ss.hits.asdpheno_annotated.csv")

write.csv(annotated.results,file="ss.hits.asd_annotated_v2.csv")

################################################
#placenta principal components maternal vs fetal
################################################

setwd('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/')

#load noob data
load("/dcl01/NDEpi/data/Projects/InProgress/kbakulsk/kbakulsk/EARLI/450k-round2/Noob-beta-All-Samples-EARLI-both-rounds.rda")

#load pd
load("/dcl01/NDEpi/data/MasterCohortData/EARLI/EARLIBiologicalSamples/450k/EARLI_450k_all1155samples_pd_with_ancestry.rda")
load("/dcl01/NDEpi/data/MasterCohortData/EARLI/EARLIBiologicalSamples/450k/earli.450k.pd.N=1225.20161117.rda")

pd.placenta <- pd[pd$Tissue=='placenta',]
pd.placenta$predictedSex <- pd.ancest[rownames(pd.placenta),'predictedSex']

beta.placenta <- noob.beta[,rownames(pd.placenta)]
rm(noob.beta,pd,pd.ancest)

#principal components
#pc.placenta <- prcomp(t(beta.placenta),center=T,scale.=F)
#save(pc.placenta,file='princomp-placenta-withsexchr.rda')
load('princomp-placenta-withsexchr.rda')

library(RColorBrewer)

pdf("PCA-placenta.pdf")

screeplot(pc.placenta, col="dodgerblue", xlab="Principal Components of Raw Beta Values", main="", cex.lab=1.3)

myColors <- brewer.pal(8, 'Set1')
palette(myColors)

plot.new()
par(mar=c(0, 0, 0, 0))
legend("bottom", levels(as.factor(pd.placenta$predictedSex)), fill=myColors, title="Principal Components by Sex")
pairs(pc.placenta$x[,1:8], col=as.factor(pd.placenta$predictedSex), labels=c("PC1", "PC2", "PC3", "PC4","PC5", "PC6", "PC7", "PC8"), pch=1, cex=0.5)

plot.new()
par(mar=c(0, 0, 0, 0))
legend("bottom", levels(as.factor(pd.placenta$Phenotype)), fill=myColors, title="Principal Components by Sex")
pairs(pc.placenta$x[,1:8], col=as.factor(pd.placenta$Phenotype), labels=c("PC1", "PC2", "PC3", "PC4","PC5", "PC6", "PC7", "PC8"), pch=1, cex=0.5)

dev.off()

#correlation of PCs to variables
num.pcs <- 10
vars <- c("Maternal/Fetal", "OffspringSex")
pvals <- matrix(0,nrow=length(vars),ncol=num.pcs)
colnames(pvals) <- colnames(pc.placenta$x)[1:num.pcs]
rownames(pvals) <- vars

for(pc in colnames(pvals)){
  pvals[1,pc] <- t.test(pc.placenta$x[,pc]~pd.placenta$Phenotype)$p.value
  pvals[2,pc] <- t.test(pc.placenta$x[,pc]~pd.placenta$predictedSex)$p.value
}

pvals

#output
# PC1          PC2          PC3          PC4          PC5
# Maternal/Fetal 1.878493e-06 1.159255e-01 6.322813e-01 2.161231e-02 0.0062262033
# OffspringSex   8.971537e-05 4.742107e-08 1.270237e-46 1.481137e-10 0.0001666712
# PC6       PC7          PC8         PC9      PC10
# Maternal/Fetal 0.0013189350 0.1042597 8.501319e-08 0.001981499 0.9817727
# OffspringSex   0.0009157649 0.3090217 2.878219e-01 0.117214923 0.8974489
###########



################################################
#enrichment using missMethyl
################################################
library(missMethyl)
library(org.Hs.eg.db)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
source("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/ASD_sens_functions.R")

# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/Results_CordBlood_ASD_NoSV-20171026.rda")
# #load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_36MonthSRS/Results_CordBlood-SRS_NoSV-20170824.rda")
# #load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_12MonthAOSI/Results_CordBlood-AOSI_NoSV-20170824.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/Results_MomBloodEarlyPreg-ASD_NoSV-20170824.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/Results_MomBloodLatePreg-ASD_NoSV-20170824.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/Results_PlacentaChild-ASD_NoSV-20170825.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/Results_PlacentaMom-ASD_NoSV-20170825.rda")

load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/asd-v-td-sv/res.sv.cord.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/asd-v-td-sv/res.sv.plm.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/asd-v-td-sv/res.sv.plc.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/asd-v-td-sv/res.sv.lpreg.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/asd-v-td-sv/res.sv.epreg.rda")

sfari <- read.csv("/dcl01/NDEpi/data/Projects/InProgress/kbakulsk/SFARI-Gene_genes_export28-08-2017.csv")
sfari.genes <- as.character(sfari$gene.symbol[!is.na(sfari$gene.symbol)])
sfari.genes <- mapIds(org.Hs.eg.db, keys=sfari.genes, column='ENTREZID', keytype='SYMBOL', multiVals='first')

cord.enrich <- gset.enrich(res.sv.cord,sfari.genes)
mbloodearly.enrich <- gset.enrich(res.sv.epreg,sfari.genes)
mbloodlate.enrich <- gset.enrich(res.sv.lpreg,sfari.genes)
plc.enrich <- gset.enrich(res.sv.plc,sfari.genes)
plm.enrich <- gset.enrich(res.sv.plm,sfari.genes)
saveRDS(cord.enrich,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/missmeth-cord.rda')
saveRDS(mbloodearly.enrich,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/missmeth-mbearly.rda')
saveRDS(mbloodlate.enrich,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/missmeth-mblate.rda')
saveRDS(plc.enrich,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/missmeth-plc.rda')
saveRDS(plm.enrich,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/missmeth-plm.rda')


#random set
ann <- flattenAnn('450K')
ann.g <- unique(ann$symbol)

# nperm <- 10
# #set seeds for sake of reproducability of plot
# seeds <- c(446,666,1453,1718,4353,325683,21460,24986,42,10)
# random.permutations <- matrix(0,ncol=nperm,nrow=15)
# rownames(random.permutations) <- as.character(c(0.0001, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04 , 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.99))
# colnames(random.permutations) <- paste('p',1:nperm,sep='')
# 
# set.seed(seeds[i])
set.seed(1453)
rando <- sample(ann.g,nrow(sfari))
rando <- mapIds(org.Hs.eg.db, keys=rando, column='ENTREZID', keytype='SYMBOL', multiVals='first')
rand.cord <- gset.enrich(res.sv.cord,rando)
# random.permutations[,i] <- rand.cord$logp


#
library(org.Hs.eg.db)
egGO2ALLEGS <- tryCatch(getFromNamespace("org.Hs.egGO2ALLEGS","org.Hs.eg.db"), error=function(e) FALSE)
universe <- intersect(AnnotationDbi::Lkeys(egGO2ALLEGS),genes.lung$universe)
AnnotationDbi::Lkeys(egGO2ALLEGS) <- universe
EG.GO <- AnnotationDbi::toTable(egGO2ALLEGS)[,c("gene_id","go_id","Ontology")]
d <- duplicated(EG.GO[,c("gene_id", "go_id")])
EG.GO <- EG.GO[!d, ]


#####################
#trouble shooting
sig <- rownames(res.nsv.asd[res.nsv.asd$P.Value<0.03,])
all <- rownames(res.nsv.asd)
gsameth(sig,all,collection=sfari.genes,array.type='450K')
out <- getMappedEntrezIDs(sig.cpg=sig,all.cpg=all,array.type='450K')
sorted.eg.sig <- out$sig.eg
eg.universe <- out$universe
freq_genes <- out$freq
test.de <- out$de
collection <- list(collection=sfari.genes)
collection <- lapply(collection, as.character)
# Make sure gene set collections don't have any NAs
collection <- lapply(collection, function(x) x[!is.na(x)])
# Make sure only collections with geneids present in universe are included
inUniv <- lapply(collection, function(x) sum(eg.universe %in% x))
collection <- collection[inUniv!=0]
library(limma)
estimatePWF <- function(D,bias)
{
  prior.prob <- bias
  o <- order(bias)
  prior.prob[o] <- tricubeMovingAverage(D[o],span=0.5)
  prior.prob
}
pwf <- estimatePWF(D=test.de,bias=as.vector(freq_genes))
InSet <- eg.universe %in% collection[[1]]
pw.red <- sum(pwf[InSet])/872
pw.white <- sum(pwf[!InSet])/(length(eg.universe)-872)
odds <- pw.red/pw.white

BiasedUrn::pWNCHypergeo(426,872,length(eg.universe)-872,length(sorted.eg.sig),odds,lower.tail=FALSE,precision=1E-20) 
phyper(q=319-0.5,m=length(sorted.eg.sig),n=length(eg.universe)-length(sorted.eg.sig),k=872,lower.tail=FALSE) 
#yay tests
########################################
library(RColorBrewer)
col <- brewer.pal(5, 'Set1')
pdf('sfari-enrichment-all-asd-models.pdf')
plot(enrich.asd$cutoffs,enrich.asd$p.values,type='n',xlab='EWAS P-Value Cut-off',ylab='-Log10 P-Value',main='SFARI Enrichment',las=1)
lines(enrich.asd$cutoffs,enrich.asd$p.values,col=col[1])
lines(enrich.epreg$cutoffs,enrich.epreg$p.values,col=col[2])
lines(enrich.lpreg$cutoffs,enrich.lpreg$p.values,col=col[3])
lines(enrich.plc$cutoffs,enrich.plc$p.values,col=col[4])
lines(enrich.plm$cutoffs,enrich.plm$p.values,col=col[5])
legend('topright',legend=c('Cord Blood','Maternal Blood Early Pregnancy','Maternal Blood Late Pregnancy','Fetal Placenta','Maternal Placenta'),
       col=col,lty=1, text.col=col, title='Tissue', title.col='black')
dev.off()

#######################################
## pick out cpgs in SFARI genes
library(missMethyl)
library(org.Hs.eg.db)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
source("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/ASD_sens_functions.R")

load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/Results_CordBlood_ASD_NoSV-20171026.rda")
sfari <- read.csv("/dcl01/NDEpi/data/Projects/InProgress/kbakulsk/SFARI-Gene_genes_export28-08-2017.csv")

ann <- missMethyl:::.flattenAnn('450K')
ann.sfari <- ann[ann$cpg %in% rownames(res.nsv.asd),]
ann.sfari <- ann.sfari[ann.sfari$symbol %in% sfari$gene.symbol,]



########################################
#placenta results for share
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/")

load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/Results_PlacentaChild-ASD_NoSV-20170825.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/Results_PlacentaMom-ASD_NoSV-20170825.rda")

names(res.nsv.plc)[1] <- c('Beta.Coef')
names(res.nsv.plm)[1] <- c('Beta.Coef')

res.nsv.plc <- res.nsv.plc[,1:6]
res.nsv.plm <- res.nsv.plm[,1:6]

#chr and pos
data(Locations)
plc <- data.frame(cbind(Locations[rownames(res.nsv.plc),1:2],res.nsv.plc))
plm <- data.frame(cbind(Locations[rownames(res.nsv.plm),1:2],res.nsv.plm))
plc$cpg <- rownames(plc)
plm$cpg <- rownames(plm)

#get std errors
fit.plc <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/plc.fit.rda")
fit.plm <- readRDS("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/plm.fit.rda")
std.error.plc <- sqrt(fit.plc$s2.post)*fit.plc$stdev.unscaled
std.error.plm <- sqrt(fit.plm$s2.post)*fit.plm$stdev.unscaled
plc$std.error <- std.error.plc[rownames(plc),2]
plm$std.error <- std.error.plm[rownames(plm),2]

plc.fmt <- plc[,c(10,1:2,4,11,6:8,5)]
plm.fmt <- plm[,c(10,1:2,4,11,6:8,5)]

setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/")
write.csv(plc.fmt,row.names=F,file='placenta_fetalside.csv')
write.csv(plm.fmt,row.names=F,file='placenta_maternalside.csv')


#############################################
### Volcano Plot
#############################################
source("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/ASD_sens_functions.R")

load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/asd-v-td-sv/res.sv.cord.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/asd-v-td-sv/res.sv.epreg.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/asd-v-td-sv/res.sv.lpreg.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/asd-v-td-sv/res.sv.plc.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/asd-v-td-sv/res.sv.plm.rda")

# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/Results_CordBlood_ASD_NoSV-20171026.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/Results_MomBloodEarlyPreg-ASD_NoSV-20170824.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/Results_MomBloodLatePreg-ASD_NoSV-20170824.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/Results_PlacentaChild-ASD_NoSV-20170825.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/Results_PlacentaMom-ASD_NoSV-20170825.rda")


setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Volcano")
# plot.volcano <- function(ss.hits,name,TIS='',pheno='ASD'){
#   #png(paste0(name,'.png'), width=10, height=10, units='in', res=300, type='cairo')
#   par(mar=c(5,5,3,3))
#   betas <- ss.hits$logFC
#   pvals <- -log(ss.hits$P.Value,10)
#   pcol <- ifelse(pvals < -log(0.05,10), 'darkslategrey', ifelse(betas>0,'cadetblue4','cadetblue3'))
#   plot(betas,pvals, las=1, xlab=paste0('Difference in Methylation ',pheno,' compared to TD'), ylab=expression('-Log'[10]*'(P-Value)'),main=TIS,
#        xlim=c(-0.4,0.4),ylim=c(0,8),cex.axis=2.1,cex.lab=2.1,cex.main=2.1,cex=0.7,pch=21,col=pcol)
#   abline(h=-log(0.05,10),col='blue')
#   abline(h=7,col='red')
#   #dev.off()
# }
# plot.volcdis <- function(ss.hits,name,TIS=''){
#   png(paste0(name,'.png'), width=10, height=10, units='in', res=300, type='cairo')
#   par(mar=c(5,5,3,3))
#   betas <- ss.hits$logFC
#   pctmeth <- ss.hits$AveExpr
#   pvals <- -log(ss.hits$P.Value,10)
#   pcol <- ifelse(pvals < -log(0.05,10), 'darkslategrey', ifelse(betas>0,'cadetblue4','cadetblue3'))
#   plot(pctmeth,pvals, las=1, xlab='Average Methylation', ylab='-Log 10 P-Value',main=TIS,
#        xlim=c(0,1),ylim=c(0,8),cex.axis=2.1,cex.lab=2.1,cex.main=2.1,cex=0.7,pch=21,col=pcol)
#   abline(h=-log(0.05,10),col='blue')
#   abline(h=5,col='red')
#   dev.off()
# }

plot.volcano <- function(ss.hits,main=NULL,pheno='ASD'){
  par(mar=c(5,5,3,3))
  
  #get betas and pvals from single sites results
  betas <- ss.hits$logFC*100
  pvals <- -log(ss.hits$P.Value,10)
  
  #assign colors
  pcol <- ifelse(pvals < -log(0.05,10), 'darkslategrey', ifelse(betas>0,'cadetblue4','cadetblue3'))
  
  #get percent hyper/hypo for sites above threshold
  total <- table(pcol)['cadetblue4'] + table(pcol)['cadetblue3']
  pct <- c(round(table(pcol)['cadetblue4']/total * 100,1),round(table(pcol)['cadetblue3']/total * 100,1))
  
  #plotting stuff
  plot(betas,pvals, las=1, main=main, 
       xlab=paste0('Difference in Methylation ',pheno,' compared to TD'), 
       ylab=expression('-Log'[10]*'(P-Value)'),
       xlim=c(-40,40),ylim=c(0,7),cex.axis=2.1,cex.lab=2.1,cex.main=2.1,cex=0.7,pch=21,col=pcol)
  abline(h=-log(0.05,10),col='blue')
  abline(h=7,col='red')
  text(30,6,paste0(pct[1],'%'),col='cadetblue4',cex=2.2)
  text(-30,6,paste0(pct[2],'%'),col='cadetblue3',cex=2.2)
}



# plot.volcano(res.nsv.asd,name='cord','Cord Blood')
# plot.volcano(res.nsv.epreg,name='mblood-epreg','Maternal Blood - Early Pregnancy')
# plot.volcano(res.nsv.lpreg,name='mblood-lpreg','Maternal Blood - Late Pregnancy')
# plot.volcano(res.nsv.plc,name='placenta-fetal','Placenta - Fetal Side')
# plot.volcano(res.nsv.plm,name='placenta-mom','Placenta - Maternal Side')


tiff("Volcano-Plots-Combined.tiff", width=18, height=10, units='in',res=300,type='cairo')
layout(matrix(c(1,1,2,2,3,3,4,5,5,6,6,7),2,6, byrow=T), widths=rep(1,6), heights=c(1,1))
  plot.volcano(res.sv.cord,'Cord Blood')
  plot.volcano(res.sv.epreg,'Maternal Blood - Early Pregnancy')
  plot.volcano(res.sv.lpreg,'Maternal Blood - Late Pregnancy')
  plot.new()
  plot.volcano(res.sv.plc,'Placenta - Fetal Side')
  plot.volcano(res.sv.plm,'Placenta - Maternal Side')
  plot.new()
dev.off()

plot.volcdis(res.sv.cord,name='sv-cord-ave','Cord Blood')

hy.table <- matrix(0,nrow=5,ncol=4)
colnames(hy.table) <- c('N Hypo','Pct Hypo','N Hyper','Pct Hyper')
rownames(hy.table) <- c('cord','epreg','lpreg','plc','plm')

# hits <- list(res.nsv.asd,res.nsv.epreg,res.nsv.lpreg,res.nsv.plc,res.nsv.plm)
hits <- list(res.sv.cord,res.sv.epreg,res.sv.lpreg,res.sv.plc,res.sv.plm)

counts <- sapply(hits,function(X){
  pcut <- X[X$P.Value<0.05,]
  hyp <- table(pcut$logFC > 0)
})
hy.table[,1] <- counts[1,]
hy.table[,3] <- counts[2,]
hy.table[,2] <- round((hy.table[,1] / (hy.table[,1]+hy.table[,3]))*100,2)
hy.table[,4] <- round((hy.table[,3] / (hy.table[,1]+hy.table[,3]))*100,2)

write.csv(hy.table,file='sv-hypo-hyper.csv')


load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/ntd-v-td-sv/res.sv.cord.rda")
tiff('Cord NTD vs TD.tiff', width=10, height=10, units='in', res=300, type='cairo')
plot.volcano(res.sv.cord,name='sv-cord-ntd','Cord Blood',pheno='non-TD')
dev.off()

hits <- list(res.sv.cord)
counts <- sapply(hits,function(X){
  pcut <- X[X$P.Value<0.05,]
  hyp <- table(pcut$logFC > 0)
})


### ntd vs asd beta
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/ntd-v-td-sv/res.sv.cord.rda")
res.ntd <- res.sv.cord
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/asd-v-td-sv/res.sv.cord.rda")
res.asd <- res.sv.cord

#order
res.ntd <- res.ntd[rownames(res.asd),]

setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/")
pdf("asd and ntd beta v beta.pdf")
smoothScatter(res.asd$logFC,res.ntd$logFC,xlab='ASD vs TD Betas',ylab='Non-TD vs TD Betas',colramp=colorRampPalette(c('white','slateblue3','slateblue4')), bandwidth=0.003)
dev.off()
cor(res.asd$logFC,res.ntd$logFC)

#############################################
### Top Hits
#############################################
source("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/ASD_sens_functions.R")

load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/asd-v-td-sv/res.sv.cord.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/asd-v-td-sv/res.sv.epreg.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/asd-v-td-sv/res.sv.lpreg.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/asd-v-td-sv/res.sv.plc.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/asd-v-td-sv/res.sv.plm.rda")

setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/asd-v-td-sv/")
top.cord <- res.sv.cord[res.sv.cord$P.Value<1e-4,]
top.cord <- fmt.tophits(top.cord)
write.csv(top.cord,file='cord.tophits.csv')



#############################################
#global
#############################################
source("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/ASD_sens_functions.R")

setwd('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/')
pd.cord <- readRDS('pdsv.cord.td.rds')
beta.cord <- readRDS('noob.cord.rds')

setwd('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/')
pd.epreg <- readRDS('pdsv.epreg.td.rds')
beta.epreg <- readRDS('noob.epreg.rds')

setwd('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/')
pd.lpreg <- readRDS('pdsv.lpreg.td.rds')
beta.lpreg <- readRDS('noob.lpreg.rds')

setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/")
pd.plc <- readRDS('pdsv.plc.td.rds')
beta.plc <- readRDS('noob.plc.rds')

setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/")
pd.plm <- readRDS('pdsv.plm.td.rds')
beta.plm <- readRDS('noob.plm.rds')

source("/dcl01/NDEpi/data/Projects/InProgress/jdou/Folate/EARLI_folate_intake_functions.R")

global.vars <- function(pd, beta){
  pd$meanDNAm <- colMeans(beta)*100
  pd$meanDNAm.sea <- colMeans(genomic.region(beta,region='OpenSea',anno='450k'))*100
  pd$meanDNAm.shore <- colMeans(genomic.region(beta,region='Shore',anno='450k'))*100
  pd$meanDNAm.shelf <- colMeans(genomic.region(beta,region='Shelf',anno='450k'))*100
  pd$meanDNAm.island <- colMeans(genomic.region(beta,region='Island',anno='450k'))*100
  pd
}
pd.cord <- global.vars(pd.cord, beta.cord[,rownames(pd.cord)])
pd.epreg <- global.vars(pd.epreg, beta.epreg[,rownames(pd.epreg)])
pd.lpreg <- global.vars(pd.lpreg, beta.lpreg[,rownames(pd.lpreg)])
pd.plc <- global.vars(pd.plc, beta.plc[,rownames(pd.plc)])
pd.plm <- global.vars(pd.plm, beta.plm[,rownames(pd.plm)])

#function expects 0 and 1
td1 <- function(pd){
  pd$td1 <- ifelse(pd$td=='TD',0,1)
  return(pd)
}
pd.cord <- td1(pd.cord)
pd.epreg <- td1(pd.epreg)
pd.lpreg <- td1(pd.lpreg)
pd.plc <- td1(pd.plc)
pd.plm <- td1(pd.plm)

cord.gl <- global.island(autosomes(beta.cord[,rownames(pd.cord)]),pd.cord,'td1')
epreg.gl <- global.island(autosomes(beta.epreg[,rownames(pd.epreg)]),pd.epreg,'td1')
lpreg.gl <- global.island(autosomes(beta.lpreg[,rownames(pd.lpreg)]),pd.lpreg,'td1')
plc.gl <- global.island(autosomes(beta.plc[,rownames(pd.plc)]),pd.plc,'td1')
plm.gl <- global.island(autosomes(beta.plm[,rownames(pd.plm)]),pd.plm,'td1')

globalRegions <- function(pd, x, adj=NULL, summary.name){
  #set up blank results
  regions <- c('meanDNAm','meanDNAm.sea','meanDNAm.shore','meanDNAm.shelf','meanDNAm.island')
  fits <- list()
  fits[regions] <- NA
  results <- matrix(NA, nrow=5, ncol=4)
  rownames(results) <- regions
  colnames(results) <- c('est', 'std_err', 't', 'p')
  results <- data.frame(results)
  
  #do linear model for each region
  for(region in regions){
    eqn <- paste0(region, '~', x, '+', adj)
    fits[[region]] <- lm(as.formula(eqn), data=pd)
    results[region,] <- summary(fits[[region]])$coefficients[2,]
  }
  
  results$lower <- results$est - 1.96*results$std_err
  results$upper <- results$est + 1.96*results$std_err
  
  return(results)
}

cord.gl <- globalRegions(pd.cord, 'td1', 'sv1+sv2+sv3')
epreg.gl <- globalRegions(pd.epreg, 'td1', 'sv1+sv2+sv3')
lpreg.gl <- globalRegions(pd.lpreg, 'td1', 'sv1+sv2+sv3')
plc.gl <- globalRegions(pd.plc, 'td1', 'sv1+sv2+sv3')
plm.gl <- globalRegions(pd.plm, 'td1', 'sv1+sv2')

# setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/")
# library(metafor)
# 
# pdf("Global-5tissue.pdf")
# global <- rbind(cord.gl,epreg.gl,lpreg.gl,plc.gl,plm.gl)
# par(xpd=NA,mar=c(4,4,2,2))
# forest(global[,2], ci.lb=global$lower, ci.ub=global$upper,
#        annotate=FALSE, slab=global$feature,
#        rows=c(32:28,25:21,18:14,11:7,4:0), ylim=c(0,36), refline=NULL,
#        xlab="Difference in mean percent methylation",alim=c(-2.0,1.0),xlim=c(-2.3,1.0))
# segments(0,-1,0,34,lty='dashed')
# par(font=1)
# #text(0.25, 37, "")
# text(0, c(33,26,19,12,5), pos=4, c("Cord Blood","Maternal Blood Early Preg","Maternal Blood Late Preg","Placenta Fetal","Placenta Maternal"))
# dev.off()
# 
# 
# ### modeling per person DNAm
# extract <- function(lm){
#   coef <- lm$coefficients['tdASD','Estimate']
#   se <- lm$coefficients['tdASD','Std. Error']
#   p <- lm$coefficients['tdASD','Pr(>|t|)']
#   c(coef,coef+1.96*se,coef-1.96*se,p)
# }


setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/")
library(metafor)

pdf("Global-5tissue-model_sv_adjusted.pdf")
global <- rbind(cord.gl,epreg.gl,lpreg.gl,plc.gl,plm.gl)
global <- data.frame(global)
global$feature <- c("All","Sea","Shelf","Shore","Island")
par(xpd=NA,mar=c(4,4,2,2))
forest(global$est, ci.lb=global$lower, ci.ub=global$upper,
       annotate=FALSE, slab=global$feature,
       rows=c(32:28,25:21,18:14,11:7,4:0), ylim=c(0,36), refline=NULL,
       xlab="Difference in mean percent methylation",alim=c(-1.5,1.0),xlim=c(-1.5,1.0))
segments(0,-1,0,34,lty='dashed')
par(font=1)
#text(0.25, 37, "")
text(0, c(33,26,19,12,5), pos=4, c("Cord Blood","Maternal Blood Early Preg","Maternal Blood Late Preg","Placenta Fetal","Placenta Maternal"))
dev.off()

# feature    meandiff        upper        lower           p  N
# 1      All  0.05564510  0.248258364 -0.136968156 0.573019401 76
# 2      Sea  0.03990633  0.286990411 -0.207177754 0.752507989 76
# 3    Shelf  0.03762165  0.283049276 -0.207805969 0.764712819 76
# 4    Shore  0.06994990  0.300622273 -0.160722478 0.554162458 76
# 5   Island  0.06850068  0.242735879 -0.105734522 0.443515387 76

# 6      All  0.06163219  0.157824156 -0.034559783 0.212932686 83
# 7      Sea  0.09726751  0.284608042 -0.090073017 0.311997484 83
# 8    Shelf  0.07095306  0.282228523 -0.140322408 0.512328075 83
# 9    Shore  0.05525369  0.151851036 -0.041343649 0.265676479 83
# 10  Island  0.02274199  0.154679261 -0.109195289 0.736387895 83

# 11     All  0.17815152  0.337355459  0.018947584 0.034649662 42
# 12     Sea  0.17782236  0.396908579 -0.041263859 0.120153975 42
# 13   Shelf  0.11377008  0.311068067 -0.083527902 0.265658780 42
# 14   Shore  0.19323622  0.379697345  0.006775099 0.049463218 42
# 15  Island  0.18689125  0.463743698 -0.089961195 0.193916904 42

# 16     All -0.26966320 -0.057127556 -0.482198847 0.016418584 53
# 17     Sea -0.38601041 -0.120152694 -0.651868121 0.006496455 53
# 18   Shelf -0.42961094 -0.142442686 -0.716779185 0.005143901 53
# 19   Shore -0.24560660 -0.007856044 -0.483357155 0.048473782 53
# 20  Island -0.10567716  0.131786980 -0.343141303 0.387416329 53

# 21     All  0.03414375  0.277233182 -0.208945680 0.784245781 53
# 22     Sea  0.06829151  0.333516799 -0.196933772 0.616052464 53
# 23   Shelf -0.09400003  0.111777366 -0.299777426 0.374984420 53
# 24   Shore  0.01897817  0.346640777 -0.308684427 0.910079831 53
# 25  Island  0.04533626  0.406401843 -0.315729327 0.806631067 53


### GAMP
library(GAMP)
#td+PC1+PC2+predictedSex+Maternal.Age+Round+Gran+nRBC
cdf.test <- function(beta,pd,var,adj){
  adj.mat <- model.matrix(as.formula(paste0('~',adj)))
  # adj.mat <- model.matrix(~pd$Maternal.Age+
  #                           factor(pd$Round)+
  #                           factor(pd$predictedSex)+
  #                           pd$PC1+
  #                           pd$PC2+
  #                           pd$Gran+
  #                           pd$nRBC)
  beta <- beta[,rownames(pd)]
  
  #set up matrix for results
  res <- matrix(0,nrow=5,ncol=2)
  rownames(res) <- c('All','OpenSea','Shelf','Shore','Island')
  colnames(res) <- c('unadj','adj')
  
  res['All','unadj'] <- TestCDFs(beta, var, X=NULL, outcomeType='D') 
  res['OpenSea','unadj'] <- TestCDFs(genomic.region(beta,'OpenSea'), var, X=NULL, outcomeType='D') 
  res['Shelf','unadj'] <- TestCDFs(genomic.region(beta,'Shelf'), var, X=NULL, outcomeType='D') 
  res['Shore','unadj'] <- TestCDFs(genomic.region(beta,'Shore'), var, X=NULL, outcomeType='D') 
  res['Island','unadj'] <- TestCDFs(genomic.region(beta,'Island'), var, X=NULL, outcomeType='D') 
  
  if(is.null(adj)){
    res[,'adj'] <- NA
  }else{
    res['All','adj'] <- TestCDFs(beta, var, X=adj.mat, outcomeType='D') 
    res['OpenSea','adj'] <- TestCDFs(genomic.region(beta,'OpenSea'), var, X=adj.mat, outcomeType='D') 
    res['Shelf','adj'] <- TestCDFs(genomic.region(beta,'Shelf'), var, X=adj.mat, outcomeType='D') 
    res['Shore','adj'] <- TestCDFs(genomic.region(beta,'Shore'), var, X=adj.mat, outcomeType='D') 
    res['Island','adj'] <- TestCDFs(genomic.region(beta,'Island'), var, X=adj.mat, outcomeType='D')
  }
  
  return(res)
}

beta.cord <- autosomes(beta.cord[,rownames(pd.cord)])
cord.gamp <- cdf.test(beta.cord,pd.cord,pd.cord$td,'pd$sv1+pd$sv2+pd$sv3')
#         unadj       adj
# All     0.1967010 0.5699555
# OpenSea 0.1728253 0.4180258
# Shelf   0.1790831 0.3112726
# Shore   0.2454935 0.6510004
# Island  0.4776278 0.7498439

beta.epreg <- autosomes(beta.epreg[,rownames(pd.epreg)])
epreg.gamp <- cdf.test(beta.epreg,pd.epreg,pd.epreg$td,'pd$sv1+pd$sv2+pd$sv3')
#         unadj       adj
# All     0.2914015 0.6340562
# OpenSea 0.3236482 0.4020719
# Shelf   0.3925153 0.4854970
# Shore   0.1848151 0.7891814
# Island  0.3685294 0.7705457

beta.lpreg <- autosomes(beta.lpreg[,rownames(pd.lpreg)])
lpreg.gamp <- cdf.test(beta.lpreg,pd.lpreg,pd.lpreg$td,'pd$sv1+pd$sv2+pd$sv3')
#         unadj       adj
# All     0.6543893 0.2661575
# OpenSea 0.7682701 0.4679605
# Shelf   0.8030803 0.7048855
# Shore   0.4056164 0.1585462
# Island  0.5155863 0.2665552

beta.plc <- autosomes(beta.plc[,rownames(pd.plc)])
plc.gamp <- cdf.test(beta.plc,pd.plc,pd.plc$td,'pd$sv1+pd$sv2+pd$sv3')
#         unadj         adj
# All     1.0000000 0.014136198
# OpenSea 0.7902639 0.004982003
# Shelf   0.8118266 0.008259423
# Shore   1.0000000 0.019178217
# Island  0.7664668 0.417194947

beta.plm <- autosomes(beta.plm[,rownames(pd.plm)])
plm.gamp <- cdf.test(beta.plm,pd.plm,as.numeric(pd.plm$td)-1,'pd$sv1+pd$sv2')
#         unadj       adj
# All     0.4123999 0.2623669
# OpenSea 0.4411303 0.2666772
# Shelf   0.2368604 0.1418353
# Shore   0.4025310 0.2703912
# Island  0.3341060 0.4737558

library(openxlsx)
write.xlsx(list(cord_blood=cord.gamp,
                mat_blood_epreg=epreg.gamp,
                mat_blood_lpreg=lpreg.gamp,
                placenta_fetal=plc.gamp,
                placenta_maternal=plm.gamp),
           file="Cumulative Density Function Differences.xlsx",
           row.names=T)

################################
# correlation plot
################################
library(corrplot)
setwd('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/')

# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/Results_CordBlood_ASD_NoSV-20171026.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/Results_MomBloodEarlyPreg-ASD_NoSV-20170824.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/Results_MomBloodLatePreg-ASD_NoSV-20170824.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/Results_PlacentaChild-ASD_NoSV-20170825.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/Results_PlacentaMom-ASD_NoSV-20170825.rda")

load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/asd-v-td-sv/res.sv.cord.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/asd-v-td-sv/res.sv.epreg.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/asd-v-td-sv/res.sv.lpreg.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/asd-v-td-sv/res.sv.plc.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/asd-v-td-sv/res.sv.plm.rda")

# coefficients <- cbind(cord=res.nsv.asd[,1],
#                       mombloodearly=res.nsv.epreg[rownames(res.nsv.asd),1],
#                       mombloodlate=res.nsv.lpreg[rownames(res.nsv.asd),1],
#                       placentafetal=res.nsv.plc[rownames(res.nsv.asd),1],
#                       placentamom=res.nsv.plm[rownames(res.nsv.asd),1])
coefficients <- cbind(cord=res.sv.cord[,1],
                      mombloodearly=res.sv.epreg[rownames(res.sv.cord),1],
                      mombloodlate=res.sv.lpreg[rownames(res.sv.cord),1],
                      placentafetal=res.sv.plc[rownames(res.sv.cord),1],
                      placentamom=res.sv.plm[rownames(res.sv.cord),1])
M <- cor(coefficients)
rownames(M) <- c('Cord Blood','Mom Blood Early Preg','Mom Blood Late Preg','Placenta Fetal','Placenta Maternal')
colnames(M) <- rownames(M)


pdf('sv-corrplot-5tissue.pdf')
corrplot(M,type='upper',method='ellipse',addCoef.col='black',diag=F)
dev.off()

### invariant probes

#cord blood
setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/")
pd.cord <- readRDS('pdsv.cord.td.rds')
beta.cord <- readRDS('noob.cord.rds')
beta.cord <- beta.cord[,rownames(pd.cord)]

#early preg mblood
setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/")
pd.epreg <- readRDS('pdsv.epreg.td.rds')
beta.epreg <- readRDS('noob.epreg.rds')
beta.epreg <- beta.epreg[,rownames(pd.epreg)]

#late preg mblood
setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/")
pd.lpreg <- readRDS('pdsv.lpreg.td.rds')
beta.lpreg <- readRDS('noob.lpreg.rds')
beta.lpreg <- beta.lpreg[,rownames(pd.lpreg)]

#placenta child
setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/")
pd.plc <- readRDS('pdsv.plc.td.rds')
beta.plc <- readRDS('noob.plc.rds')
beta.plc <- beta.plc[,rownames(pd.plc)]

#placenta mom
setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/")
pd.plm <- readRDS('pdsv.plm.td.rds')
beta.plm <- readRDS('noob.plm.rds')
beta.plm <- beta.plm[,rownames(pd.plm)]

#varainces of betas
probe.var <- matrix(-1,nrow=nrow(beta.cord),ncol=5)
colnames(probe.var) <- c('cord','epreg','lpreg','plc','plm')
rownames(probe.var) <- rownames(beta.cord)
probe.var <- data.frame(probe.var)

probe.var$cord <- apply(beta.cord,1,var)
probe.var$epreg <- apply(beta.epreg,1,var)
probe.var$lpreg <- apply(beta.lpreg,1,var)
probe.var$plc <- apply(beta.plc,1,var)
probe.var$plm <- apply(beta.plm,1,var)

setwd('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/')
save(probe.var,file='probe.var.rda')

pdf('probe.var.pdf')
hist(probe.var$cord,breaks=20)
hist(probe.var$epreg,breaks=20)
hist(probe.var$lpreg,breaks=20)
hist(probe.var$plc,breaks=20)
hist(probe.var$plm,breaks=20)
dev.off()



################################
# results tables
################################
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(Locations)
data(Other)
setwd('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/')

load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/asd-v-td-sv/res.sv.cord.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/asd-v-td-sv/res.sv.epreg.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/asd-v-td-sv/res.sv.lpreg.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/asd-v-td-sv/res.sv.plc.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/asd-v-td-sv/res.sv.plm.rda")

asd_td_res <- res.sv.cord
asd_td_res <- asd_td_res[,c(1,4)]
colnames(asd_td_res) <- c("cord_blood_effect_estimate", "cord_blood_p_value")
asd_td_res$chr <- Locations[rownames(asd_td_res),'chr']
asd_td_res$pos <- Locations[rownames(asd_td_res),'pos']
asd_td_res$gene <- Other[rownames(asd_td_res),'UCSC_RefGene_Name']
asd_td_res$gene <- paste0(' ',asd_td_res$gene)
asd_td_res <- asd_td_res[,c(3,4,5,1,2)]
asd_td_res$maternal_blood_early_preg_effect_estimate <- res.sv.epreg[rownames(asd_td_res),'logFC']
asd_td_res$maternal_blood_early_preg_p_value <- res.sv.epreg[rownames(asd_td_res),'P.Value']
asd_td_res$maternal_blood_late_preg_effect_estimate <- res.sv.lpreg[rownames(asd_td_res),'logFC']
asd_td_res$maternal_blood_late_preg_p_value <- res.sv.lpreg[rownames(asd_td_res),'P.Value']
asd_td_res$placenta_fetal_late_preg_effect_estimate <- res.sv.plc[rownames(asd_td_res),'logFC']
asd_td_res$placenta_fetal_late_preg_p_value <- res.sv.plc[rownames(asd_td_res),'P.Value']
asd_td_res$placenta_maternal_late_preg_effect_estimate <- res.sv.plm[rownames(asd_td_res),'logFC']
asd_td_res$placenta_maternal_late_preg_p_value <- res.sv.plm[rownames(asd_td_res),'P.Value']
write.csv(asd_td_res, file="Supplementary Table 3 - ASD vs TD Single Site Results.csv")

rm(res.sv.cord, res.sv.epreg, res.sv.lpreg, res.sv.plm, res.sv.plc)

load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/ntd-v-td-sv/res.sv.cord.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/ntd-v-td-sv/res.sv.epreg.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/ntd-v-td-sv/res.sv.lpreg.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/ntd-v-td-sv/res.sv.plc.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/ntd-v-td-sv/res.sv.plm.rda")

ntd_td_res <- res.sv.cord
ntd_td_res <- ntd_td_res[,c(1,4)]
colnames(ntd_td_res) <- c("cord_blood_effect_estimate", "cord_blood_p_value")
ntd_td_res$chr <- Locations[rownames(ntd_td_res),'chr']
ntd_td_res$pos <- Locations[rownames(ntd_td_res),'pos']
ntd_td_res$gene <- Other[rownames(ntd_td_res),'UCSC_RefGene_Name']
ntd_td_res$gene <- paste0(' ',ntd_td_res$gene)
ntd_td_res <- ntd_td_res[,c(3,4,5,1,2)]
ntd_td_res$maternal_blood_early_preg_effect_estimate <- res.sv.epreg[rownames(ntd_td_res),'logFC']
ntd_td_res$maternal_blood_early_preg_p_value <- res.sv.epreg[rownames(ntd_td_res),'P.Value']
ntd_td_res$maternal_blood_late_preg_effect_estimate <- res.sv.lpreg[rownames(ntd_td_res),'logFC']
ntd_td_res$maternal_blood_late_preg_p_value <- res.sv.lpreg[rownames(ntd_td_res),'P.Value']
ntd_td_res$placenta_fetal_late_preg_effect_estimate <- res.sv.plc[rownames(ntd_td_res),'logFC']
ntd_td_res$placenta_fetal_late_preg_p_value <- res.sv.plc[rownames(ntd_td_res),'P.Value']
ntd_td_res$placenta_maternal_late_preg_effect_estimate <- res.sv.plm[rownames(ntd_td_res),'logFC']
ntd_td_res$placenta_maternal_late_preg_p_value <- res.sv.plm[rownames(ntd_td_res),'P.Value']
write.csv(ntd_td_res, file="Supplementary Table 4 - Non-TD vs TD Single Site Results.csv")

################################
# meQTL
################################
library(missMethyl)
library(annotate)
library(org.Hs.eg.db)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
setwd("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/")
source("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/ASD_sens_functions.R")

load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/asd-v-td-sv/res.sv.cord.rda")

# cord <- read.csv(gzfile("/dcl01/NDEpi/data/Projects/InProgress/kbakulsk/cordmeqtls.csv.gz"),header=TRUE)
# meqtl <- names(table(cord$CpG_ID))
# save(meqtl,file='meqlt.rda')
load('meqlt.rda')

#top cpgs enriched for meqtls?
cpg.meqtl <- function(p,rand=F){
  top <- res.nsv.asd[res.nsv.asd$P.Value<p,]
  if(rand){
    top <- res.nsv.asd[sample(x=1:nrow(res.nsv.asd),size=nrow(top)),]
  }
  
  tab <- matrix(0,nrow=2,ncol=2)
  tab[1,c(2,1)] <- table(rownames(top) %in% meqtl)
  tab[2,1] <- table(meqtl %in% rownames(top))[1]
  tab[2,2] <- nrow(res.nsv.asd) - tab[2,1] - tab[1,1] - tab[1,2]
  return(tab)
}

p1e4 <- cpg.meqtl(0.0001)
p1e3 <- cpg.meqtl(0.001)
p1e2 <- cpg.meqtl(0.01)
p3e2 <- cpg.meqtl(0.03)
p5e2 <- cpg.meqtl(0.05)



#meqtl genes in sfari genes?
universe <- getMappedEntrezIDs(rownames(res.nsv.asd),array.type='450K')
universe <- universe$universe
universe <- getSYMBOL(universe,data='org.Hs.eg')
  
meqtl.gene <- getMappedEntrezIDs(meqtl,array.type='450K')
meqtl.gene <- getSYMBOL(meqtl.gene$sig.eg,data='org.Hs.eg')

sfari<-read.csv("/dcl01/NDEpi/data/Projects/InProgress/kbakulsk/SFARI-Gene_genes_export28-08-2017.csv")

tab <- matrix(0,nrow=2,ncol=2)
tab[1,c(2,1)] <- table(meqtl.gene %in% sfari$gene.symbol)
tab[2,1] <- table(sfari$gene.symbol %in% meqtl.gene)[1]
tab[2,2] <- length(universe) - tab[2,1] - tab[1,1] - tab[1,2]

tab
chi <- chisq.test(tab)
chi; chi$expected

#enrichment restricted to meqtls
cpg.restrict <- res.sv.cord[rownames(res.sv.cord) %in% meqtl,]
enrich.asd.me <- sfari.enrich(cpg.restrict, sparse.p=T, cut=0.05)
enrich.asd.full <- sfari.enrich(res.sv.cord, sparse.p=T, cut=0.05)

library(RColorBrewer)
col <- brewer.pal(5, 'Set1')
pdf('sfari-enrichment-cord-meqtl.pdf')
par(mar=c(4.5,4.5,3.8,2))
plot(enrich.asd.full$cutoffs,enrich.asd.full$p.values,type='n',xlab='EWAS P-Value Cut-off',ylab=expression('-Log'[10]*'(P-Value)'),main='SFARI Enrichment',las=1)
axis(1, at=enrich.asd.full$cutoffs, labels=F, tck=0.04, lwd.ticks=1.2)
lines(enrich.asd.full$cutoffs,enrich.asd.full$p.values,col=col[1])
lines(enrich.asd.me$cutoffs,enrich.asd.me$p.values,col=col[2])
legend('topright',legend=c('Full','meQTL'),
       col=col,lty=1, text.col=col, title.col='black')
dev.off()

#sfari vs meqtl vs cord
sfari.gene <- as.character(sfari$gene.symbol)
meqtl.gene <- getMappedEntrezIDs(meqtl,array.type='450K')
meqtl.gene <- getSYMBOL(meqtl.gene$sig.eg,data='org.Hs.eg')
meqtl.gene <- as.character(meqtl.gene)


annot<-read.csv("/dcl01/NDEpi/data/Projects/InProgress/kbakulsk/kbakulsk/NCSv/NoControls/Resid/annot-loc-symbol-entrez09-24-13")
dat.annot<-merge(res.nsv.asd, annot, by.x="row.names",by.y="X", sort=F)
dat.na<-dat.annot[!(is.na(dat.annot$entrezids)),]
dat.o<-dat.na[order(dat.na$entrezids, dat.na$P.Value),]
dat.o$index<-ave(rep(NA, nrow(dat.o)), dat.o$entrezids, FUN=seq_along) 
dat.low.p<-dat.o[dat.o$index==1,]
dat.low.p <- dat.low.p[dat.low.p$P.Value<0.03,]
cord.gene <- as.character(dat.low.p$symbols)

library(VennDiagram)
library(RColorBrewer)
col <- brewer.pal(3, 'Set1')
venn.diagram(list(SFARI=sfari.gene,meQTL=meqtl.gene,cord=cord.gene),filename='venn-diagram-sfari-meqtl-cord.tiff')


################################
# pathway
################################
library(missMethyl)
setwd('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/')

# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/Results_CordBlood_ASD_NoSV-20171026.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/Results_MomBloodEarlyPreg-ASD_NoSV-20170824.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/Results_MomBloodLatePreg-ASD_NoSV-20170824.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/Results_PlacentaChild-ASD_NoSV-20170825.rda")
# load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/Results_PlacentaMom-ASD_NoSV-20170825.rda")

load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/asd-v-td-sv/res.sv.cord.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_EarlyPreg_ASD/asd-v-td-sv/res.sv.epreg.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Mom_LatePreg_ASD/asd-v-td-sv/res.sv.lpreg.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaChild_ASD/asd-v-td-sv/res.sv.plc.rda")
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/PlacentaMom_ASD/asd-v-td-sv/res.sv.plm.rda")

load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/ntd-v-td-sv/res.sv.cord.rda")


library('IlluminaHumanMethylation450kanno.ilmn12.hg19')

p <- 0.001

top.cord <- as.character(rownames(res.sv.cord[res.sv.cord$P.Value < p,]))
genes.cord <- getMappedEntrezIDs(top.cord,all.cpg=rownames(res.sv.cord),array.type='450K')
path.cord <- gometh(sig.cpg=top.cord,all.cpg=rownames(res.sv.cord))
path.cord <- path.cord[order(path.cord$P.DE),]
path.cord <- path.cord[path.cord$Ont=='BP',]
write.csv(path.cord,file='/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/asd-v-td-sv/cord.asd.path.csv')

top.epreg <- as.character(rownames(res.sv.epreg[res.sv.epreg$P.Value < p,]))
path.epreg <- gometh(sig.cpg=top.epreg,all.cpg=rownames(res.sv.epreg))
path.epreg <- path.epreg[order(path.epreg$P.DE),]
path.epreg <- path.epreg[path.epreg$Ont=='BP',]

top.lpreg <- as.character(rownames(res.sv.lpreg[res.sv.lpreg$P.Value < p,]))
path.lpreg <- gometh(sig.cpg=top.lpreg,all.cpg=rownames(res.sv.lpreg))
path.lpreg <- path.lpreg[order(path.lpreg$P.DE),]
path.lpreg <- path.lpreg[path.lpreg$Ont=='BP',]

top.plc <- as.character(rownames(res.sv.plc[res.sv.plc$P.Value < p,]))
path.plc <- gometh(sig.cpg=top.plc,all.cpg=rownames(res.sv.plc))
path.plc <- path.plc[order(path.plc$P.DE),]
path.plc <- path.plc[path.plc$Ont=='BP',]

top.plm <- as.character(rownames(res.sv.plm[res.sv.plm$P.Value < p,]))
path.plm <- gometh(sig.cpg=top.plm,all.cpg=rownames(res.sv.plm))
path.plm <- path.plm[order(path.plm$P.DE),]
path.plm <- path.plm[path.plm$Ont=='BP',]

#rank order
path.cord$rank <- c(1:nrow(path.cord))
path.epreg$rank <- c(1:nrow(path.epreg))
path.lpreg$rank <- c(1:nrow(path.lpreg))
path.plc$rank <- c(1:nrow(path.plc))
path.plm$rank <- c(1:nrow(path.plm))

#order by cord ranks
path.epreg <- path.epreg[rownames(path.cord),]
path.lpreg <- path.lpreg[rownames(path.cord),]
path.plc <- path.plc[rownames(path.cord),]
path.plm <- path.plm[rownames(path.cord),]

#sum ranks
path <- path.cord[,c(1,5,7)]
names(path)[c(2,3)] <- c('cord.p','cord.rank')
path$epreg.p <- path.epreg$P.DE
path$epreg.rank <- path.epreg$rank
path$lpreg.p <- path.lpreg$P.DE
path$lpreg.rank <- path.lpreg$rank
path$plc.p <- path.plc$P.DE
path$plc.rank <- path.plc$rank
path$plm.p <- path.plm$P.DE
path$plm.rank <- path.plm$rank
path$sum.rank <- path$cord.rank + path$epreg.rank + path$lpreg.rank + path$plc.rank + path$plm.rank
path <- path[order(path$sum.rank),]

#save
write.csv(path,file='sv-path-ranksum.csv')


#########################################################
# boxplots
#########################################################
/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/
setwd('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/')
load("/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/asd-v-td-sv/res.sv.cord.rda")

pd.cord <- readRDS('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/pdsv.cord.td.rds')
beta.cord <- readRDS('/dcl01/NDEpi/data/Projects/InProgress/jdou/ASD_EARLI_sens/Cord_ASD/noob.cord.rds')
beta.cord <- beta.cord[,rownames(pd.cord)]

#top 6 CpGs
top6 <- rownames(res.sv.cord)[1:6]

#set up colors
col <- c('skyblue3','indianred4')
colight1 <- rgb(t(col2rgb(col)),maxColorValue = 255, alpha=220)
colight2 <- rgb(t(col2rgb(col)),maxColorValue = 255, alpha=20)

tiff('boxplots.tiff',height=10, width=18, units='in', res=300, type='cairo')
#layout(matrix(c(1,2,3,4,5,6),2,3, byrow=T), widths=rep(1,3), heights=c(1,1))
par(mfrow=c(2,3))
boxplot(beta.cord[top6[1],]*100~pd.cord$td, ylab='Percent Methylation', main=top6[1], cex=2,cex.main=2,cex.lab=1.2,cex.axis=1.5, col=colight2, medcol=colight1, whiskcol=colight1, staplecol=colight1, boxcol=colight1, outcol='white')
stripchart(beta.cord[top6[1],]*100~pd.cord$td, vertical=T, method='jitter',add=T,pch=20,col=col,cex=2)
boxplot(beta.cord[top6[2],]*100~pd.cord$td, ylab='Percent Methylation', main=top6[2], cex=2,cex.main=2,cex.lab=1.2,cex.axis=1.5, col=colight2, medcol=colight1, whiskcol=colight1, staplecol=colight1, boxcol=colight1, outcol='white')
stripchart(beta.cord[top6[2],]*100~pd.cord$td, vertical=T, method='jitter',add=T,pch=20,col=col,cex=2)
boxplot(beta.cord[top6[3],]*100~pd.cord$td, ylab='Percent Methylation', main=top6[3], cex=2,cex.main=2,cex.lab=1.2,cex.axis=1.5, col=colight2, medcol=colight1, whiskcol=colight1, staplecol=colight1, boxcol=colight1, outcol='white')
stripchart(beta.cord[top6[3],]*100~pd.cord$td, vertical=T, method='jitter',add=T,pch=20,col=col,cex=2)
boxplot(beta.cord[top6[4],]*100~pd.cord$td, ylab='Percent Methylation', main=top6[4], cex=2,cex.main=2,cex.lab=1.2,cex.axis=1.5, col=colight2, medcol=colight1, whiskcol=colight1, staplecol=colight1, boxcol=colight1, outcol='white')
stripchart(beta.cord[top6[4],]*100~pd.cord$td, vertical=T, method='jitter',add=T,pch=20,col=col,cex=2)
boxplot(beta.cord[top6[5],]*100~pd.cord$td, ylab='Percent Methylation', main=top6[5], cex=2,cex.main=2,cex.lab=1.2,cex.axis=1.5, col=colight2, medcol=colight1, whiskcol=colight1, staplecol=colight1, boxcol=colight1, outcol='white')
stripchart(beta.cord[top6[5],]*100~pd.cord$td, vertical=T, method='jitter',add=T,pch=20,col=col,cex=2)
boxplot(beta.cord[top6[6],]*100~pd.cord$td, ylab='Percent Methylation', main=top6[6], cex=2,cex.main=2,cex.lab=1.2,cex.axis=1.5, col=colight2, medcol=colight1, whiskcol=colight1, staplecol=colight1, boxcol=colight1, outcol='white')
stripchart(beta.cord[top6[6],]*100~pd.cord$td, vertical=T, method='jitter',add=T,pch=20,col=col,cex=2)
dev.off()
