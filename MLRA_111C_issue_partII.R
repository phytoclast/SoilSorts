# load required libraries
library(BiodiversityR)
library(soilDB)
library(plyr)
library(foreign)
library(vegan)
library(cluster)
library(ape)
######################################----

mlra_testunits <- read.dbf("testmlra/MLRA_TESTUNITS.dbf")
mlra_testcomb <- read.dbf("testmlra/MLRATESTCOMB.dbf")

t07 <- read.dbf("testmlra/ZonalTestt07.dbf")
t01 <- read.dbf("testmlra/ZonalTestt01.dbf")
MAP <- read.dbf("testmlra/ZonalTestPrecip.dbf")
MAAT <- read.dbf("testmlra/ZonalTestMAAT.dbf")
winterP <- read.dbf("testmlra/ZonalTestwinterP.dbf")
s <- read.dbf("testmlra/s.dbf")
s <- unique(subset(s, select= -c(LRU, Site)))
t07$t07 <- t07$MEAN
t01$t01 <- t01$MEAN
MAP$MAP <- MAP$MEAN
MAP$MAPMIN <- MAP$MIN
MAP$MAPMAX <- MAP$MAX
MAAT$MAAT <- MAAT$MEAN 
MAAT$MAATMAX <- MAAT$MAX 
MAAT$MAATMIN <- MAAT$MIN
winterP$winterP <- winterP$MEAN

MLRAclim <- merge(t07[,c('Unit','t07')], t01[,c('Unit','t01')], by='Unit')
MLRAclim <- merge(MLRAclim, MAAT[,c('Unit','MAAT', 'MAATMAX','MAATMIN')], by='Unit')
MLRAclim <- merge(MLRAclim, MAP[,c('Unit','MAP', 'MAPMAX','MAPMIN')], by='Unit')
MLRAclim <- merge(MLRAclim, winterP[,c('Unit','winterP','COUNT')], by='Unit')


MLRAsoil <- merge(mlra_testunits[,c('Unit','Value')], mlra_testcomb[,c('MLRA_TESTU','MapunitRas','Count')], by.x='Value', by.y='MLRA_TESTU')
MLRAsoil <- merge(MLRAsoil, s, by.x = 'MapunitRas', by.y = 'lmapunitii')

#by subgroup----
SubgrpTotal <- aggregate(MLRAsoil[,c('Count')], by=list(MLRAsoil$Unit), FUN='sum')
colnames(SubgrpTotal) <- c('Unit','totalCount')
MLRAsoil2 <- merge(MLRAsoil, SubgrpTotal, by='Unit')
MLRAsoil2$percent <- MLRAsoil2$Count/MLRAsoil2$totalCount*100
MLRAsoil2 <- subset(MLRAsoil2, !is.na(percent) & !is.na(taxsubgrp) & !is.na(Unit))
MLRAsoil2$Unit1 <-  as.character(paste0('M',MLRAsoil2$Unit))
lrumatrix <- makecommunitydataset(MLRAsoil2, row = 'Unit1', column = 'taxsubgrp', value = 'percent')
#----count soils
MLRAsoil2$mlra <- strsplit(as.character(MLRAsoil2$Unit), split='-')
MLRAsoil2$mlra2 <- as.character(lapply(MLRAsoil2$mlra, "[[",1))

selected <- subset(MLRAsoil2, mlra2 %in% c('94A','94C','96','97','98','99','139'))
MLRAsoilsubgroups <- aggregate(selected$Count, by=list(selected$mlra2, selected$taxsubgrp), FUN='sum')
colnames(MLRAsoilsubgroups) <- c('MLRA','subgroup','percent')
MLRAsoilcompname <- aggregate(selected$Count, by=list(selected$mlra2, selected$compname), FUN='sum')
colnames(MLRAsoilcompname) <- c('MLRA','compname','percent')
write.csv(MLRAsoilsubgroups, 'MLRAsoilsubgroups.csv',row.names = F)
write.csv(MLRAsoilcompname, 'MLRAsoilcompname.csv',row.names = F)
#----

lrudist <- vegdist(lrumatrix,method='bray', na.rm=T)
lrutree <- agnes(lrudist, method='average')

plot(as.phylo(as.hclust(lrutree)), main='Units by subgroup',label.offset=0.125, direction='right', font=1, cex=0.85)

lrutree2 <- nj(lrudist)
plot(as.phylo((lrutree2)), main='Units by subgroup',label.offset=0.125, direction='right', font=1, cex=0.85)
#distmatrix
matrixtable <- as.data.frame(as.matrix(lrudist))

matrixtable <- matrixtable[c('M98.1', 'M98.2', 'M111C', 'M98.111B','M98.111C','M111B.1','M111B.2'),
                           c('M98.1', 'M98.2', 'M111C', 'M98.111B','M98.111C','M111B.1','M111B.2')]


#regroup 
MLRAsoil$Unit <- as.character(MLRAsoil$Unit)
MLRAsoil$Unit2 <- MLRAsoil$Unit
MLRAsoil$Unit2 <- ifelse(MLRAsoil$Unit %in% c('98-1', '98-2', '98-3'), '98', MLRAsoil$Unit2)
MLRAsoil$Unit2 <- ifelse(MLRAsoil$Unit %in% c('98-111C', '98-111B'), '98-111', MLRAsoil$Unit2)
MLRAsoil$Unit2 <- ifelse(MLRAsoil$Unit %in% c('97','97-98', '97-110'), '97', MLRAsoil$Unit2)
MLRAsoil$Unit2 <- ifelse(MLRAsoil$Unit %in% c('96','96-98','96-94A'), '96', MLRAsoil$Unit2)

SubgrpTotal <- aggregate(MLRAsoil[,c('Count')], by=list(MLRAsoil$Unit2), FUN='sum')
colnames(SubgrpTotal) <- c('Unit2','totalCount')
MLRAsoil2 <- merge(MLRAsoil, SubgrpTotal, by='Unit2')
MLRAsoil2 <- subset(MLRAsoil2, !Unit2 %in% c('110-98','98-110-1','98-111C-2'))
MLRAsoil2$percent <- MLRAsoil2$Count/MLRAsoil2$totalCount*100
MLRAsoil2 <- subset(MLRAsoil2, !is.na(percent) & !is.na(taxsubgrp) & !is.na(Unit2))
MLRAsoil2$Unit1 <-  as.character(paste0('M',MLRAsoil2$Unit2))
lrumatrix <- makecommunitydataset(MLRAsoil2, row = 'Unit1', column = 'taxsubgrp', value = 'percent')

lrudist <- vegdist(lrumatrix,method='bray', na.rm=T)
lrutree <- agnes(lrudist, method='average')

plot(as.phylo(as.hclust(lrutree)), main='Units by subgroup',label.offset=0.125, direction='right', font=1, cex=0.85)

#by series----
SubgrpTotal <- aggregate(MLRAsoil[,c('Count')], by=list(MLRAsoil$Unit), FUN='sum')
colnames(SubgrpTotal) <- c('Unit','totalCount')
MLRAsoil2 <- merge(MLRAsoil, SubgrpTotal, by='Unit')
MLRAsoil2$percent <- MLRAsoil2$Count/MLRAsoil2$totalCount*100
MLRAsoil2$Unit1 <-  as.character(paste0('M',MLRAsoil2$Unit))
MLRAsoil2 <- subset(MLRAsoil2, !is.na(percent) & !is.na(compname) & !is.na(Unit1))
lrumatrix <- makecommunitydataset(MLRAsoil2, row = 'Unit1', column = 'compname', value = 'percent')

lrudist <- vegdist(lrumatrix,method='bray', na.rm=T)
lrutree <- agnes(lrudist, method='average')

plot(as.phylo(as.hclust(lrutree)), main='Units by compname',label.offset=0.125, direction='right', font=1, cex=0.85)
#distmatrix
matrixtable <- as.data.frame(as.matrix(lrudist))

matrixtable <- matrixtable[c('M98.1', 'M98.2', 'M111C', 'M98.111B','M98.111C','M111B.1','M111B.2'),
                           c('M98.1', 'M98.2', 'M111C', 'M98.111B','M98.111C','M111B.1','M111B.2')]
#top series
MLRAsoil3<- aggregate(MLRAsoil2[,c('percent')], by=list(MLRAsoil2$Unit, MLRAsoil2$compname), FUN='sum')
colnames(MLRAsoil3) <- c('Unit','compname','percent')
MLRAsoil3 <- MLRAsoil3[with(MLRAsoil3, order(Unit, -percent)),]
MLRAsoil3$rownumber <-  1:nrow(MLRAsoil3)
MLRArnk <- aggregate(MLRAsoil3$rownumber, by=list(MLRAsoil3$Unit), FUN='min')
colnames(MLRArnk) <- c('Unit', 'rowmin')
MLRAsoil3 <- merge(MLRAsoil3, MLRArnk, by='Unit')
MLRAsoil3$rank <- MLRAsoil3$rownumber-MLRAsoil3$rowmin+1
topten <- MLRAsoil3[MLRAsoil3$rank <= 10,]
topten$Unit <- as.character(topten$Unit)
topten$compname <- as.character(topten$compname)
units <- unique(topten$Unit)

topten2 <- topten[topten$Unit %in% units[1] &  topten$rank %in% c(1:10),c('compname')]
for (n in 2:length(unique(topten$Unit))){
toptenn <- topten[topten$Unit %in% units[n] &  topten$rank %in% c(1:10),c('compname')]
topten2 <- cbind(topten2,  toptenn)
}

colnames(topten2) <- paste0('M',units)
write.csv(topten2,'testmlra/toptencompbyMLRA.csv')
topten2a <- topten2[,c('M98-1', 'M98-2', 'M111C', 'M98-111B','M98-111C','M111B-1','M111B-2')]
#clim
rownames(MLRAclim) <- paste0('M',MLRAclim$Unit)
gowerdist <- vegdist(MLRAclim[,2:ncol(MLRAclim)], method='gower')

gwoertree <- agnes(gowerdist, method='average')
plot(as.phylo(as.hclust(gwoertree)), main='Units by climate',label.offset=0.125, direction='right', font=1, cex=0.85)
matrixtable <- as.data.frame(as.matrix(gowerdist))
combdist <- gowerdist+lrudist
combtree <- agnes(combdist, method='average')
plot(as.phylo(as.hclust(combtree)), main='Units by subgroup and climate',label.offset=0.125, direction='right', font=1, cex=0.85)

#clim regroup

MLRAclim$Unit <- as.character(MLRAclim$Unit)
MLRAclim$Unit2 <- MLRAclim$Unit
MLRAclim$Unit2 <- ifelse(MLRAclim$Unit %in% c('98-1', '98-2', '98-3'), '98', MLRAclim$Unit2)
MLRAclim$Unit2 <- ifelse(MLRAclim$Unit %in% c('98-111C', '98-111B'), '98-111', MLRAclim$Unit2)
MLRAclim$Unit2 <- ifelse(MLRAclim$Unit %in% c('97','97-98', '97-110'), '97', MLRAclim$Unit2)
MLRAclim$Unit2 <- ifelse(MLRAclim$Unit %in% c('96','96-98','96-94A'), '96', MLRAclim$Unit2)
MLRAclim2 <- aggregate(MLRAclim[c('t07','t01','MAAT','MAATMAX','MAATMIN','MAP','MAPMAX','MAPMIN','winterP')], by=list(MLRAclim$Unit2), FUN = 'mean', weight=MLRAclim$COUNT)
colnames(MLRAclim2) <- c('Unit2','t07','t01','MAAT','MAATMAX','MAATMIN','MAP','MAPMAX','MAPMIN','winterP')
rownames(MLRAclim2) <- paste0('M',MLRAclim2$Unit2)
MLRAclim2 <- subset(MLRAclim2, !Unit2 %in% c('110-98','98-110-1','98-111C-2'))
gowerdist <- vegdist(MLRAclim2[,2:ncol(MLRAclim2)], method='gower')

gwoertree <- agnes(gowerdist, method='average')
plot(as.phylo(as.hclust(gwoertree)), main='Units by climate',label.offset=0.125, direction='right', font=1, cex=0.85)
matrixtable <- as.data.frame(as.matrix(gowerdist))
combdist <- gowerdist+lrudist
combtree <- agnes(combdist, method='average')
plot(as.phylo(as.hclust(combtree)), main='Units by subgroup and climate',label.offset=0.125, direction='right', font=1, cex=0.85)

