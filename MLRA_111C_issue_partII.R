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
s <- read.dbf("testmlra/s.dbf")

t07$t07 <- t07$MEAN
t01$t01 <- t01$MEAN
MAP$MAP <- MAP$MEAN
MAAT$MAAT <- MAAT$MEAN

MLRAclim <- merge(t07[,c('Unit','t07')], t01[,c('Unit','t01')], by='Unit')
MLRAclim <- merge(MLRAclim, MAAT[,c('Unit','MAAT')], by='Unit')
MLRAclim <- merge(MLRAclim, MAP[,c('Unit','MAP')], by='Unit')

MLRAsoil <- merge(mlra_testunits[,c('Unit','Value')], mlra_testcomb[,c('MLRA_TESTU','MapunitRas','Count')], by.x='Value', by.y='MLRA_TESTU')
MLRAsoil <- merge(MLRAsoil, s, by.x = 'MapunitRas', by.y = 'lmapunitii')

#by subgroup----
SubgrpTotal <- aggregate(MLRAsoil[,c('Count')], by=list(MLRAsoil$Unit), FUN='sum')
colnames(SubgrpTotal) <- c('Unit','totalCount')
MLRAsoil <- merge(MLRAsoil, SubgrpTotal, by='Unit')
MLRAsoil$percent <- MLRAsoil$Count/MLRAsoil$totalCount*100
MLRAsoil <- subset(MLRAsoil, !is.na(percent) & !is.na(taxsubgrp) & !is.na(Unit))
MLRAsoil$Unit1 <-  as.character(paste('M',MLRAsoil$Unit))
lrumatrix <- makecommunitydataset(MLRAsoil, row = 'Unit1', column = 'taxsubgrp', value = 'percent')

lrudist <- vegdist(lrumatrix,method='bray', na.rm=T)
lrutree <- agnes(lrudist, method='average')

lrudisttab <- as.data.frame(as.matrix(lrudist))
plot(as.phylo(as.hclust(lrutree)), main='Relationships among MLRA Units based on soil subgroup composition',label.offset=0.125, direction='right', font=1, cex=0.85)

lrutree2 <- nj(lrudist)
plot(as.phylo((lrutree2)), main='Relationships among MLRA Units based on soil subgroup composition',label.offset=0.125, direction='right', font=1, cex=0.85)
rownames(MLRAclim) <- MLRAclim$Unit
gowerdist <- vegdist(MLRAclim[,2:5], method='gower')

gwoertree <- agnes(gowerdist, method='average')
plot(as.phylo(as.hclust(gwoertree)), main='Relationships among MLRA Units based on soil subgroup composition',label.offset=0.125, direction='right', font=1, cex=0.85)


