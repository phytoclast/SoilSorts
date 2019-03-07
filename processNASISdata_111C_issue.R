# load required libraries
library(BiodiversityR)
library(soilDB)
library(plyr)
library(foreign)
library(vegan)
library(cluster)
library(ape)
######################################----

#---- 

mlra_areas <- read.delim("data/mlra_areatotal.txt")
mlra_areas$km2 <- floor(mlra_areas$km2 )
mlra_areas1 <- mlra_areas[order(mlra_areas$AREASYMBOL, -mlra_areas$km2),]
mlraranks <- cbind(ddply(mlra_areas1, c("AREASYMBOL"), summarise, rank = rank(-km2)), mlra_areas1[,c(1:3, 5:6)])

mlrarank2 <- merge(unique(mlraranks[mlraranks$MLRA %in% c('94A', '94C', '96', '97', '98', '99', '139'),c('AREASYMBOL', 'LKEY')]),mlraranks[mlraranks$rank == 1,c('AREASYMBOL', 'LRU', 'km2')], by= 'AREASYMBOL')

mlrarank2 <- merge(unique(mlrarank2), mlraranks[mlraranks$rank == 2,c('AREASYMBOL', 'LRU', 'km2')], by= 'AREASYMBOL', all.x = T)
mlrarank2$majority <- ifelse(is.na(mlrarank2$km2.y),1, mlrarank2$km2.x/(mlrarank2$km2.x + mlrarank2$km2.y))
mlrarank2 <- subset(mlrarank2, select = -c(km2.x, km2.y))
colnames(mlrarank2) <- c('AREASYMBOL', 'LKEY', 'PrimLRU', 'SecondLRU', 'majority')

#mlra97 <- mlra_areas[grepl('97', mlra_areas$LRU),]
#unique(mlra97[,4])
#mlra98 <- mlra_areas[grepl('98', mlra_areas$LRU),]
#unique(mlra98[,4])
#mlra96 <- mlra_areas[grepl('96', mlra_areas$LRU),]
#unique(mlra96[,4])
#mlra94 <- mlra_areas[grepl('94', mlra_areas$LRU),]
#unique(mlra94[,4])
#mlra99 <- mlra_areas[grepl('99', mlra_areas$LRU),]
#unique(mlra99[,4])
#mlra111 <- mlra_areas[grepl('111BA', mlra_areas$LRU),]
#unique(mlra111[,4])
#mlra110 <- mlra_areas[grepl('110', mlra_areas$LRU),]
#unique(mlra110[,4])
#mlra139 <- mlra_areas[grepl('139', mlra_areas$LRU),]
#unique(mlra139[,4])
#mlra140 <- mlra_areas[grepl('140', mlra_areas$LRU),]
#unique(mlra140[,4])
#mlra101 <- mlra_areas[grepl('101', mlra_areas$LRU),]
#unique(mlra101[,4])
#mlra127 <- mlra_areas[grepl('127', mlra_areas$LRU),]
#unique(mlra127[,4])
#extra101 <- c('NY063', 'NY029', 'NY121', 'NY003', 'NY073', 'NY037', 'PA105', 'PA083', 'PA047', 'PA065', 'PA031', 'OH019', 'OH081','WV029', 'OH033')
#mlraextra <- mlra_areas[mlra_areas$AREASYMBOL %in% extra101,]
#unique(mlraextra[,4])
######################################----
#comp <- get_component_data_from_NASIS_db(SS=T)
#saveRDS(comp, file='data/comp.RDS')
comp <- readRDS("data/comp.RDS")

#cm <- get_comonth_from_NASIS_db(fill = TRUE, SS=F)
#saveRDS(cm, file='data/cm.RDS')
cm <- readRDS("data/cm.RDS")

#pm <- get_component_copm_data_from_NASIS_db(SS=F)
#saveRDS(pm, file='data/pm.RDS')
pm <- readRDS(file='data/pm.RDS')

#mu <- get_component_correlation_data_from_NASIS_db(SS=T)
#saveRDS(mu, file='data/mu.RDS')
mu <- readRDS(file='data/mu.RDS')

#fc <- fetchNASIS(from='components', fill=TRUE, SS=TRUE, rmHzErrors=FALSE)
#saveRDS(fc, file='data/fc.RDS')
fc <- readRDS(file='data/fc.RDS')


fc.site <- site(fc)
fc.hz <- horizons(fc)
fc.hz$rock <- 'Not'
fc.hz$rock <- ifelse(fc.hz$texture %in% c('UWB','WB','BR')|grepl('Cr', fc.hz$hzname), 'BR','Soil')

#-fix pmorigin

pm <- merge(comp[,c('coiid', 'compname')], pm, by='coiid') # Need component names
pm$pmorigin <- as.character(pm$pmorigin) # factors screw with ifelse replacements, better to convert to character
pm$pmkind <- as.character(pm$pmkind) # factors screw with ifelse replacements, better to convert to character
pm$pmgenmod <- as.character(pm$pmgenmod) # factors screw with ifelse replacements, better to convert to character
pm$pmorigin <- ifelse(is.na(pm$pmorigin), ifelse(grepl('stone',pm$pmkind)|grepl('chert',pm$pmkind)|grepl('igneous',pm$pmkind)|grepl('shale',pm$pmkind)|grepl('quartzite',pm$pmkind),pm$pmkind,pm$pmorigin),pm$pmorigin) #check other field for bedrock
pm$pmorigin <- ifelse(is.na(pm$pmorigin), ifelse(grepl('stone',pm$pmorigin)|grepl('chert',pm$pmgenmod)|grepl('igneous',pm$pmgenmod)|grepl('shale',pm$pmgenmod)|grepl('quartz',pm$pmgenmod),pm$pmgenmod,pm$pmorigin),pm$pmorigin) #check other field for bedrock
fc.hz.rock <-  unique(as.data.frame(fc.hz[fc.hz$rock %in% 'BR',c('coiid')])) #identify all soils containing bedrock
colnames(fc.hz.rock) <- 'coiid'
pmmissing <- pm[is.na(pm$pmorigin),]
pmmissing <- merge(pmmissing, fc.hz.rock, by='coiid') #restrict to coiids that lack data for the pmorigins field

pmmissing2 <- c(unique(pmmissing$compname))
xOSD <- fetchOSD(pmmissing2[1:50], extended = TRUE) #Get OSD pmorigin for only bedrock soils lacking data for this field, only 50 at a time
OSD.pmorigin <- xOSD$pmorigin
OSD.pmkind <- xOSD$pmkind
if(length(pmmissing2)>50){ 
  for(i in 2:floor(length(pmmissing2)/50+1)){ #loop through by 50 and add to list so as not to overwhelm internet query
    xOSD <- fetchOSD(pmmissing2[(i*50-49):(i*50)], extended = TRUE)
    OSD.pmorigin1 <- xOSD$pmorigin
    OSD.pmorigin <- rbind(OSD.pmorigin, OSD.pmorigin1)
    OSD.pmkind1 <- xOSD$pmkind
    OSD.pmkind <- rbind(OSD.pmkind, OSD.pmkind1)
  }}
OSD.pmorigin$q_param <- tolower(OSD.pmorigin$q_param) #change case to match later comparison
OSD.pmkind$q_param <- tolower(OSD.pmkind$q_param) #change case to match later comparison
OSD.pm <- rbind(OSD.pmorigin, OSD.pmkind[grepl('stone', OSD.pmkind$q_param)|grepl('shale', OSD.pmkind$q_param)|grepl('igneous', OSD.pmkind$q_param),])#add in pmkind only if it contains bedrock info
pm$compname <- toupper(pm$compname)
pm <- merge(pm, OSD.pm[,c('series','q_param')], by.x = 'compname', by.y = 'series', all.x = TRUE)
pm$pmorigin <- ifelse(is.na(pm$pmorigin)|(pm$pmorigin %in% 'interbedded sedimentary')|(pm$pmorigin %in% 'mixed'), ifelse(grepl('stone',pm$q_param)|grepl('chert',pm$q_param)|grepl('igneous',pm$q_param)|grepl('shale',pm$q_param)|grepl('quartz',pm$q_param),pm$q_param, pm$pmorigin),pm$pmorigin) #final check of OSD for bedrock

pm <- unique(pm[!is.na(pm$pmorigin),c("coiid", "pmorigin")])

pmlist <- ddply(pm[grepl('stone', pm$pmorigin) | grepl('chert', pm$pmorigin) | grepl('conglomerate', pm$pmorigin) | grepl('igneous', pm$pmorigin)  | grepl('coal', pm$pmorigin) | grepl('shale', pm$pmorigin) | grepl('quartz',pm$pmorigin),], "coiid", function(pmorigin) head(pmorigin,1))
colnames(pmlist) <- c("coiid", "pmorigin2")
selectpm2 <- pmlist[pmlist$coiid %in% 272581,]

fc.hz <- unique(subset(fc.hz, select = -c(structgrpname,texture, hzname)))
#selectC<- fc.site[fc.site$compname %in% 'Capac',]
#fetched.hz <- get_component_horizon_data_from_NASIS_db()
#selecthz <- merge(selectC[,c('compname','coiid')], fc.hz, by='coiid')

#comp$region <- "r12"
#write.table(unique(comp[,c('dmuiid', 'region')]), file='C:/workspace2/r12.txt', sep = '\t', row.names = F, na = "") 
r12 <- read.delim("data/r12.txt")
mur12 <- merge(mu, r12,  by='dmuiid', all.x = T)
mur12$region <- ifelse(mur12$region %in% 'r12', 'r12', 'out')
mur12$lmapunitiid <- as.character(mur12$lmapunitiid)
write.dbf(mur12[,c('lmapunitiid', 'region')], file='output/mur12.dbf')

spodic <- read.delim("data/Spodics.txt")
#### areas <- get_mapunit_from_NASIS() need to get legend areas
#lrunew <- read.csv("C:/workspace2/Export_Output.txt")
#lrunew<- lrunew[,c(3,5,9)]
#colnames(lrunew) <- c('Count', 'MUKEY', 'LRU')
#lrunew$km2 <- lrunew$Count*30*30/10000/100
#write.table(lrunew, file='C:/workspace2/LRUMUKEY2.txt', sep = '\t', row.names = F, na = "" )
#sort(unique(s$LRU))
LRUMUKEY <- read.delim("data/LRUMUKEY2.txt")
#make table of series by LRU----
LRULMUKEY <- merge(LRUMUKEY, mu[,c('dmuiid', 'lmapunitiid')], by.x = 'MUKEY', by.y = 'lmapunitiid')
LRULLMUKEY <- aggregate(LRULMUKEY[,c('km2')], by=list(LRULMUKEY$LRU, LRULMUKEY$dmuiid), FUN='sum')
colnames(LRULLMUKEY) <- c('LRU','dmuiid', 'km2')
LRUtotal <- aggregate(LRULLMUKEY[,c('km2')], by=list(LRULLMUKEY$LRU), FUN='sum')
colnames(LRUtotal) <- c('LRU','total')
LRUtotal$total <- LRUtotal$total + 500 #don't allow spurious pieces too much credit
LRUMUKEYtotal <- merge(LRULLMUKEY, LRUtotal,  by='LRU')
LRUMUKEYtotal$lrupercent <- LRUMUKEYtotal$km2/LRUMUKEYtotal$total
MUKEYmax <- aggregate(LRUMUKEYtotal[,c('lrupercent')], by=list(LRUMUKEYtotal$dmuiid), FUN='max')
colnames(MUKEYmax) <- c('dmuiid','dmumax')

LRUdmumax <- merge(LRUMUKEYtotal, MUKEYmax, by='dmuiid')
LRUdmumax$affinity <- LRUdmumax$lrupercent/LRUdmumax$dmumax

LRUdmumax$zone <- 'Not'

LRUdmumaxzone <- aggregate(LRUdmumax[,c('affinity')], by=list(LRUdmumax$dmuiid, LRUdmumax$zone), FUN='max')
colnames(LRUdmumaxzone) <- c('dmuiid', 'zone', 'maxzone')
LRUdmumax <- merge(LRUdmumax, LRUdmumaxzone, by=c('dmuiid', 'zone'))
preLRUdmumax <- unique(LRUdmumax[,c('LRU', 'dmuiid')]) #find only DMU and LRU combos
LRUduplicates <- aggregate(preLRUdmumax[,'LRU'], by=list(preLRUdmumax$dmuiid), FUN='length') #count duplicates
colnames(LRUduplicates) <- c('dmuiid','dup')
LRUdupout <- merge(LRUdmumax[LRUdmumax$dmumax - LRUdmumax$lrupercent <= 0.05,], LRUduplicates,  by='dmuiid') #output table merging counts with dominant LRU 
LRUdupoutmu <- merge( mu[,c('lmapunitiid', 'dmuiid', 'muname')], LRUdupout, by='dmuiid')
LRUSeries <- merge(LRUdupoutmu[,c('dmuiid','LRU','km2')], fc.site[,c('dmuiid','compname')], by.x = 'dmuiid', by.y = 'dmuiid')
LRUSeriestotal <- aggregate(LRUSeries[,c('km2')], by=list(LRUSeries$LRU), FUN='sum')
colnames(LRUSeriestotal) <- c('LRU','totalkm2')
LRUSeries <- merge(LRUSeries[LRUSeries$LRU !=' '& !is.na(LRUSeries$LRU),], LRUSeriestotal, by='LRU')
LRUSeries$percent <- LRUSeries$km2/LRUSeries$totalkm2*100
LRUSeries$LRU <-  as.character(LRUSeries$LRU)
LRUselect <- LRUSeries[LRUSeries$LRU %in% c("110XY", "111A", "111BA", "111BB", "111CY", "111D", "111E",  "139A", "139B", "94AA", "94AB", "94B",  "94C",  "96A",  "96B",  "97A1", "97A2", "97B1", "97B2", "98A1", "98A2", "98A3", "98A4", "98A5", "98B",  "99A",  "99B"),]
lrumatrix <- makecommunitydataset(LRUselect, row = 'LRU', column = 'compname', value = 'percent')

lrudist <- vegdist(lrumatrix,method='bray', na.rm=T)
lrutree <- agnes(lrudist, method='average')

lrudisttab <- as.data.frame(as.matrix(lrudist))
plot(as.phylo(as.hclust(lrutree)), main='Relationships among LRU based on soil series composition',label.offset=0.125, direction='right', font=1, cex=0.85)

#by subgroup----
LRUSubgrp <- merge(LRUdupoutmu[,c('dmuiid','LRU','km2')], fc.site[!is.na(fc.site$taxsubgrp),c('dmuiid','taxsubgrp')], by.x = 'dmuiid', by.y = 'dmuiid')
LRUSeriestotal <- aggregate(LRUSubgrp[,c('km2')], by=list(LRUSubgrp$LRU), FUN='sum')
colnames(LRUSeriestotal) <- c('LRU','totalkm2')
LRUSubgrp <- merge(LRUSubgrp[LRUSubgrp$LRU !=' '& !is.na(LRUSubgrp$LRU),], LRUSeriestotal, by='LRU')
LRUSubgrp$percent <- LRUSubgrp$km2/LRUSubgrp$totalkm2*100
LRUSubgrp$LRU <-  as.character(LRUSubgrp$LRU)
LRUselect <- LRUSubgrp[LRUSubgrp$LRU %in% c("110XY", "111A", "111BA", "111BB", "111CY", "111D", "111E",  "139A", "139B", "94AA", "94AB", "94B",  "94C",  "96A",  "96B",  "97A1", "97A2", "97B1", "97B2", "98A1", "98A2", "98A3", "98A4", "98A5", "98B",  "99A",  "99B"),]
lrumatrix <- makecommunitydataset(LRUselect, row = 'LRU', column = 'taxsubgrp', value = 'percent')

lrudist <- vegdist(lrumatrix,method='bray', na.rm=T)
lrutree <- agnes(lrudist, method='average')

lrudisttab <- as.data.frame(as.matrix(lrudist))
plot(as.phylo(as.hclust(lrutree)), main='Relationships among LRU based on soil subgroup composition',label.offset=0.125, direction='right', font=1, cex=0.85)


