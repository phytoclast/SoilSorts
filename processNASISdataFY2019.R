# load required libraries
#devtools::install_github("ncss-tech/soilDB", dependencies=FALSE, upgrade_dependencies=FALSE, force=T)
#remotes::install_github("ncss-tech/soilDB", dependencies=FALSE, upgrade=FALSE, build=FALSE) #20190101
library(soilDB)
library(plyr)
library(foreign)

######################################----

#---- 
#mu <- get_component_correlation_data_from_NASIS_db(SS=F)
#saveRDS(mu, file='fy2020/mu.RDS')
mu <- readRDS(file='fy2020/mu.RDS')

mlra_areas <- read.delim("fy2020/MLRA_LMUKEY.txt")
listoflmuiid <- unique(mu[,c('lmapunitiid', 'nationalmusym','muname')])
listoflmuiid <- merge(listoflmuiid, mlra_areas, by.x = 'lmapunitiid', by.y = 'NE_Soil')

lmuiidagg <- aggregate(listoflmuiid[,'km2'], by=list(listoflmuiid$nationalmusym, listoflmuiid$muname, listoflmuiid$MLRARSYM ), FUN='sum')
colnames(lmuiidagg) <- c('nationalmusym', 'muname', 'MLRARSYM', 'km2')

lmuiidaggtotalmu <- aggregate(lmuiidagg[,'km2'], by=list(lmuiidagg$nationalmusym), FUN='sum')
colnames(lmuiidaggtotalmu) <- c('nationalmusym', 'totalmu')
lmuiidagg <- merge(lmuiidagg, lmuiidaggtotalmu, by='nationalmusym')
lmuiidagg$percent <- lmuiidagg$km2/lmuiidagg$totalmu * 100

lmuiidaggtotalmlra <- aggregate(lmuiidagg[,'km2'], by=list(lmuiidagg$MLRARSYM), FUN='sum')
colnames(lmuiidaggtotalmlra) <- c('MLRARSYM', 'totalmlra')
lmuiidagg <- merge(lmuiidagg, lmuiidaggtotalmlra, by='MLRARSYM')
lmuiidagg$preaffinity <- lmuiidagg$km2/lmuiidagg$totalmlra * 100


lmuiidaggtotalmlra <- aggregate(lmuiidagg[,'preaffinity'], by=list(lmuiidagg$nationalmusym), FUN='sum')
colnames(lmuiidaggtotalmlra) <- c('nationalmusym', 'totalpreaffinity')
lmuiidagg <- merge(lmuiidagg, lmuiidaggtotalmlra, by='nationalmusym')
lmuiidagg$mlra_affinity <- lmuiidagg$preaffinity/lmuiidagg$totalpreaffinity * 100

write.csv(lmuiidagg, file='output/mlraaffinity2020.csv', row.names = F)
write.dbf(listoflmuiid, "output/listoflmuiid.dbf") 

######################################----
#comp <- get_component_data_from_NASIS_db(SS=F)
#saveRDS(comp, file='fy2020/comp.RDS')
comp <- readRDS("fy2020/comp.RDS")

#cm <- get_comonth_from_NASIS_db(fill = TRUE, SS=F)
#saveRDS(cm, file='fy2020/cm.RDS')
cm <- readRDS("fy2020/cm.RDS")

#pm <- get_component_copm_data_from_NASIS_db(SS=F)
#saveRDS(pm, file='fy2020/pm.RDS')
pm <- readRDS(file='fy2020/pm.RDS')



#fc <- fetchNASIS(from='components', fill=TRUE, SS=F, rmHzErrors=FALSE)
#saveRDS(fc, file='fy2020/fc.RDS')
fc <- readRDS(file='fy2020/fc.RDS')


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
OSD.pmorigin$q_param <- tolower(OSD.pmorigin$pmorigin) #change case to match later comparison
OSD.pmkind$q_param <- tolower(OSD.pmkind$pmkind) #change case to match later comparison
OSD.pmkind <- subset(OSD.pmkind, select= -c(pmkind))
OSD.pmorigin <- subset(OSD.pmorigin, select= -c(pmorigin))
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
Appalachia <- c("124", "126", "127", "130A", "136", "139A", "139B", "140", "141", "142", "143", "144A", "144B", "145", "146", "147", "148", "149A", "149B", "153C", "153D")
Plains <- c("101", "108", "110", "111A", "111B", "111C", "111D", "111E", "114A", "90A", "93B","94AA", "94AB", "94B",  "94C",  "95B",  "96A",  "96B",  "97A", "97B", "98A", "98A2", "98B", "99A", "99B")
Northern <- c("94AA", "94AB", "94B",  "94C",  "95B",  "96A",  "96B")
Frigid <- c("94AA", "94AB", "94B",  "94C",  "95B")
LRUMUKEY$LRU <- as.character(LRUMUKEY$LRU)
LRUMUKEY$LRU <- ifelse(grepl('97B',LRUMUKEY$LRU),'97B',LRUMUKEY$LRU) #cleanup the MLRA labels
LRUMUKEY$LRU <- ifelse(grepl('97A',LRUMUKEY$LRU),'97A',LRUMUKEY$LRU)
LRUMUKEY$LRU <- ifelse(grepl('98B',LRUMUKEY$LRU),'98B',LRUMUKEY$LRU)
LRUMUKEY$LRU <- ifelse(grepl('98A',LRUMUKEY$LRU) & !grepl('98A2',LRUMUKEY$LRU),'98A',LRUMUKEY$LRU)
LRUMUKEY$LRU <- ifelse(grepl('111C',LRUMUKEY$LRU),'111C',LRUMUKEY$LRU)
LRUMUKEY$LRU <- ifelse(grepl('111B',LRUMUKEY$LRU),'111B',LRUMUKEY$LRU)
LRUMUKEY$LRU <- ifelse(grepl('110',LRUMUKEY$LRU),'110',LRUMUKEY$LRU)
LRUMUKEY$LRU <- ifelse(grepl('108',LRUMUKEY$LRU),'108',LRUMUKEY$LRU)
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
LRUdmumax$zone <- ifelse(LRUdmumax$LRU %in% Northern, ifelse(LRUdmumax$LRU %in% Frigid,'Frigid','North'),'South')
LRUdmumaxzone <- aggregate(LRUdmumax[,c('affinity')], by=list(LRUdmumax$dmuiid, LRUdmumax$zone), FUN='max')
colnames(LRUdmumaxzone) <- c('dmuiid', 'zone', 'maxzone')
LRUdmumax <- merge(LRUdmumax, LRUdmumaxzone, by=c('dmuiid', 'zone'))
LRUdmumax <- LRUdmumax[LRUdmumax$affinity >=0.2 & LRUdmumax$affinity/LRUdmumax$maxzone > 0.99,] #take out the minority LRUs

preLRUdmumax <- unique(LRUdmumax[,c('LRU', 'dmuiid')]) #find only DMU and LRU combos
LRUduplicates <- aggregate(preLRUdmumax[,'LRU'], by=list(preLRUdmumax$dmuiid), FUN='length') #count duplicates
colnames(LRUduplicates) <- c('dmuiid','dup')
LRUdupout <- merge(LRUdmumax[LRUdmumax$dmumax - LRUdmumax$lrupercent <= 0.05,], LRUduplicates,  by='dmuiid') #output table merging counts with dominant LRU 
LRUdupoutmu <- merge( mu[,c('lmapunitiid', 'dmuiid', 'muname')], LRUdupout, by='dmuiid')

####################### 2018-11-16 - add more to make only one ecosite per zone....

LRUdupoutmu$lmapunitiid <- as.character(LRUdupoutmu$lmapunitiid)
write.dbf(LRUdupoutmu, file='output/LRUdupout.dbf')
rm(LRULLMUKEY, LRUdmumax, LRULMUKEY, LRUduplicates, LRUdupout, MUKEYmax, preLRUdmumax)
rm(mlraranks, mlra_areas1, mlrarank2)
# fix null values in particle size for mucky soils ----
fc.hz$fragvoltot_r <- ifelse(fc.hz$om_r > 25, ifelse(is.na(fc.hz$fragvoltot_r), 0,fc.hz$fragvoltot_r ), fc.hz$fragvoltot_r)
fc.hz$sandtotal_r <- ifelse(fc.hz$om_r > 25, ifelse(is.na(fc.hz$sandtotal_r), 0,fc.hz$sandtotal_r ), fc.hz$sandtotal_r)
fc.hz$silttotal_r <- ifelse(fc.hz$om_r > 25, ifelse(is.na(fc.hz$silttotal_r), 0,fc.hz$silttotal_r ), fc.hz$silttotal_r)
fc.hz$claytotal_r <- ifelse(fc.hz$om_r > 25, ifelse(is.na(fc.hz$claytotal_r), 0,fc.hz$claytotal_r ), fc.hz$claytotal_r)
#top 50 cm ----
fc.hz.notnull <- fc.hz[!is.na(fc.hz$sandtotal_r),]
fc.hz.notnull$h50b <- ifelse(fc.hz.notnull$hzdepb_r > 50, 50, fc.hz.notnull$hzdepb_r)
fc.hz.notnull$h50t <- ifelse(fc.hz.notnull$hzdept_r > 50, 50, fc.hz.notnull$hzdept_r)
fc.hz.notnull$h50k <- fc.hz.notnull$h50b - fc.hz.notnull$h50t
fc.hz.notnull$h50sandk <- fc.hz.notnull$h50k * fc.hz.notnull$sandtotal_r
fc.hz.notnull$pH <- ifelse(!is.na(fc.hz.notnull$ph1to1h2o_r) & fc.hz.notnull$ph1to1h2o_r > 0, fc.hz.notnull$ph1to1h2o_r,
                           ifelse(!is.na(fc.hz.notnull$ph01mcacl2_r) & fc.hz.notnull$ph01mcacl2_r > 0,fc.hz.notnull$ph01mcacl2_r,NA))
fc.hz.notnull$h50phk <- fc.hz.notnull$pH * fc.hz.notnull$h50k

fc.hz.sum <- aggregate(fc.hz.notnull[,c('h50sandk', 'h50phk', 'h50k')], by=list(fc.hz.notnull$coiid), FUN = 'sum' )

colnames(fc.hz.sum) <- c('coiid', 'T50_sand', 'T50_pH', 'h50k')
x <- (as.data.frame(fc.hz.sum[,1]))
colnames(x)<- 'coiid'
fc.hz50.total <- cbind(x, fc.hz.sum[,2:ncol(fc.hz.sum)]/fc.hz.sum$h50k)
#top 150 cm ----
fc.hz.notnull$h150b <- ifelse(fc.hz.notnull$hzdepb_r > 150, 150, fc.hz.notnull$hzdepb_r)
fc.hz.notnull$h150t <- ifelse(fc.hz.notnull$hzdept_r > 150, 150, fc.hz.notnull$hzdept_r)
fc.hz.notnull$h150k <- fc.hz.notnull$h150b - fc.hz.notnull$h150t
fc.hz.notnull$h150sandk <- fc.hz.notnull$h150k * fc.hz.notnull$sandtotal_r
fc.hz.notnull$h150clayk <- fc.hz.notnull$h150k * fc.hz.notnull$claytotal_r
fc.hz.notnull$h150omk <- fc.hz.notnull$h150k * fc.hz.notnull$om_r

fc.hz.sum <- aggregate(fc.hz.notnull[,c('h150sandk', 'h150clayk', 'h150omk', 'h150k')], by=list(fc.hz.notnull$coiid), FUN = 'sum' )

colnames(fc.hz.sum) <- c('coiid','T150_sand', 'T150_clay', 'T150_OM', 'h150k' )
x <- (as.data.frame(fc.hz.sum[,1]))
colnames(x)<- 'coiid'
fc.hz150.total <- cbind(x, fc.hz.sum[,2:ncol(fc.hz.sum)]/fc.hz.sum$h150k)

#depths ----
fc.hz.bedrock <- fc.hz[fc.hz$rock %in% 'BR',]
fc.hz.bedrock.min <- aggregate(fc.hz.bedrock[,c('hzdept_r')], by=list(fc.hz.bedrock$coiid), FUN='min')
colnames(fc.hz.bedrock.min) <- c('coiid','rockdepth')

fc.hz.carb <- fc.hz[fc.hz$caco3_r >0 & !is.na(fc.hz$caco3_r),]
fc.hz.carb.min <- aggregate(fc.hz.carb[,c('hzdept_r')], by=list(fc.hz.carb$coiid), FUN='min')
colnames(fc.hz.carb.min) <- c('coiid','carbdepth')

fc.hz.humic <- fc.hz[!is.na(fc.hz$om_r) & (fc.hz$om_r >= 2.5 | (fc.hz$om_r >= 2.0 & fc.hz$hzdept_r < 1)) ,]
fc.hz.mineral <- fc.hz[fc.hz$om_r < 2 & !is.na(fc.hz$om_r) ,]
fc.hz.humic.max <- aggregate(fc.hz.humic[,c('hzdepb_r')], by=list(fc.hz.humic$coiid), FUN='max')
colnames(fc.hz.humic.max) <- c('coiid','humicmax')
fc.hz.humic.min <- aggregate(fc.hz.mineral[,c('hzdept_r')], by=list(fc.hz.mineral$coiid), FUN='min')
colnames(fc.hz.humic.min) <- c('coiid','humicmin')
fc.hz.humic.max <- merge(fc.hz.humic.max, fc.hz.humic.min, by='coiid', all.x= TRUE)
fc.hz.humic.max$humicmin <- ifelse(is.na(fc.hz.humic.max$humicmin), fc.hz.humic.max$humicmax, fc.hz.humic.max$humicmin )
fc.hz.humic.max$humicdepth <- 0
fc.hz.humic.max$humicdepth <- pmin(fc.hz.humic.max$humicmin, fc.hz.humic.max$humicmax)
fc.hz.humic.max <- subset(fc.hz.humic.max, select = -c(humicmax, humicmin))
#select_s <- fc.site[fc.site$coiid %in% 245287,]
#selecthz <- fc.hz[fc.hz$coiid %in% c(239794 ),]
cm.flood <- cm[!is.na(cm$flodfreqcl) & !(cm$flodfreqcl %in% 'none'),]
cm.flood$flood <- 'flood'
cm.flood <- unique(cm.flood[,c('coiid', 'flood')])
#merge ----
compsorts <- merge(comp[,c('coiid', 'majcompflag','comppct_r' ,'compname', 'taxorder', 'taxsubgrp', 'taxclname', 'drainagecl' ,'slope_r')], fc.hz.carb.min, by='coiid', all.x = T )
compsorts <- merge(compsorts,fc.hz.bedrock.min, by='coiid', all.x = T )
compsorts$carbdepth <- ifelse(is.na(compsorts$carbdepth), 250, compsorts$carbdepth)
compsorts$rockdepth <- ifelse(is.na(compsorts$rockdepth), 250, compsorts$rockdepth)
compsorts <- merge(compsorts,fc.hz150.total[,c('coiid', 'T150_sand', 'T150_clay',   'T150_OM')] , by='coiid', all.x = T )
compsorts <- merge(compsorts,fc.hz50.total[,c('coiid', 'T50_sand', 'T50_pH')] , by='coiid', all.x = T )
compsorts <- merge(compsorts, cm.flood, by='coiid', all.x = T )
compsorts$flood <- ifelse(is.na(compsorts$flood), 'none', 'flood')
compsorts <- merge(compsorts, fc.hz.humic.max, by='coiid', all.x = T )
compsorts$humicdepth <- ifelse(is.na(compsorts$humicdepth), 0, compsorts$humicdepth)

compsorts <- merge(compsorts,fc.site[,c('coiid', 'pmkind', 'pmorigin', 'landform_string')] , by='coiid', all.x = T )
compsorts <- merge(comp[,c('dmuiid','coiid')],compsorts, by='coiid' )
compsorts <- merge(mu[,c('lmapunitiid','dmuiid','nationalmusym','muname')], compsorts, by.x = 'dmuiid', by.y = 'dmuiid', all.y = F )
#fix pmorigin
compsorts <- merge(compsorts, pmlist, by='coiid', all.x=TRUE)
compsorts$pmorigin <- as.character(compsorts$pmorigin)
compsorts$pmorigin <- ifelse(is.na(compsorts$pmorigin) | compsorts$pmorigin %in% 'NA', compsorts$pmorigin2, compsorts$pmorigin)
compsorts <- subset(compsorts, select = -c(pmorigin2))
#add data from comp table
s <- merge(compsorts, comp[,c('coiid','hydricrating')], by.x = 'coiid', by.y = 'coiid', all.y = F )
s <- merge(s, unique(spodic[spodic$Bhs %in% 'yes',c('Bhs', 'compname')]), by='compname', all.x = T)
s$Bhs <- ifelse(is.na(s$Bhs),'no','yes')
s$slope_r <- ifelse(is.na(s$slope_r), 0,s$slope_r)

s$Water_Table <- 250
s$Water_Table <- as.numeric(s$Water_Table)
s$Water_Table <- ifelse(s$hydricrating %in% "no" & !is.na(s$hydricrating), 
                        ifelse(s$drainagecl %in% "moderately well" & !is.na(s$drainagecl), 100,
                               ifelse(grepl('poorly', s$drainagecl)& !is.na(s$drainagecl), 50, s$Water_Table)),
                        ifelse(s$hydricrating %in% "yes" & !is.na(s$hydricrating),
                               ifelse(s$drainagecl %in% "very poorly" & !is.na(s$drainagecl),0,25),
                               ifelse(s$drainagecl %in% "very poorly" & !is.na(s$drainagecl),0,
                                      ifelse(s$drainagecl %in% "poorly" & !is.na(s$drainagecl),25,
                                             ifelse(s$drainagecl %in% "somewhat poorly"& !is.na(s$drainagecl),50,
                                                    ifelse(s$drainagecl %in% "moderately well" & !is.na(s$drainagecl),100,
                                                           250))))))



rm(compsorts, fc.hz.bedrock,  fc.hz.bedrock.min, fc.hz.carb, fc.hz50.total, fc.hz.carb.min, fc.hz.notnull, fc.hz.sum, fc.hz150.total, cm.flood, fc.hz.humic.max, fc.hz.humic.max, fc.hz.humic, fc.hz.mineral )


summarydrcl <- aggregate(s[,c('drainagecl')], by=list(s$drainagecl, s$hydricrating,s$Water_Table), FUN = 'length')

s <- merge(s, unique(LRUdupoutmu[,c('lmapunitiid', 'LRU')]), by.x='lmapunitiid', by.y = 'lmapunitiid', all.x = F)

#0 clear sites______________________________________________________________
s$Site<-"Not"



#1 Lake and beaches______________________________________________________________
s$Site<-ifelse(s$Site %in% "Not",
               ifelse(s$muname %in% c("Water","water") & s$compname %in% c("Water","water"),"Water",
                      ifelse(s$compname %in% c("Beaches","Lake beach","Lake bluffs", "Dune land")| s$muname %in% c("Beach and Dune sand","Stony lake beaches","Dune land","Sand dunes"),"Shoreline Complex","Not"))
               ,s$Site)
#2 Bedrock Floodplains Mucks_____________________________________________________________                
s$Site<-ifelse(s$Site %in% "Not",
               ifelse(s$muname %in% c("flooded") | s$compname %in% c("Alluvial land")|grepl("flood",s$landform_string)|s$flood %in% 'flood' | grepl('flood', s$landform_string),
                      ifelse(s$Water_Table>=50,"F1 Floodplain","F2 Wet Floodplain"),
                      ifelse((s$Water_Table<0 & s$T150_OM >=20)|((grepl("Histic",s$taxsubgrp)|grepl("ists",s$taxsubgrp)|grepl("ists",s$taxclname))& s$Water_Table<=0),"Muck",
                             ifelse(s$rockdepth < 150 & s$Water_Table >100,"Bedrock","Not")))
               ,s$Site)
#3 Sand-Till---------------_____________________________________________________________  
s$Site<-ifelse(s$Site %in% "Not",
               ifelse(
                 (((s$T150_sand >= 80 & s$T50_sand >= 70)|(s$T50_sand >= 80)|(s$T150_clay < 20 & s$T50_pH < 6 & grepl('ult', s$taxsubgrp))) & !(s$LRU %in% Appalachia)), "Sandy",
                 ifelse((s$LRU %in% Appalachia) & 
                   (
                   (s$T50_pH < 5.5 & s$carbdepth > 150 & !grepl('oll',s$taxsubgrp)) |
                     (s$T50_pH < 6 & s$carbdepth > 100 & grepl('frag',s$taxsubgrp)) |
                     grepl('ult', s$taxsubgrp)|grepl('dys', s$taxsubgrp)
                 ),
                 "Acid","Calcareous"))
                 ,s$Site)

#4A Sandy________________________________________________________
s$Site<-ifelse(s$Site %in% "Sandy",
               ifelse(s$Water_Table<=25,"S1 Wet Sandy Depression",
                      ifelse(s$Water_Table<=50,"S2 Sandy Depression",
                             ifelse(s$Water_Table<=100,"S3 Moist Sandy Plains",
                                    ifelse(s$slope_r >=15,"S5 Sandy Slopes","S4 Sandy Plains"))))
               ,s$Site)
#spodic
s$Site<-ifelse(s$Site %in% c("S5 Sandy Slopes") & s$LRU %in% Northern,
               ifelse(grepl('Spod', s$taxorder)|grepl('Spod', s$taxsubgrp)|grepl('spod', s$taxorder)|grepl('spod', s$taxsubgrp),
                      ifelse(s$Bhs %in% 'yes', "SRS5 Rich Spodic Slopes", "SS5 Spodic Slopes"),s$Site)
               ,s$Site)

s$Site<-ifelse(s$Site %in% c("S4 Sandy Plains") & s$LRU %in% Northern,
               ifelse(grepl('Spod', s$taxorder)|grepl('Spod', s$taxsubgrp)|grepl('spod', s$taxorder)|grepl('spod', s$taxsubgrp),
                      ifelse(s$Bhs %in% 'yes', "SRS4 Rich Spodic Plains", "SS4 Spodic Plains"),s$Site)
               ,s$Site)
#rich
s$Site<-ifelse(s$Site %in% c("S5 Sandy Slopes") & (s$LRU %in% Northern |s$LRU %in% c('97A', '97B')),
               ifelse(s$T150_sand < 80 | s$T50_pH >= 6 | s$carbdepth < 100, "SR5 Rich Sandy Slopes", s$Site)
               ,s$Site)

s$Site<-ifelse(s$Site %in% c("S4 Sandy Plains") & s$LRU %in% Northern,
               ifelse(s$T150_sand < 80 | s$T50_pH >= 6 | s$carbdepth < 100, "SR4 Rich Sandy Plains", s$Site)
               ,s$Site)
#acidic
s$Site<-ifelse(s$Site %in% c("S1 Wet Sandy Depression") & !(s$LRU %in% c('99A', '99B')),
               ifelse(grepl('Spod', s$taxorder)|grepl('Spod', s$taxsubgrp)|grepl('spod', s$taxorder)|grepl('spod', s$taxsubgrp)|grepl('ult', s$taxsubgrp)|(!grepl('oll', s$taxsubgrp) & s$T50_pH < 5.5), "SA1 Wet Acidic Sandy Depression",s$Site), s$Site)

s$Site<-ifelse(s$Site %in% c("S2 Sandy Depression") & !(s$LRU %in% c('99A', '99B')),
               ifelse(grepl('Spod', s$taxorder)|grepl('Spod', s$taxsubgrp)|grepl('spod', s$taxorder)|grepl('spod', s$taxsubgrp)|grepl('ult', s$taxsubgrp)|(!grepl('oll', s$taxsubgrp) & s$T50_pH < 5.5), "SA2 Acidic Depression",s$Site), s$Site)

s$Site<-ifelse(s$Site %in% c("S3 Moist Sandy Plains") & !(s$LRU %in% c('99A', '99B')),
               ifelse(grepl('Spod', s$taxorder)|grepl('Spod', s$taxsubgrp)|grepl('spod', s$taxorder)|grepl('spod', s$taxsubgrp)|grepl('ult', s$taxsubgrp)|(!grepl('oll', s$taxsubgrp) & s$T50_pH < 5.5), "SA3 Moist Acidic Plains",s$Site), s$Site)


#4B Tills________________________________________________________
s$Site<-ifelse(s$Site %in% "Calcareous",
               ifelse(s$Water_Table<=25,"C1 Wet Loamy Depression",
                      ifelse(s$Water_Table<=50,"C2 Loamy Depression",
                             ifelse(s$Water_Table<=100,"C3 Moist Loamy Plains",
                                    ifelse(s$slope_r >=15,"C5 Loamy Slopes", "C4 Loamy Plains"))))
               ,s$Site)
#prairie & (s$LRU %in% c('98A', '98A1', '98B', '108',  '110', '111C', '111B'))
s$Site<-ifelse(s$Site %in% c("C2 Loamy Depression") ,
              ifelse(grepl('oll', s$taxorder) | (s$humicdepth >=25 & grepl('ult', s$taxsubgrp )), "CP2 Moist Loamy Prairie",s$Site),s$Site)

s$Site<-ifelse(s$Site %in% c("C3 Moist Loamy Plains")  ,
               ifelse(grepl('oll', s$taxorder) | (s$humicdepth >=25 & grepl('ult', s$taxsubgrp )), "CP3 Moist Loamy Prairie",s$Site),s$Site)

s$Site<-ifelse(s$Site %in% c("C4 Loamy Plains") ,
               ifelse(grepl('oll', s$taxorder) | (s$humicdepth >=25 & grepl('ult', s$taxsubgrp )), "CP4 Loamy Prairie",s$Site),s$Site)


#4C Acid Tills________________________________________________________  
s$Site<-ifelse(s$Site %in% "Acid",
               ifelse(s$Water_Table<=25,"A1 Wet Acid Depression",
                      ifelse(s$Water_Table<=50,"A2 Acid Depression",
                             ifelse(s$Water_Table<=100,"A3 Moist Acid Plains",
                                    ifelse(s$slope_r >=15,"A5 Acid Slopes", "A4 Acid Plains"))))
               ,s$Site)

#5 Outwash________________________________________________________
s$Site<-ifelse(s$Site %in% "Outwash",
               ifelse(s$Water_Table<=25,"O1 Wet Outwash Depression",
                      ifelse(s$Water_Table<=50,"O2 Outwash Depression",
                             ifelse(s$Water_Table<=100,"O3 Moist Outwash",
                                    ifelse(s$slope_r>=15,"O5 Outwash Slopes","O4 Outwash"))))
               ,s$Site)




#6 Bedrock________________________________________________________
s$Site<-ifelse(s$Site %in% "Bedrock",
               ifelse(grepl("sandstone", s$pmorigin)|grepl("sandstone", s$pmkind)|(grepl("shale", s$pmkind) & s$carbdepth > 200 & s$T50_pH < 5.5)|(is.na(s$pmorigin) & s$T50_pH < 6 & s$T150_sand >= 70)|(is.na(s$pmorigin) & s$LRU %in% c('x99B', 'x139A', 'x139B')),
                      ifelse(s$Water_Table<=25,"BS1 Wet Sandstone Depression",
                             ifelse(s$Water_Table<=50,"BS2 Sandstone Depression",
                                    ifelse(s$Water_Table<=100,"BS3 Moist Sandstone",
                                           ifelse(s$slope_r>=15,"BS5 Sandstone Slopes",
                                                  ifelse(s$rockdepth<=100,"BS6 Shallow Sandstone","BS4 Sandstone Residuum"))))),
                      ifelse(grepl("limestone", s$pmorigin)|grepl("dolomite", s$pmorigin)|grepl("limestone", s$pmkind)|grepl("dolomite", s$pmkind)|(grepl("shale", s$pmkind) & s$carbdepth < 150 & s$T50_pH >= 6)| (is.na(s$pmorigin) & grepl('oll', s$taxsubgrp))|(is.na(s$pmorigin) & s$LRU %in% c('x99A', 'x94AB', 'x94C', 'x110')),
                             ifelse(s$Water_Table<=25,"BL1 Wet Limestone Depression",
                                    ifelse(s$Water_Table<=50,"BL2 Limestone Depression",
                                           ifelse(s$Water_Table<=100,"BL3 Moist Limestone",
                                                  ifelse(s$slope_r>=15,"BL5 Limestone Slopes",
                                                         ifelse(s$rockdepth<=100,"BL6 Shallow Limestone","BL4 Limestone Residuum"))))),
                             ifelse(grepl("shale", s$pmorigin)|grepl("unknown", s$pmorigin)|grepl("shale", s$pmkind)|grepl("unknown", s$pmkind),
                                    ifelse(s$Water_Table<=25,"BH1 Wet Shale Depression",
                                           ifelse(s$Water_Table<=50,"BH2 Shale Depression",
                                                  ifelse(s$Water_Table<=100,"BH2 Moist Shale",
                                                         ifelse(s$slope_r>=15,"BH5 Shale Slopes",
                                                                ifelse(s$rockdepth<=100,"BH6 Shallow Shale","BH4 Shale Residuum")))))
                                    ,s$Site)))
               ,s$Site)
#7 Muck________________________________________________________
s$Site<-ifelse(s$Site %in% "Muck",
               ifelse(!is.na(s$T50_pH), ifelse(s$T50_pH < 5.0, "M1 Acidic Peaty Depression", "M2 Mucky Depression"), ifelse(grepl('dysic',s$taxclname), "M1 Acidic Peaty Depression", "M2 Mucky Depression")),s$Site)
 

sitesummary <- aggregate(s[,c('Site')], by = list(s$LRU, s$Site), FUN = 'length')
colnames(sitesummary) <- c('LRU', 'Site', 'Rows')
sitesummary <- sitesummary[sitesummary$LRU %in% c('94AA', '94AB', '94C', '96A', '96B', '97A', '97B', '98A', '98A2', '98B', '99A', '99B', '139A', '139B'),]

write.csv(sitesummary, file= 'output/sitesummary.csv', row.names = F, na = "")          
#Export----


s$lmapunitiid <- as.character(s$lmapunitiid) #ensures that this is a string field, since it joins with gSSURGO mukey (would need to add duplicate row and overwrite contents with a string if exporting as text format instead)
write.dbf(s[!grepl('Urban', s$compname) & (s$majcompflag == 1 | grepl('Urban', s$muname)),], "output/s.dbf")  #exports to dbf excluding minor components unless they are in urban map units, but only showing non urban components.



ss<- s[1:2,] #duplicate a couple of rows
ss$lmapunitiid <- 'A' #ensure that this field will be text (necessary if exporting in text format)
sss <- rbind(ss,s[!grepl('Urban', s$compname) & (s$majcompflag == 1 | grepl('Urban', s$muname)),]) #combine duplicated rows with main data frame
write.table(sss, file='output/s.txt', sep = '\t', row.names = F, na = "") #export ensuring that null data is left blank.

mlraassign <- read.csv("data/mlra139.csv")

mlraassign2 <- merge(s[s$majcompflag == 1 ,], mlraassign, by.x = c('LRU', 'Site'), by.y = c('LRU', 'Site'))
mlraassign2 <- merge(mlraassign2, r12, by = c('dmuiid'))
mlraassign2$lmapunitiid <- as.character(mlraassign2$lmapunitiid) 
write.dbf(mlraassign2, "output/mlraassign2.dbf") 

write.csv(unique(mlraassign2[!is.na(mlraassign2$SiteID) & mlraassign2$SiteID != "", c('coiid','SiteID')]), file='output/mlraassign2.csv', row.names = F)



ecositecooiid <- unique(mlraassign2[!is.na(mlraassign2$SiteID) & mlraassign2$SiteID != "", c('coiid','SiteID')])

dupcoiid <- aggregate(ecositecooiid[,c('SiteID')], by=list(ecositecooiid$coiid), FUN='length')
colnames(dupcoiid) <- c('coiid','Count')
mlraassign3 <- merge(mlraassign2, dupcoiid, by='coiid')
mlraassign3 <- merge(mlraassign3, unique(mu[,c('dmuiid', 'mutype')]), by=c('dmuiid'))
mlraassign3$fix <- 'no'
mlraassign3$fix <- ifelse(mlraassign3$Count > 1 & mlraassign3$mutype %in% 'mlra map unit', 'yes', 'no')
write.dbf(mlraassign3, "output/mlraassign2.dbf") 

write.dbf(mlraassign3[mlraassign3$fix %in% 'yes',], "output/fix.dbf") 

write.csv(unique(mlraassign3[!is.na(mlraassign3$SiteID) & mlraassign3$SiteID != "" & mlraassign3$Count >1 , c('coiid','SiteID', 'nationalmusym','muname' )]), file='output/ecositecoiidduplicates.csv', row.names = F)

write.csv(unique(mlraassign3[!is.na(mlraassign3$SiteID) & mlraassign3$SiteID != "" & mlraassign3$Count >1 , c('nationalmusym','muname', 'mutype')]), file='output/coiidduplicatesmuonly.csv', row.names = F)

