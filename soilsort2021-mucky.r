library(soilDB)
library(aqp)
library(plyr)
library(foreign)
###2021-12-20


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#Andrew Brown wrote this function ----
library(soilDB)

# works for many different tables (see below for list)
get_group_from_table  <- function(target_table = 'datamapunit'){
  soilDB::dbQueryNASIS(soilDB::NASIS(), sprintf(
    "SELECT * FROM %s 
   INNER JOIN nasisgroup ON %s.grpiidref = nasisgroup.grpiid", target_table, target_table))
}


#establish dominant MRLA per lmukey ----
mlra <- read.csv('fy2021-refresh/mlragrid.csv')
mlra.mu <- read.csv('fy2021-refresh/combgrid.csv')
tgs_zonal <- read.csv('fy2021-refresh/Tgs_zonal2.csv')
mlra <- subset(mlra, Count > 100000)
mu <- readRDS(file='fy2021-refresh/mu.RDS')
mu_tgs <- merge(unique(mu[,c('lmapunitiid', 'dmuiid')]), tgs_zonal, by.x =  'lmapunitiid', by.y = 'VALUE')
mu_tgs$tgswt <- mu_tgs$MEAN*mu_tgs$COUNT
mu_tgs <- aggregate(list(tgswt=mu_tgs$tgswt, wt=mu_tgs$COUNT), by=list(dmuiid =mu_tgs$dmuiid ), FUN='sum')
mu_tgs$tgs <- mu_tgs$tgswt/mu_tgs$wt
mu <-  merge(mu, mu_tgs[,c('dmuiid','tgs')], by='dmuiid', all.x = TRUE)
mlra$MLRA.total <- mlra$Count
mlra$MLRA.value <- mlra$Value

mlra.mu <- merge(mlra[,c('MLRA.value','MLRARSYM', 'MLRA.total')], mlra.mu[,c('Value', 'Count', 'Feature_CONU2', 'northeast_gssurgo30m')], by.x='MLRA.value', by.y='Feature_CONU2')
mlra.mu <- merge(mlra.mu, unique(mu[,c('lmapunitiid', 'dmuiid')]), by.x='northeast_gssurgo30m', by.y='lmapunitiid')
mlra.mu$percentmlra <- mlra.mu$Count/mlra.mu$MLRA.total*100
mlra.mu.agg <- aggregate(list(pmlra = mlra.mu$percentmlra), by=list(dmuiid = mlra.mu$dmuiid, MLRARSYM = mlra.mu$MLRARSYM), FUN='sum')
mlra.mu <- merge(mlra.mu, mlra.mu.agg, by=c('dmuiid', 'MLRARSYM'))
mlra.mu.totals <- aggregate(list(mu.total=mlra.mu.agg$pmlra), by= list(dmu = mlra.mu.agg$dmuiid), FUN ='sum')
mlra.mu <- merge(mlra.mu, mlra.mu.totals, by.x='dmuiid', by.y='dmu')
mlra.mu$affinity <- mlra.mu$pmlra/mlra.mu$mu.total * 100
mlra.mu$MLRA <- mlra.mu$MLRARSYM
mlra.best <- aggregate(list(best=mlra.mu$affinity), by= list(dmuiid = mlra.mu$dmuiid), FUN ='max')
mlra.mu <- merge(mlra.mu, mlra.best, by.x='dmuiid', by.y='dmuiid')
mlra.mu$ofbest <- mlra.mu$affinity / mlra.mu$best *100

mlra.mu <-  unique(subset(mlra.mu, ofbest > 0, select =c(affinity, ofbest, dmuiid, MLRA)))

#establish dominant LRU within MLRA ----
lru.label <- read.csv('fy2021-refresh/lrufeaturelabels.csv')
lru.label <- subset(lru.label, select=c(Value,LRU))

lru.mu <- read.csv('fy2021-refresh/lrusoilcomb.csv')
lru.mu <- subset(lru.mu, select=c(northeast_gssurgo30m,Feature_CONU3,Count))

lru.mu <- merge(lru.mu, lru.label, by.x = 'Feature_CONU3', by.y = 'Value' )
names(lru.mu) <- c('Feature_CONU3', 'lmukey', 'Count', 'LRU'); lru.mu$Feature_CONU3 <- NULL
lru.mu <- merge(lru.mu, mu[,c('lmapunitiid','dmuiid')], by.x='lmukey', by.y='lmapunitiid')
lrulink <- as.data.frame(cbind(
  listofmlra = c('94A',
                 '94A',
                 '94C',
                 '96',
                 '96',
                 '97',
                 '97',
                 '98',
                 '98',
                 '98',
                 '99',
                 '99',
                 '139',
                 '139')
  ,
  listoflru = c('94AA',
                '94AB',
                '94C',
                '96A',
                '96B',
                '97A',
                '97B',
                '98A2',
                '98A1',
                '98B',
                '99A',
                '99B',
                '139A',
                '139B')
))

lru.mu <- merge(lru.mu, lrulink, by.x='LRU', by.y='listoflru')
lru.mu <- aggregate(list(Count=lru.mu$Count), by= list(LRU = lru.mu$LRU, MLRA = lru.mu$listofmlra, dmuiid = lru.mu$dmuiid), FUN ='sum')
lru.totals <- aggregate(list(lrutotal=lru.mu$Count), by= list(LRU = lru.mu$LRU), FUN ='sum')
lru.mu <- merge(lru.mu, lru.totals, by='LRU')
lru.mu$percentlru <- lru.mu$Count/lru.mu$lrutotal*100
lru.mu.total <- aggregate(list(percentlrutotal=lru.mu$percentlru), by= list(dmuiid = lru.mu$dmuiid, MLRA = lru.mu$MLRA), FUN ='sum')
lru.mu <- merge(lru.mu, lru.mu.total, by=c('dmuiid','MLRA'))
lru.mu$lruaffinity <- lru.mu$percentlru/lru.mu$percentlrutotal*100
lru.best <- aggregate(list(best=lru.mu$lruaffinity), by= list(MLRA = lru.mu$MLRA, dmuiid = lru.mu$dmuiid), FUN ='max')
lru.mu <- merge(lru.mu, lru.best, by=c('dmuiid','MLRA'))
lru.mu <- subset(lru.mu, lruaffinity/best >= 0, select=c(lruaffinity,dmuiid, MLRA, LRU))

mlra.mu <- merge(mlra.mu, lru.mu, by=c('dmuiid', 'MLRA'), all.x = T)
######################################----
#usergroup <- get_group_from_table("datamapunit")
#saveRDS(usergroup, 'fy2021-refresh/usergroup20220119.RDS')
usergroup <- readRDS("fy2021-refresh/usergroup20220119.RDS")
# comp <- get_component_data_from_NASIS_db(SS=F)
# saveRDS(comp, file='fy2021-refresh/comp.RDS')
comp <- readRDS("fy2021-refresh/comp.RDS")

# cm <- get_comonth_from_NASIS_db(fill = TRUE, SS=F)
# saveRDS(cm, file='fy2021-refresh/cm.RDS')
cm <- readRDS("fy2021-refresh/cm.RDS")

# pm <- get_component_copm_data_from_NASIS_db(SS=F)
# saveRDS(pm, file='fy2021-refresh/pm.RDS')
pm <- readRDS(file='fy2021-refresh/pm.RDS')

# gm <- get_component_cogeomorph_data_from_NASIS_db(SS=F)
# saveRDS(gm, file='fy2021-refresh/gm.RDS')
gm <- readRDS(file='fy2021-refresh/gm.RDS')

# leg <- get_legend_from_NASIS(SS=F)
# saveRDS(leg, file='fy2021-refresh/leg.RDS')
leg <- readRDS(file='fy2021-refresh/leg.RDS')

# mu <- get_component_correlation_data_from_NASIS_db(SS=F)
# saveRDS(mu, file='fy2021-refresh/mu.RDS')
#mu <- readRDS(file='fy2021-refresh/mu.RDS')

#fc <- fetchNASIS(from='components', fill=TRUE, SS=F, rmHzErrors=FALSE)
#saveRDS(fc, file='fy2021-refresh/fc94A-139.RDS')
#saveRDS(fc, file='fy2021-refresh/fc108A-115A.RDS')
#saveRDS(fc, file='fy2021-refresh/fc101-149B.RDS')
#saveRDS(fc, file='fy2021-refresh/fc120A-153D.RDS')
fc.A <- readRDS(file='fy2021-refresh/fc94A-139.RDS')
fc.B <- readRDS(file='fy2021-refresh/fc108A-115A.RDS')
fc.C <- readRDS(file='fy2021-refresh/fc101-149B.RDS')
fc.D <- readRDS(file='fy2021-refresh/fc120A-153D.RDS')
redundant <- fc.A$coiid
narrow <- subset(fc.B, !coiid %in% redundant)
fc <- aqp::combine(fc.A, narrow)
redundant <- fc$coiid
narrow <- subset(fc.C, !coiid %in% redundant)
fc <- aqp::combine(fc, narrow)
redundant <- fc$coiid
narrow <- subset(fc.D, !coiid %in% redundant)
fc <- aqp::combine(fc, narrow)

#fc <- aqp::rebuildSPC(fc)

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
#Enable if need new OSD data ----
# xOSD <- fetchOSD(pmmissing2[1:50], extended = TRUE) #Get OSD pmorigin for only bedrock soils lacking data for this field, only 50 at a time
# OSD.pmorigin <- xOSD$pmorigin
# OSD.pmkind <- xOSD$pmkind
# 
# if(length(pmmissing2)>50){ 
#   for(i in 2:floor(length(pmmissing2)/50+1)){ #loop through by 50 and add to list so as not to overwhelm internet query
#     xOSD <- fetchOSD(pmmissing2[(i*50-49):(i*50)], extended = TRUE)
#     OSD.pmorigin1 <- xOSD$pmorigin
#     OSD.pmorigin <- rbind(OSD.pmorigin, OSD.pmorigin1)
#     OSD.pmkind1 <- xOSD$pmkind
#     OSD.pmkind <- rbind(OSD.pmkind, OSD.pmkind1)
#   }}
# saveRDS(OSD.pmorigin,'fy2021-refresh/OSD.pmorigin.RDS')
# saveRDS(OSD.pmkind,'fy2021-refresh/OSD.pmkind.RDS')
# ----
OSD.pmorigin <- readRDS('fy2021-refresh/OSD.pmorigin.RDS')
OSD.pmkind <- readRDS('fy2021-refresh/OSD.pmkind.RDS')

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
spodic <- read.delim("fy2021-refresh/Spodics.txt")

# fix null values in particle size for mucky soils ----
fc.hz$fragvoltot_r <- ifelse(is.na(fc.hz$fragvoltot_r), 0,fc.hz$fragvoltot_r )
fc.hz$sandtotal_r <- ifelse(fc.hz$om_r > 25, ifelse(is.na(fc.hz$sandtotal_r), 0,fc.hz$sandtotal_r ), fc.hz$sandtotal_r)
fc.hz$silttotal_r <- ifelse(fc.hz$om_r > 25, ifelse(is.na(fc.hz$silttotal_r), 0,fc.hz$silttotal_r ), fc.hz$silttotal_r)
fc.hz$claytotal_r <- ifelse(fc.hz$om_r > 25, ifelse(is.na(fc.hz$claytotal_r), 0,fc.hz$claytotal_r ), fc.hz$claytotal_r)
#top 50 cm ----
fc.hz.notnull <- fc.hz[!is.na(fc.hz$sandtotal_r),]
fc.hz.notnull$h50b <- ifelse(fc.hz.notnull$hzdepb_r > 50, 50, fc.hz.notnull$hzdepb_r)
fc.hz.notnull$h50t <- ifelse(fc.hz.notnull$hzdept_r > 50, 50, fc.hz.notnull$hzdept_r)
fc.hz.notnull$h50k <- fc.hz.notnull$h50b - fc.hz.notnull$h50t
fc.hz.notnull$h50sandk <- fc.hz.notnull$h50k * fc.hz.notnull$sandtotal_r
fc.hz.notnull$h50clayk <- fc.hz.notnull$h50k * fc.hz.notnull$claytotal_r
fc.hz.notnull$h50siltk <- fc.hz.notnull$h50k * fc.hz.notnull$silttotal_r
fc.hz.notnull$h50fragk <- fc.hz.notnull$h50k * fc.hz.notnull$fragvoltot_r
fc.hz.notnull$h50omk <- fc.hz.notnull$h50k * fc.hz.notnull$om_r
fc.hz.notnull$h50awck <- fc.hz.notnull$h50k * fc.hz.notnull$awc_r
fc.hz.notnull$pH <- ifelse(!is.na(fc.hz.notnull$ph1to1h2o_r) & fc.hz.notnull$ph1to1h2o_r > 0, fc.hz.notnull$ph1to1h2o_r,
                           ifelse(!is.na(fc.hz.notnull$ph01mcacl2_r) & fc.hz.notnull$ph01mcacl2_r > 0,fc.hz.notnull$ph01mcacl2_r,NA))
fc.hz.notnull$h50phk <- fc.hz.notnull$pH * fc.hz.notnull$h50k

fc.hz.sum <- aggregate(fc.hz.notnull[,c('h50fragk','h50sandk','h50siltk','h50clayk', 'h50omk', 'h50awck', 'h50phk', 'h50k')], by=list(fc.hz.notnull$coiid), FUN = 'sum' )

colnames(fc.hz.sum) <- c('coiid', 'T50_frag','T50_sand', 'T50_silt', "T50_clay", 'T50_OM','T50_AWC','T50_pH', 'h50k')
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
fc.hz.notnull$h150siltk <- fc.hz.notnull$h150k * fc.hz.notnull$silttotal_r
fc.hz.notnull$h150fragk <- fc.hz.notnull$h150k * fc.hz.notnull$fragvoltot_r
fc.hz.notnull$h150awck <- fc.hz.notnull$h150k * fc.hz.notnull$awc_r
fc.hz.sum <- aggregate(fc.hz.notnull[,c('h150fragk','h150sandk', 'h150siltk','h150clayk', 'h150omk', 'h150awck', 'h150k')], by=list(fc.hz.notnull$coiid), FUN = 'sum' )

colnames(fc.hz.sum) <- c('coiid','T150_frag', 'T150_sand', 'T150_silt', 'T150_clay', 'T150_OM', 'T150_AWC', 'h150k' )
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
#salt
fc.hz.ec.max <- aggregate(list(ec = fc.hz[!is.na(fc.hz$ec_r),]$ec_r), by=list(coiid=fc.hz[!is.na(fc.hz$ec_r),]$coiid), FUN='max') 
#flood
cm.flood <- cm[!is.na(cm$flodfreqcl) & !(cm$flodfreqcl %in% 'none'),]
cm.flood$flood <- 'flood'
cm.flood <- unique(cm.flood[,c('coiid', 'flood')])
#merge ----
compsorts <- merge(comp[,c('coiid', 'majcompflag','comppct_r' ,'compname', 'taxorder', 'taxsubgrp', 'taxclname', 'drainagecl' ,'slope_r')], fc.hz.carb.min, by='coiid', all.x = T )
compsorts <- merge(compsorts,fc.hz.bedrock.min, by='coiid', all.x = T )
compsorts$carbdepth <- ifelse(is.na(compsorts$carbdepth), 250, compsorts$carbdepth)
compsorts$rockdepth <- ifelse(is.na(compsorts$rockdepth), 250, compsorts$rockdepth)
compsorts <- merge(compsorts,fc.hz50.total[,c('coiid', 'T50_frag','T50_sand','T50_silt','T50_clay','T50_OM', 'T50_AWC', 'T50_pH')] , by='coiid', all.x = T )
compsorts <- merge(compsorts,fc.hz150.total[,c('coiid', 'T150_frag','T150_sand','T150_silt', 'T150_clay','T150_OM', 'T150_AWC')] , by='coiid', all.x = T )
compsorts$T150_AWC <- ifelse(compsorts$rockdepth > 150,150,compsorts$rockdepth)*compsorts$T150_AWC
compsorts$T50_AWC <- ifelse(compsorts$rockdepth > 50,50,compsorts$rockdepth)*compsorts$T50_AWC
compsorts <- merge(compsorts, cm.flood, by='coiid', all.x = T )
compsorts$flood <- ifelse(is.na(compsorts$flood), 'none', 'flood')
compsorts <- merge(compsorts, fc.hz.humic.max, by='coiid', all.x = T )
compsorts$humicdepth <- ifelse(is.na(compsorts$humicdepth), 0, compsorts$humicdepth)
compsorts <- merge(compsorts, fc.hz.ec.max, by='coiid', all.x = T )
compsorts <- merge(compsorts,fc.site[,c('coiid', 'pmkind', 'pmorigin', 'landform_string')] , by='coiid', all.x = T )
compsorts <- merge(comp[,c('dmuiid','coiid')],compsorts, by='coiid' )
compsorts <- merge(mu[,c('areasymbol','areaname','lmapunitiid','dmuiid','muiid','nationalmusym','muname','tgs')], compsorts, by.x = 'dmuiid', by.y = 'dmuiid', all.y = F )
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
s$Water_Table <- ifelse(s$hydricrating %in% "no", 
                        ifelse(s$drainagecl %in% "moderately well", 100,
                               ifelse(s$drainagecl %in% c("somewhat poorly"), 50,
                                      ifelse(s$drainagecl %in% c("poorly","very poorly","subaqueous"), 25, s$Water_Table))), #for urban soils that lost their hydric status, we will continue to recognize them as former wet ES.
                               ifelse(s$hydricrating %in% "yes",
                                      ifelse(s$drainagecl %in% c("very poorly","subaqueous"),0,25),
                                      ifelse(s$drainagecl %in% c("very poorly","subaqueous"),0,
                                             ifelse(s$drainagecl %in% "poorly",25,
                                                    ifelse(s$drainagecl %in% "somewhat poorly",50,
                                                           ifelse(s$drainagecl %in% "moderately well",100,
                                                                  250))))))
                        


rm(compsorts, fc.hz.bedrock,  fc.hz.bedrock.min, fc.hz.carb, fc.hz50.total, fc.hz.carb.min, fc.hz.notnull, fc.hz.sum, fc.hz150.total, cm.flood, fc.hz.humic.max, fc.hz.humic, fc.hz.mineral )


summarydrcl <- aggregate(s[,c('drainagecl')], by=list(s$drainagecl, s$hydricrating,s$Water_Table), FUN = 'length')

s <- merge(mlra.mu, s, by.x = 'dmuiid', by.y = 'dmuiid')




s.reserve <- s

s <- s.reserve
#s <- subset(s, MLRA %in% '98')



#soil sort for MLRA 139 ----
#0 clear sites______________________________________________________________
s$Site<-"Not"
#___________________________________                
s$Site<-ifelse(s$muname %in% c("flooded") | s$compname %in% c("Alluvial land")|
                        (grepl("flood",s$landform_string)& (grepl("fluv",s$taxsubgrp)|grepl("psam",s$taxsubgrp)|grepl("cumu",s$taxsubgrp)))|
                        s$flood %in% 'flood', 'Floodplain',
                      ifelse((s$Water_Table<0 & s$T150_OM >=20)|((grepl("^histic",s$taxsubgrp)|grepl("Histosols",s$compname)|grepl("ists",s$taxsubgrp)|grepl("ists",s$taxclname))& s$Water_Table<=0),"Muck",
                             ifelse(s$rockdepth < 150 & s$Water_Table >100,"Bedrock",s$Site)))
               

#7 Muck________________________________________________________
s$Site<-ifelse(s$Site %in% "Muck",
               ifelse(!is.na(s$T50_pH) & s$T50_pH < 5.0, "Acidic Peaty Depression",ifelse(grepl('dysic',s$taxclname), "Acidic Peaty Depression", "Mucky Depression"))
               ,s$Site)

mlradmus <- unique(s[s$MLRA %in% c('98','97','99') & s$Site %in% 'Mucky Depression' & s$majcompflag %in% 'TRUE',]$dmuiid)
s1 <- subset(s, dmuiid %in% mlradmus)
usergrouplmuid <- merge(s1, usergroup[,c('grpname','dmuiid')], by='dmuiid')
write.csv(usergrouplmuid, 'fy2021-refresh/mlramucks.csv', row.names = F)

