library(soilDB)
library(aqp)
library(plyr)
library(foreign)
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


#load MLRA ownership 
mlraowner <- read.csv('data/s.groupmlra.1.csv')
mlraowner <- subset(mlraowner, officematch %in% 1)
#usergroup <- get_group_from_table("datamapunit")
#saveRDS(usergroup, 'fy2021-refresh/usergroup20220408.RDS')
usergroup <- readRDS("fy2021-refresh/usergroup20220408.RDS")
mlra.mu <- merge(mlra.mu, usergroup[,c('grpname','dmuiid')], by='dmuiid')
mlra.mu <- merge(mlra.mu, mlraowner[,c('grpname','MLRARSYM')], by=c('MLRARSYM','grpname'))


#calculate total
mlra.mu.agg <- aggregate(list(pmlra = mlra.mu$percentmlra), by=list(dmuiid = mlra.mu$dmuiid, MLRARSYM = mlra.mu$MLRARSYM), FUN='sum')
mlra.mu <- merge(mlra.mu, mlra.mu.agg, by=c('dmuiid', 'MLRARSYM'))
mlra.mu.totals <- aggregate(list(mu.total=mlra.mu.agg$pmlra), by= list(dmu = mlra.mu.agg$dmuiid), FUN ='sum')
mlra.mu <- merge(mlra.mu, mlra.mu.totals, by.x='dmuiid', by.y='dmu')
mlra.mu$affinity <- mlra.mu$pmlra/mlra.mu$mu.total * 100
mlra.mu$MLRA <- mlra.mu$MLRARSYM
mlra.best <- aggregate(list(best=mlra.mu$affinity), by= list(dmuiid = mlra.mu$dmuiid), FUN ='max')
mlra.mu <- merge(mlra.mu, mlra.best, by.x='dmuiid', by.y='dmuiid')
mlra.mu$ofbest <- mlra.mu$affinity / mlra.mu$best *100

mlra.mu <-  unique(subset(mlra.mu, ofbest > 99, select =c(dmuiid, MLRA)))

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
lru.mu <- subset(lru.mu, lruaffinity/best >= 0.99, select=c(dmuiid, MLRA, LRU))

mlra.mu <- merge(mlra.mu, lru.mu, by=c('dmuiid', 'MLRA'), all.x = T)
######################################----
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
s <- merge(s, usergroup[,c('grpname','dmuiid')], by='dmuiid')



#### MLRA 111 Parent Material key ----

#0 clear sites --------------------------------------------------------- 
s$Site111 <-"Not"



#1 Lake and beaches --------------------------------------------------------- 
s$Site111<-ifelse(s$Site111 %in% "Not",
                  ifelse(s$muname %in% c("Water","water") & s$compname %in% c("Water","water"),"Water",
                         ifelse(s$compname %in% c("Beaches","Lake beach","Lake bluffs", "Dune land")| s$muname %in% c("Beach and Dune sand","Stony lake beaches","Dune land","Sand dunes"),"Shoreline Complex","Not"))
                  ,s$Site111)
#2  Mucks ---------------------------------------------------------               
s$Site111 <- ifelse(s$Site111 %in% "Not", ifelse(grepl("Histosols",s$compname)|grepl("ists",s$taxsubgrp)|grepl("ists",s$taxclname)|grepl("organic",s$pmkind),"Muck",s$Site111),s$Site111)

s$Site111 <- ifelse(s$Site111 %in% "Muck", ifelse(grepl("Terric",s$compname)|grepl("Terric",s$taxsubgrp)|grepl("Terric",s$taxclname),"Mineral Muck",s$Site111),s$Site111)
s$Site111 <- ifelse(s$Site111 %in% "Muck", ifelse(grepl("Limnic",s$compname)|grepl("Limnic",s$taxsubgrp)|grepl("Limnic",s$taxclname),"Limnic Muck", "Deep Muck"),s$Site111)


#3 lacustrine Parent Material  --------------------------------------------------------- 
s$Site111<-ifelse(s$Site111 %in% "Not",
                  ifelse(grepl("lacustrine",s$pmkind)|grepl("lake plain",s$pmkind), 'lacustrine', s$Site111), s$Site111)

s$Site111<-ifelse(s$Site111 %in% "lacustrine",
                  ifelse(s$drainagecl %in% c("very poorly", "poorly"), 'Lacustrine Flatwood', 'Lacustrine Forest'), s$Site111)

#3 alluvium Parent Material  --------------------------------------------------------- 
s$Site111<-ifelse(s$Site111 %in% "Not",
                  ifelse(grepl("alluvium",s$pmkind), 'alluvium', s$Site111), s$Site111)

s$Site111<-ifelse(s$Site111 %in% "alluvium",
                  ifelse(grepl("Mollisols",s$compname)|grepl("olls",s$taxsubgrp)|grepl("olls",s$taxclname), 
                         ifelse((grepl("aquic",s$taxsubgrp) & !grepl("oxy",s$taxsubgrp))|(grepl("Aquic",s$taxclname) & !grepl("Oxy",s$taxclname)),'Wet Alluvium Floodplain', 'Dry Alluvium Floodplain'),
                         ifelse(s$drainagecl %in% c("very poorly", "poorly"), 'Wet Alluvium Forest', 'Dry Alluvium Forest'))
                  , s$Site111)

#4 bedrock Parent Material  --------------------------------------------------------- 
s$Site111<-ifelse(s$Site111 %in% "Not",
                  ifelse(s$rockdepth < 200, 'bedrock', s$Site111), s$Site111)

s$Site111<-ifelse(s$Site111 %in% "bedrock",
                  ifelse(grepl("Mollisols",s$compname)|grepl("olls",s$taxsubgrp)|grepl("olls",s$taxclname), 'Dark Bedrock Prairie',
                         ifelse(s$drainagecl %in% c("very poorly", "poorly", "somewhat poorly"), 'Mesic Bedrock Forest', 'Dry Bedrock Forest'))
                  , s$Site111)

#5 outwash Parent Material  --------------------------------------------------------- 
s$Site111<-ifelse(s$Site111 %in% "Not",
                  ifelse(grepl("outwash",s$pmkind)|grepl("fluv",s$pmkind), 'outwash', s$Site111), s$Site111)

s$Site111<-ifelse(s$Site111 %in% "outwash",
                  ifelse(grepl("Mollisols",s$compname)|grepl("olls",s$taxsubgrp)|grepl("olls",s$taxclname), 
                         ifelse((grepl("aquic",s$taxsubgrp) & !grepl("oxy",s$taxsubgrp))|
                                  (grepl("Aquic",s$taxclname) & !grepl("Oxy",s$taxclname))| s$drainagecl %in% c("very poorly", "poorly"),
                                'Wet Outwash Mollisol', 'Dry Outwash Intergrade'),
                         ifelse(s$drainagecl %in% c("very poorly", "poorly", "somewhat poorly"), 'Outwash Upland', 'Dry Outwash Upland'))
                  , s$Site111)

#5 till Parent Material  --------------------------------------------------------- 
s$Site111<-ifelse(s$Site111 %in% "Not",
                  ifelse(grepl("till",s$pmkind)|grepl("drift",s$pmkind), 'till', s$Site111), s$Site111)

s$Site111<-ifelse(s$Site111 %in% "till",
                  ifelse(s$drainagecl %in% c("very poorly", "poorly"),'Till Depression',
                         ifelse(s$drainagecl %in% c("very poorly", "poorly"),'Wet Till Ridge','Till Ridge'))
                  , s$Site111)

#111A ----

ES111A.list <- list(
  F111AY001IN = c('Palms', 'Linwood', 'Adrian'),
  R111AY002IN = c('Warners', 'Wallkill', 'Muskego', 'Martisco', 'Edwards', 'Benadum'),
  F111AY003IN = c('Houghton', 'Carlisle'),
  F111AY004IN = c('Aetna', 'Algiers', 'Banlic', 'Bartle', 'Beaucoup', 'Bellcreek', 'Ceresco', 'Cohoctah', 'Euclid', 'Evansville', 'Henshaw', 'Holton', 'Orrville', 'Piopolis', 'Rockmill', 'Saranac', 'Shoals', 'Sloan', 'Southwest', 'Stendal', 'Vincennes', 'Wakeland', 'Washtenaw', 'Wilhite', 'Henshaw Variant', 'Shoals Variant', 'Sloan Variant'), # add 'Henshaw Variant', 'Shoals Variant', 'Sloan Variant'
  F111AY005IN = c('Wirt', 'Wilbur', 'Uniontown', 'Tremont', 'Stonelick', 'Steff', 'Skidmore', 'Sardinia', 'Rossburg', 'Ross', 'Romeo', 'Pekin', 'Otwell', 'Oldenburg', 'Moundhaven', 'Millstone', 'Medway', 'Lobdell', 'Lash', 'Lanier', 'Landes', 'Kinn', 'Haymond', 'Gessie', 'Genesee', 'Elkinsville', 'Eel', 'Dearborn', 'Cuba', 'Clifty', 'Chagrin', 'Beckville', 'Armiesburg', 'Abscota Variant', 'Medway Variant', 'Riverwash', 'Ross Variant'),# add 'Abscota Variant', 'Medway Variant', 'Riverwash', 'Ross Variant'
  F111AY006IN = c('Wetzel', 'Reesville', 'Nappanee', 'Haskins', 'Condit', 'Blount'),
  F111AY007IN = c('Treaty', 'Sidell', 'Pewamo', 'Millsdale', 'Marengo', 'Kokomo', 'Cyclone', 'Brookston'),
  F111AY008IN = c('Westboro', 'Vigo', 'Sugarvalley', 'Schaffer', 'Randolph', 'Pyrmont', 'Fincastle', 'Crosby', 'Cobbsfork', 'Bennington', 'Avonburg', 'Randolph variant'), #add 'Randolph variant'
  F111AY009IN = c('Xenia', 'Wynn', 'Williamstown', 'Wapahani', 'Tuscola', 'Thrifton', 'Tarlton', 'Strawn', 'Senachwine', 'Ryker', 'Russell', 'Rossmoyne', 'Rockfield', 'Ritchey', 'Rawson', 'Rainsville', 'Nabb', 'Morningsun', 'Morley', 'Mississinewa', 'Milton', 'Miamian', 'Miami', ' Lybrand', 'Loudonville', 'Losantville', 'Lewisburg', 'Hickory', 'Hennepin', 'Glynwood', 'Cliftycreek', 'Cincinnati', 'Celina', 'Cardington', 'Cana', 'Bonnell', 'Blocher', 'Birkbeck', 'Ava', 'Amanda', 'Alexandria', 'Kendallville', 'Cana Variant', 'Celina Variant', 'Miamian Variant', 'Milton variant'),# add 'Kendallville', 'Cana Variant', 'Celina Variant', 'Miamian Variant', 'Milton variant'
  R111AY010IN = c('Raub', 'Odell', 'Dana', 'Corwin'),
  F111AY011IN = c('Zipp', 'Peoga', 'Paulding', 'Patton', 'Montgomery', 'Minster', 'Milford', 'Latty'),
  F111AY012IN = c('Wyatt', 'Tyler', 'Shinrock', 'Omulga', 'Mentor', 'McGary', 'Markland', 'Lauer', 'Haubstadt', 'Glenford', 'Fulton', 'Dubois', 'Del Rey', 'Del Rey Variant'), # add 'Del Rey', 'Del Rey Variant'
  F111AY013IN = c('Sable', 'Ragsdale', 'Muren', 'Iva', 'Alford'),
  F111AY014IN = c('Whitaker', 'Waynetown', 'Thackery', 'Taggart', 'Starks', 'Sleeth', 'Shadeland', 'Savona', 'Roby', 'Rainsboro', 'Libre', 'Ionia', 'Homer', 'Haney', 'Fitchville', 'Digby', 'Thackery Variant'), # add 'Thackery Variant'
  F111AY015IN = c('Williamsburg', 'Wawaka', 'Spargus', 'Sisson', 'Shelocta', 'Rush', 'Pike', 'Parke', 'Ockley', 'Negley', 'Muncie', 'Mountpleasant', 'Martinsville', 'Gallman', 'Fox', 'Eldean', 'Chili', 'Chetwynd', 'Casco', 'Camden', 'Belmore'),
  R111AY016IN = c('Westland', 'Sebewa', 'Rensselaer', 'Pella', 'Millgrove', 'Mahalasville', 'Lippincott', 'Kane', 'Dunham', 'Drummer', 'Crane'),
  R111AY017IN = c('Wea', 'Waupecan', 'Warsaw', 'Tippecanoe', 'Rodman', 'Plattville', 'Nineveh', 'Lorenzo', 'Lickcreek', 'Eldean', 'Donnelsville', 'Warsaw Variant', 'Wea Variant'), # add 'Warsaw Variant', 'Wea Variant'
  F111AY018IN = c('Weikert', 'Rohan', 'Opequon', 'Gasconade', 'Corydon','Fairmount', 'Rock outcrop'),#add 'Fairmount', 'Rock outcrop'
  F111AY019IN = c('Muskingum', 'Lily', 'Latham', 'Jessietown', 'Gilwood', 'Gilpin', 'Edenton', 'Eden', 'Dekalb', 'Caneyville', 'Brownstown', 'Bratton', 'Berks'),
  F111AY020IN = c('Zenas', 'Zanesville', 'Wrays', 'Westmoreland', 'Wellston', 'Wellrock', 'Tarhollow', 'Stonehead', 'Muscatatuck', 'Grayford', 'Cruze', 'Coolville', 'Carmel', 'Brownsville', 'Boston'),
  F111AY021IN = c('Lyles', 'Bobtown', 'Ayrshire', 'Aquents'),
  F111AY022IN = c('Princeton', 'Boyer', 'Bloomfield', 'Alvin')
)



for(i in 1:length(ES111A.list)){
  ES111A.df0 <- as.data.frame(merge(names(ES111A.list[i]),ES111A.list[i])); colnames(ES111A.df0) <- c('ES','Series')
  if(i == 1){ES111A <- ES111A.df0}else{ES111A <- rbind(ES111A.df0, ES111A)}};rm(ES111A.df0) 

s.111A <- merge(ES111A, s,  by.x = 'Series', by.y = 'compname', all.y = T)
s.111A <- subset(s.111A, MLRA %in% '111A')

write.csv(s.111A, 'fy2021-refresh/mlra111A_sorts.csv', row.names = F)

#111B ----
#

ES111B.list <- list(
  F111BY101IN = c('Paulding','Bono', 'Bono variant', 'Coesse', 'Colwood', 'Hoytville', 'Latty', 'Lenawee', 'Luray',  'Maumee', 'McGuffey', 'Mermill', 'Milford', 'Minster', 'Montgomery', 'Olentangy', 'Patton', 'Patton Variant', 'Paulding', 'Pella',  'Rensselaer', 'Roundhead', 'Sebring', 'Southwest', 'Toledo'),
  F111BY102IN = c('Tiro','Steamburg','Seward','Selfridge','Lykens', 'Houcktown','Berrien', 'Cygnet', 'Darroch', 'Del Rey', 'Fitchville', 'Fulton', 'Glenford', 'Houcktown', 'Jenera', 'Kibbie', 'McGary', 'Nappanee', 'Ottokee', 'Rimer', 'Seward', 'Shawtown', 'Shinrock', 'Shinrock', 'till substratum', 'Tuscola', 'Vaughnsville', 'Ypsi'),
  F111BY201IN = c('Wabasha',  'Sloan', 'Saranac', 'Glendora', 'Cohoctah', 'Ceresco', 'Bellcreek', 'Wallkill'),
  F111BY202IN = c('Tice','Henshaw', 'Harrod', 'Spencerville', 'Rossburg', 'Ross', 'Medway', 'Landes', 'Armiesburg', 'Lash', 'Lickcreek'),
  F111BY203IN = c('Washtenaw', 'Shoals', 'Orrville', 'Newark', 'Defiance',  'Algiers', 'Griffin'),
  F111BY204IN = c('Tioga', 'Stonelick', 'Nolin', 'Lobdell', 'Lindside', 'Knoxdale', 'Gessie', 'Genesee', 'Flatrock', 'Eel', 'Chagrin'),
  F111BY302IN = c('Randolph', 'Prout'),
  F111BY303IN = c('Weikert', 'Ritchey', 'Newglarus', 'Milton', 'Latham', 'Heverlo', 'Eleva', 'Brecksville', 'Biglick', 'Berks'),
  F111BY403IN = c('Morocco', 'Lamberjack', 'Whitaker', 'Sleeth', 'Riverdale',  'Jimtown', 'Homer', 'Digby', 'Bogart', 'Aptakisic'),
  F111BY404IN = c('Pacer','Matherton','Thackery', 'Sisson', 'Scioto', 'Rush', 'Riverdale', 'Oshtemo', 'Ormas', 'Ockley', 'Muncie', 'Martinsville', 'Leoni', 'Kosciusko', 'Kalamazoo', 'Ionia', 'Haney', 'Gallman', 'Fox', 'Eldean', 'Chili', 'Casco', 'Bronson', 'Boyer', 'Belmore', 'Brems', 'Chelsea','Otisville', 'Plainfield', 'Spinks'),
  F111BY501IN = c('Treaty','Stone','Wolcott', 'Wetzel', 'Pewamo', 'Pandora', 'Kokomo', 'Brookston', 'Berville', 'Barry', 'Alvada'),
  F111BY502IN = c('Metamora', 'Markton', 'Macomb', 'Locke', 'Haskins', 'Elliott', 'Crosier', 'Conover', 'Capac', 'Blount', 'Aubbeenaubbee', 'Bennington', 'Crosby'),
  F111BY503IN = c('Wawaka', 'Wapahani', 'St. Clair', 'Muncie', 'Williamstown', 'Wawasee', 'Strawn',  'Russell', 'Riddles', 'Rawson', 'Parr', 'Owosso', 'Mortimer', 'Morley', 'Mississinewa', 'Miami', 'Metea', 'Marlette', 'Lybrand', 'Kidder', 'Kendallville',  'Hillsdale', 'Hennepin', 'Glynwood', 'Cadmus', 'Amanda', 'Miamian'),
  R111BY001IN = c('Willette','McGuffey','Roundhead', 'Palms', 'Olentangy', 'Ogden',  'Linwood', 'Kerston', 'Adrian', 'Tawas'),
  R111BY002IN = c('Martisco Variant','Warners', 'Rollin', 'Muskego', 'Martisco', 'Edwards'),
  R111BY003IN = c('Napoleon', 'Houghton', 'Carlisle', 'Boots'),
  R111BY301IN = c('Castalia', 'Millsdale', 'Joliet', 'Dunbridge', 'Channahon', 'Fries'),
  R111BY401IN = c('Westland', 'Sebewa', 'Rensselaer', 'Millgrove', 'Lippincott', 'Gilford', 'Granby'),
  R111BY402IN = c('Wasepi', 'Warsaw', 'Rodman', 'Nineveh', 'Kane', 'Crane')
)
for(i in 1:length(ES111B.list)){
  ES111B.df0 <- as.data.frame(merge(names(ES111B.list[i]),ES111B.list[i])); colnames(ES111B.df0) <- c('ES','Series')
  if(i == 1){ES111B <- ES111B.df0}else{ES111B <- rbind(ES111B.df0, ES111B)}};rm(ES111B.df0) 

s.111B <- merge(ES111B, s,  by.x = 'Series', by.y = 'compname', all.y = T)
s.111B <- subset(s.111B, MLRA %in% '111B')

write.csv(s.111B, 'fy2021-refresh/mlra111B_sorts.csv', row.names = F)

#111C ----

#0 clear sites --------------------------------------------------------- 
s$Site111 <-"Not"



#1 Lake and beaches --------------------------------------------------------- 
s$Site111<-ifelse(s$Site111 %in% "Not",
                  ifelse(s$muname %in% c("Water","water") & s$compname %in% c("Water","water"),"Water",
                         ifelse(s$compname %in% c("Beaches","Lake beach","Lake bluffs", "Dune land")| s$muname %in% c("Beach and Dune sand","Stony lake beaches","Dune land","Sand dunes"),"Shoreline Complex","Not"))
                  ,s$Site111)
#2  Mucks ---------------------------------------------------------               
s$Site111 <- ifelse(s$Site111 %in% "Not", ifelse(grepl("Histosols",s$compname)|grepl("ists",s$taxsubgrp)|grepl("ists",s$taxclname)|grepl("histic",s$taxsubgrp)|grepl("organic",s$pmkind),"Muck",s$Site111),s$Site111)

s$Site111 <- ifelse(s$Site111 %in% "Muck", ifelse(grepl("terric",s$taxsubgrp)|grepl("thapto",s$taxsubgrp)|grepl("terric",s$taxclname)|grepl("epts",s$taxsubgrp) ,"R111CY012IN",s$Site111),s$Site111)
s$Site111 <- ifelse(s$Site111 %in% "Muck", ifelse(grepl("limnic",s$taxsubgrp)|grepl("limnic",s$taxclname)|grepl("epts",s$taxsubgrp) ,"R111CY011IN", "R111CY013IN"),s$Site111)


#3 lacustrine Parent Material  --------------------------------------------------------- 
# s$Site111<-ifelse(s$Site111 %in% "Not",
#                   ifelse(grepl("lacustrine",s$pmkind)|grepl("lake plain",s$pmkind), 'lacustrine', s$Site111), s$Site111)
# 
# s$Site111<-ifelse(s$Site111 %in% "lacustrine",
#                   ifelse(s$drainagecl %in% c("very poorly", "poorly"), 'Lacustrine Flatwood', 'Lacustrine Forest'), s$Site111)

#3 alluvium Parent Material  --------------------------------------------------------- 
s$Site111<-ifelse(s$Site111 %in% "Not",
                  ifelse(grepl("alluvium",s$pmkind), 'alluvium', s$Site111), s$Site111)

s$Site111<-ifelse(s$Site111 %in% "alluvium",
                  ifelse(s$drainagecl %in% c("very poorly", "poorly", "somewhat poorly"), 'F111CY014IN', 'F111CY015IN')
                  , s$Site111)

#4 bedrock Parent Material  --------------------------------------------------------- 
# s$Site111<-ifelse(s$Site111 %in% "Not",
#                   ifelse(s$rockdepth < 200, 'bedrock', s$Site111), s$Site111)
# 
# s$Site111<-ifelse(s$Site111 %in% "bedrock",
#                   ifelse(grepl("Mollisols",s$compname)|grepl("olls",s$taxsubgrp)|grepl("olls",s$taxclname), 'Dark Bedrock Prairie',
#                          ifelse(s$drainagecl %in% c("very poorly", "poorly", "somewhat poorly"), 'Mesic Bedrock Forest', 'Dry Bedrock Forest'))
#                   , s$Site111)


#6 till Parent Material  --------------------------------------------------------- 
s$Site111<-ifelse(s$Site111 %in% "Not",
                  ifelse((grepl("till",s$pmkind) & !grepl("fluv",s$pmkind) & !grepl("outwash",s$pmkind))|
                           (grepl("drift",s$pmkind) & !grepl("fluv",s$pmkind) & !grepl("outwash",s$pmkind))|
                           s$landform_string %in% 'pothole & till plain',
                         'till', s$Site111), s$Site111)

s$Site111<-ifelse(s$Site111 %in% "till",
                  ifelse(s$drainagecl %in% c("very poorly", "poorly"),'R111CY005IN',
                         ifelse(!(grepl('humic',s$taxsubgrp) | grepl('ollic',s$taxsubgrp) | grepl('olls',s$taxsubgrp)),'F111CY007IN','R111CY006IN'))
                  , s$Site111)
#7 Sand Parent Material  --------------------------------------------------------- 
s$Site111<-ifelse(s$Site111 %in% "Not",
                  ifelse(s$T50_sand >= 70 & (s$T150_sand >= 80 | grepl('sand',s$pmkind))| s$T50_sand >= 80, 'sand', s$Site111), s$Site111)

s$Site111<-ifelse(s$Site111 %in% "sand",
                  ifelse(s$drainagecl %in% c("very poorly", "poorly"),'R111CY002IN',
                         ifelse(s$drainagecl %in% c("moderately well ", "somewhat poorly"),'F111CY003IN','R111CY001IN'))
                  , s$Site111)
#5 outwash Parent Material  --------------------------------------------------------- 
s$Site111<-ifelse(s$Site111 %in% "Not",
                  ifelse(grepl("outwash",s$pmkind)|grepl("fluv",s$pmkind)|grepl("lacustrine",s$pmkind)|grepl("lake plain",s$pmkind)|grepl("outwash",s$landform_string), 'outwash', s$Site111), s$Site111)

s$Site111<-ifelse(s$Site111 %in% "outwash",
                  ifelse(s$drainagecl %in% c("very poorly", "poorly"),'R111CY008IN',
                         ifelse(s$drainagecl %in% c("moderately well ", "somewhat poorly"),'F111CY009IN','R111CY010IN'))
                  , s$Site111)


#8 missing ----
seriessites <- c('R111CY010IN:Alvin', 'not:Aquolls', 'not:Beaches', 'F111CY009IN: Darroch', 'R111CY008IN: Faxon', 'F111CY009IN: Foresman', 'R111CY006IN: Francesville', 'F111CY007IN: Haskins', 'F111CY009IN: Headlee', 'Iroquois', 'Lucas', 'Marsh', 'Maumee', 'Medaryville', 'Milford', 'Millsdale', 'Monon', 'Montgomery', 'Nappanee', 'Navunon', 'Newglarus', 'Papineau', 'Pella', 'Peotone', 'Pits', 'Radioville', 'Rockton', 'Selma', 'Simonin', 'Strole', 'Toledo', 'Udorthents', 'Wesley', 'Whiskerville', 'Whitepost', 'Wolcott')


# s$Site111<-ifelse(grepl("^F111C", Site111)|grepl("^R111C", Site111),
#                   ifelse(s$compname %in% c(),'R111CY001IN'
#                          , s$Site111), s$Site111)










write.csv(s.111D, 'fy2021-refresh/mlra111D_sorts.csv', row.names = F)





s.111C <- subset(s, MLRA %in% '111C')

write.csv(s.111C, 'fy2021-refresh/mlra111C_sorts.csv', row.names = F)

#111D ----

#0 clear sites --------------------------------------------------------- 
s$Site111 <-"Not"



#1 Lake and beaches --------------------------------------------------------- 
s$Site111<-ifelse(s$Site111 %in% "Not",
                  ifelse(s$muname %in% c("Water","water") & s$compname %in% c("Water","water"),"Water",
                         ifelse(s$compname %in% c("Beaches","Lake beach","Lake bluffs", "Dune land")| s$muname %in% c("Beach and Dune sand","Stony lake beaches","Dune land","Sand dunes"),"Shoreline Complex","Not"))
                  ,s$Site111)
#2  Mucks ---------------------------------------------------------               
s$Site111 <- ifelse(s$Site111 %in% "Not", ifelse(grepl("Histosols",s$compname)|grepl("ists",s$taxsubgrp)|grepl("ists",s$taxclname)|grepl("organic",s$pmkind),"Muck",s$Site111),s$Site111)

s$Site111 <- ifelse(s$Site111 %in% "Muck", ifelse(grepl("terric",s$taxsubgrp)|grepl("terric",s$taxclname),"R111DY001IN",s$Site111),s$Site111)
s$Site111 <- ifelse(s$Site111 %in% "Muck", ifelse(grepl("limnic",s$taxsubgrp)|grepl("timnic",s$taxclname),"R111DY001IN", "R111DY002IN"),s$Site111)

s$Site111<-ifelse(s$Site111 %in% "Not", ifelse(s$rockdepth < 150, 'bedrock', s$Site111), s$Site111)

#3 lacustrine Parent Material  --------------------------------------------------------- 
s$Site111<-ifelse(s$Site111 %in% "Not",
                  ifelse(grepl("lacustrine",s$pmkind)|grepl("lake plain",s$pmkind), 'lacustrine', s$Site111), s$Site111)

s$Site111<-ifelse(s$Site111 %in% "lacustrine",
                  ifelse(s$drainagecl %in% c("very poorly", "poorly"), 'Lacustrine Flatwood', 'Lacustrine Forest'), s$Site111)

#3 alluvium Parent Material  --------------------------------------------------------- 
s$Site111<-ifelse(s$Site111 %in% "Not",
                  ifelse(grepl("alluvium",s$pmkind) | s$flood %in% 'flood', 'alluvium', s$Site111), s$Site111)

s$Site111<-ifelse(s$Site111 %in% "alluvium",
                  ifelse(s$drainagecl %in% c("very poorly", "poorly","somewhat poorly"), 'F111DY003IN', 'F111DY004IN')
                  , s$Site111)


#4 outwash Parent Material  --------------------------------------------------------- 
s$Site111<-ifelse(s$Site111 %in% "Not",
                  ifelse(grepl("outwash",s$pmkind)|grepl("fluv",s$pmkind), 'outwash', s$Site111), s$Site111)

s$Site111<-ifelse(s$Site111 %in% "outwash",
                  ifelse(s$drainagecl %in% c("very poorly", "poorly"),'R111CY008IN',
                         ifelse(s$drainagecl %in% c("moderately well ", "somewhat poorly"),'F111CY009IN','R111CY010IN'))
                  , s$Site111)

#5 till Parent Material  --------------------------------------------------------- 
s$Site111<-ifelse(s$Site111 %in% "Not",
                  ifelse(grepl("till",s$pmkind)|grepl("drift",s$pmkind), 'till', s$Site111), s$Site111)

s$Site111<-ifelse(s$Site111 %in% "till",
                  ifelse(s$drainagecl %in% c("very poorly", "poorly"),'R111CY005IN',
                         ifelse(s$slope_r > 4 & !(grepl('humic',s$taxsubgrp) | grepl('mollic',s$taxsubgrp) | grepl('olls',s$taxsubgrp)),'F111CY007IN','R111CY006IN'))
                  , s$Site111)
#6 Sand Parent Material  --------------------------------------------------------- 
s$Site111<-ifelse(s$Site111 %in% "Not",
                  ifelse(s$T50_sand > 70, 'sand', s$Site111), s$Site111)

s$Site111<-ifelse(s$Site111 %in% "sand",
                  ifelse(s$drainagecl %in% c("very poorly", "poorly"),'R111CY002IN',
                         ifelse(s$drainagecl %in% c("moderately well ", "somewhat poorly"),'F111CY003IN','R111CY001IN'))
                  , s$Site111)

s.111D <- subset(s, MLRA %in% '111D')

#7 bedrock Parent Material  --------------------------------------------------------- 
s$Site111<-ifelse(s$Site111 %in% "Not",
                  ifelse(s$rockdepth < 200 | grepl("residuum",s$pmkind), 'bedrock', s$Site111), s$Site111)

s$Site111<-ifelse(s$Site111 %in% "bedrock",
                  ifelse(s$rockdepth < 100, ifelse(s$rockdepth < 50, 'F111DY022IN', 'F111DY023IN'),'F111DY024IN')
                  , s$Site111)




ES111D.list <- list(F111DY003IN = c('Washtenaw', 'Wakeland', 'Vincennes', 'Tice', 'Southwest', 'Sloan', 'Shoals', 'Sawabash', 'Saranac', 'Rockmill', 'Petrolia', 'Henshaw', 'Comfrey', 'Cohoctah', 'Ceresco', 'Brouillett', 'Beaucoup', 'Algiers'),
                    F111DY004IN = c('Uniontown', 'Stringley', 'Stonelick', 'Sligo', 'Rossburg', 'Ross', 'Pinevillage', 'Ouiatenon', 'Moundhaven', 'Medway', 'Lobdell', 'Lash', 'Lanier', 'Landes', 'Jules', 'Huntsville', 'Hononegah', 'Gessie', 'Genesee', 'Elkinsville', 'Eel', 'Coblen', 'Chatterton', 'Chagrin', 'Beckville', 'Battleground', 'Armiesburg', 'Allison', 'Allison', 'Armiesburg', 'Battleground', 'Beckville'),
                    F111DY005IN = c('Yeddo', 'Reesville', 'Haskins', 'Evansville', 'Clermont', 'Blount', 'Blanchester'),
                    F111DY008IN = c('Treaty', 'Selma', 'Secondcreek', 'Pewamo', 'Pella', 'Millsdale', 'Lisbon', 'Kokomo', 'Cyclone', 'Cope', 'Chalmers', 'Brookston', 'Ashkum'),
                    F111DY009IN = c('Westboro', 'Sugarvalley', 'Schaffer', 'Randolph', 'Mitiwanga', 'Hoosierville', 'Fincastle', 'Crosier', 'Crosby'),
                    F111DY010IN = c('Xenia', 'Wynn', 'Williamstown', 'Wapahani', 'Tuscola', 'Thrifton', 'Strawn', 'Senachwine', 'Russell', 'Rossmoyne', 'Rockfield', 'Ritchey', 'Riddles', 'Richardville', 'Rawson', 'Rainsville', 'Ozaukee', 'Morningsun', 'Morley', 'Milton', 'Miamian', 'Miami', 'Manlove', 'Loudonville', 'Losantville', 'Jonesboro', 'Hickory', 'Hennepin', 'Glynwood', 'Edenton', 'Crouse', 'Cincinnati', 'Celina', 'Cadiz', 'Bonnell', 'Birkbeck', 'Ava'),
                    F111DY013IN = c('Zipp', 'Patton', 'Montgomery', 'Milford', 'Belleville'),
                    F111DY014IN = c('Del Rey', 'Haubstadt'),
                    F111DY015IN = c('Whitson', 'Sable', 'Ragsdale', 'Edwardsville'),
                    F111DY016IN = c('Muren', 'Iva', 'Alford'),
                    F111DY017IN = c('Whitaker', 'Waynetown', 'Thackery', 'Taggart', 'Starks', 'Sleeth', 'Shadeland', 'Markland', 'Libre', 'Kendall', 'Homer', 'Fitchville'),
                    F111DY018IN = c('Williamsburg', 'Wawaka', 'Spinks', 'Sparta', 'Sisson', 'Silverwood', 'Rush', 'Parke', 'Oshtemo', 'Ormas', 'Ockley', 'Negley', 'Mudlavia', 'Martinsville', 'Kosciusko', 'Kendallville', 'Kalamazoo', 'Fox', 'Eldean', 'Coloma', 'Casco', 'Camden', 'Boyer', 'Angatoka'),
                    F111DY022IN = c('Weikert', 'Opequon', 'Gasconade', 'Fairmount', 'Corydon'),
                    F111DY023IN = c('Newglarus', 'Lumberton', 'Judyville', 'Gosport', 'Fairmount', 'Eden', 'Cates', 'Bratton', 'Berks', 'Adeland'),
                    F111DY024IN = c('Woolper', 'Switzerland', 'Pate', 'Morrisville', 'Loudon', 'Lawshe', 'Grayford', 'Carmel', 'Boston'),
                    F111DY025IN = c('Ayrshire'),
                    R111DY001IN = c('Wallkill', 'Palms', 'Muskego', 'Linwood', 'Adrian', 'Ackerman'),
                    R111DY002IN = c('Houghton','Carlisle'),
                    R111DY006IN = c('Lauramie', 'Conover'),
                    R111DY007IN = c('Tecumseh', 'Sidell', 'Linkville', 'Jasper'),
                    R111DY011IN = c('Wingate', 'Toronto', 'Throckmorton', 'Octagon', 'Montmorenci', 'Mellott', 'Markham'),
                    R111DY012IN = c('Williamsport', 'Varna', 'Symerton', 'Raub', 'Parr', 'Odell', 'Elliott', 'Dana', 'Corwin', 'Barce'),
                    R111DY019IN = c('Seafield', 'Mulvey', 'Millbrook', 'Longlois', 'Glenhall', 'Desker', 'Bowes', 'Billett'),
                    R111DY020IN = c('Westland', 'Sebewa', 'Rensselaer', 'Peotone', 'Millgrove', 'Mahalasville', 'Mahalaland', 'Lafayette', 'Harpster', 'Gilford', 'Gilboa', 'Free', 'Dunham', 'Drummer', 'Darroch', 'Crane', 'Brenton', 'Andres'),
                    R111DY021IN = c('Wea', 'Waupecan', 'Warsaw', 'Troxel', 'Totanang', 'Tippecanoe', 'Shipshe', 'Rodman', 'Proctor', 'Plattville', 'Lorenzo', 'Foresman', 'Elston', 'Chetwynd', 'Carmi', 'Barce'),
                    R111DY026IN = c('Ade'),
                    R111DY027IN = c('Princeton', 'Chelsea', 'Alvin', 'Bloomfield'))

for(i in 1:length(ES111D.list)){
  ES111D.df0 <- as.data.frame(merge(names(ES111D.list[i]),ES111D.list[i])); colnames(ES111D.df0) <- c('ES','Series')
  if(i == 1){ES111D <- ES111D.df0}else{ES111D <- rbind(ES111D.df0, ES111D)}};rm(ES111D.df0) 

s.111D <- merge(ES111D, s,  by.x = 'Series', by.y = 'compname', all.y = T)
s.111D <- subset(s.111D, MLRA %in% '111D')

write.csv(s.111D, 'fy2021-refresh/mlra111D_sorts.csv', row.names = F)


#111E ----

ES111E.list <- list(
  
  F111EY101OH = c('Sebring', 'Minster', 'Milford', 'Luray', 'Lenawee', 'Colwood', 'Bono', 'Canadice', 'Lorain', 'Patton'),
  F111EY102OH = c('Del Rey', 'Bixler','Tuscola', 'Shinrock', 'Mentor', 'Kibbie', 'Glenford', 'Fitchville', 'Canal', 'Caneadea', 'Elnora', 'Henshaw', 'McGary'),
  F111EY201OH = c('Sloan', 'Saranac', 'Rockmill'),
  F111EY202OH = c('Rossburg', 'Medway'),
  F111EY203OH = c('Orrville Variant', 'Shoals', 'Orrville', 'Newark', 'Killbuck', 'Holly', 'Aetna', 'Algiers', 'Beaucoup', 'Euclid'),
  F111EY204OH = c('Tioga', 'Lobdell', 'Lindside', 'Genesee', 'Chagrin', 'Eel', 'Gessie', 'Stonelick'),
  F111EY301OH = c('Smothers', 'Mitiwanga', 'Prout'),
  F111EY302OH = c('Colyer Variant', 'Rarden', 'Milton', 'Loudonville', 'Latham', 'Brecksville', 'Dekalb', 'Germano', 'Gilpin', 'Lordstown', 'Mechanicsburg', 'Steinsburg', 'Wakeman', 'Wellston', 'Tarlton'),
  F111EY403OH = c('Sleeth', 'Jimtown', 'Digby'),
  F111EY404OH = c('Wheeling', 'Spinks', 'Oshtemo', 'Ockley', 'Martinsville', 'Haney', 'Gallman', 'Fox', 'Chili',  'Belmore', 'Eldean', 'Pike', 'Thackery', 'Bogart'),
  F111EY501OH = c('Pewamo', 'Mermill', 'Marengo', 'Condit'),
  F111EY502OH = c('Tiro', 'Hyatts', 'Haskins', 'Bennington', 'Blount', 'Crosby', 'Elliott', 'Gresham', 'Mahoning'),
  F111EY503OH = c('Heverlo','Cana Variant', 'Lykens', 'Lybrand', 'Kendallville', 'Hennepin', 'Centerburg', 'Cardington', 'Amanda', 'Alexandria', 'Cincinnati', 'Corwin',  'Geeburg', 'Glynwood', 'Hickory', 'Jeneva', 'Miamian', 'Rawson', 'Thrifton'),
  R111EY001OH = c('Wallkill', 'Linwood', 'Willette'),
  R111EY002OH = c('Olentangy', 'Muskego'),
  R111EY003OH = c('Pinnebog', 'Carlisle'),
  R111EY401OH = c('Olmsted', 'Westland', 'Millgrove'),
  R111EY402OH = c('Crane','Wilmer Variant', 'Warsaw', 'Wea'),
  F111AY013IN = c('Alford')
)

for(i in 1:length(ES111D.list)){
  ES111E.df0 <- as.data.frame(merge(names(ES111E.list[i]),ES111E.list[i])); colnames(ES111E.df0) <- c('ES','Series')
  if(i == 1){ES111E <- ES111E.df0}else{ES111E <- rbind(ES111E.df0, ES111E)}};rm(ES111E.df0) 

s.111E <- merge(ES111E, s,  by.x = 'Series', by.y = 'compname', all.y = T)
s.111E <- subset(s.111E, MLRA %in% '111E')

write.csv(s.111E, 'fy2021-refresh/mlra111E_sorts.csv', row.names = F)

#114A ----

ES114A.list <- list(
  R111AY001IN = c('Adrian'),
  F114AY101IN = c('Luray', 'McGary', 'Montgomery', 'Peoga', 'Sebring', 'Zipp', 'Booker', 'Patton', 'Secondcreek'),
  F114AY102IN = c('Bartle', 'Dubois', 'Haubstadt', 'Otwell', 'Pekin', 'Fitchville', 'Henshaw'),
  F114AY103IN = c('Glenford', 'Elkinsville', 'Markland', 'Mentor', 'Millstone', 'Shircliff', 'Greybrook', 'Olephant', 'Pottersville', 'Stubenville'),
  F114AY203IN = c("Atkins","Belknap", "Birds","Bonnie", "Evansville", "Holton", "Petrolia", "Shoals", "Sloan", "Stendal", "Taggart", "Wakeland", "Whitaker", "Wilhite", 'Riverwash'),
  F114AY204IN = c('Algiers', 'Chagrin', 'Cuba', 'Eel', 'Genesee', 'Gessie', 'Hatfiled', 'Haymond', 'Jules', 'Lanier', 'Lobdell', 'Medway', 'Mondhaven', 'Newark', 'Oldenburg', 'Ross', 'Sciotoville', 'Steff', 'Stonelick', 'Wilbur', 'Wirt', 'Armiesburg', 'Huntsville', 'Lindside', 'McAdoo', 'Moundhaven', 'Nolin', 'Piankeshaw', 'Pope', 'Sligo', 'Waupecan'),
  F114AY302IN = c('Deputy', 'Gilpin', 'Jennings', 'Jessietown', 'Loudon', 'Loudonville', 'Mechanicburg', 'Rarden Scottsburg Trappist', 'Weddel', 'Whitcomb', 'Adyeville', 'Berks', 'Cana', 'Colyer', 'Ebal', 'Fairpoint', 'Johnsburg', 'Muse', 'Muskingum', 'Neotoma', 'Patricksburg', 'Shadeland', 'Tipsaw', 'Trappist', 'Tulip', 'Tuscarawas', 'Weikert', 'Wellston', 'Zanesville'),
  F114AY305IN = c('Caneyville', 'Carmel', 'Crider', 'Eden', 'Edenton', 'Juessup', 'Morrisville', 'Switzerland', 'Zenas', 'Boston', 'Bratton', 'Corydon', 'Fairmount', 'Gasconade', 'Hagerstown', 'Lawshe', 'Lowell', 'Opequon', 'Romona', 'Sees'),
  F114AY404IN = c('Chetwynd', 'Chili', 'Libre', 'Medona', 'Negley', 'Ninevah', 'Parke', 'Pike', 'Rainsboro', 'Sardinia', 'Vallonia', 'Williamsburg', 'Alvin', 'Casco', 'Dunham', 'Elston', 'Fox', 'Gallimore', 'Lyles', 'Mahalasville', 'Martinsville', 'Ockley', 'Rensselaer', 'Roby', 'Rodman', 'Sleeth','Princeton', 'Thackery', 'Wea', 'Westland'),
  F114AY501IN = c('Atlas', 'Avonburg', 'Blanchester', 'Clermont', 'Cobbsfork', 'Brookston', 'Hoosierville', 'Schaffer'),
  F114AY502IN = c('Blocher', 'Cincinnati', 'Homewood', 'Jonesboro', 'Nabb', 'Nicely', 'Rossmoyne', 'Ryker', 'Titusville'),
  F114AY504IN = c('Bonnell', 'Grayford', 'Hickory', 'Ryker', 'Weisburg', 'Arney', 'Ava', 'Camden', 'Lieber', 'Miamian', 'Milton', 'Russell', 'Shakamak'),
  F114AY802IN = c('Ayrshire', 'Bloomfield', 'Bobtown', 'Priceton', 'Cory', 'Hosmer', 'Iona', 'Iva','Jessup', 'Muren', 'Potawatomi', 'Solsberry','Stinesville', 'Alford', 'Vigo', 'Westboro', 'Whitson')
)


for(i in 1:length(ES114A.list)){
  ES114A.df0 <- as.data.frame(merge(names(ES114A.list[i]),ES114A.list[i])); colnames(ES114A.df0) <- c('ES','Series')
  if(i == 1){ES114A <- ES114A.df0}else{ES114A <- rbind(ES114A.df0, ES114A)}};rm(ES114A.df0) 

s.114A <- merge(ES114A, s,  by.x = 'Series', by.y = 'compname', all.y = T)
s.114A <- subset(s.114A, MLRA %in% '114A')

write.csv(s.114A, 'fy2021-refresh/mlra114A_sorts.csv', row.names = F)


#114B ----

ES114B.list <- list(
  
  F114BY103IN = c('Floraville', 'Lakaskia', 'McGary', 'Montgomery', 'Okaw', 'Patton', 'Peoga', 'Wagner', 'Zipp'),
  F114BY104IN = c('Bartle', 'Colp', 'Dubois', 'Elkinsville', 'Hurst', 'Greybrook', 'Haubstadt', 'Millstadt', 'Markland', 'Otwell', 'Olephant', 'Potterville', 'Redbud', 'Stubenville', 'Uniontown'),
  F114BY203IN = c('Atkins', 'Beaucoup', 'Birds', 'Bonnie', 'Coffeen', 'Colo', 'Holton', 'Orion', 'Otter', 'Petrolia', 'Shoals', 'Stendal', 'Tice', 'Titus', 'Wabash', 'Wakeland', 'Wilhite', "Algiers", "Banlic", "Bartelso", "Bellcreek", "Blackoar", "Booker", "Creal", "Driftwood", "Dupo", "Jacob", "Lawson", "Orrville", "Piopolis", "Racoon", "Radford", "Vincennes"),
  F114BY204IN = c('Armiesburg', 'Beanblossum', 'Chagrin', 'Cuba', 'Eel', 'Elkinsville', 'Genesee', 'Haymond', 'Huntsville', 'Landes', 'Lindside', 'Lobell', 'McAdoo', 'Moundhaven', 'Nolin', 'Oldenburg', 'Pekin', 'Ross', 'Steff', 'Stonelick', 'Terril', 'Wilbur', 'Wirt', "Blyton", "Dearborn", "Gessie", "Huntington", "Lobdell", "Millstone", "Nineveh", "Ridgway", "Vallonia"),
  F114BY302IN = c('Grayford', 'Ryker', 'Stinesville', 'Wellston','Blocher', 'Caneyville', 'Carmel', 'Coolville', 'Corydon', 'Crider', 'Deputy', 'Eden', 'Edenton', 'Gosport', 'Haggatt', 'Jennings', 'Jessietown', 'Kell', 'Kurtz', 'Milton', 'Nicholson', 'Rarden', 'Rohan', 'Scottsburg', 'Stonehead', 'Switzerland', 'Trappist', 'Weddel', 'Whitcomb', 'Zenas', 'Schuline'),
  F114BY403IN = c('Whitaker', 'Taggart', 'Rensselaer', 'Roby', 'Lyles', 'Kendall', 'Geff'),
  F114BY404IN = c('Ridgeway', 'Rend', 'Pike', 'Parke', 'Ockley', 'Negley', 'Martinsville', 'Gallimore', 'Chetwynd', 'Campton', 'Camden', 'Fox'),
  F114BY502IN = c('Vigo', 'Hoosierville', 'Lieber', 'Blair','Atlas'),
  F114BY503IN = c('Arney', 'Ava', 'Cincinnati', 'Elco', 'Grantfork', 'Hickory', 'Solsberry', 'Shakamak', 'Bonnell', 'Ursa'),
  F114BY801IN = c('Princeton', 'Bloomfield', 'Ayrshire', 'Alvin', 'Bobtown'),
  F114BY803IN = c('Bunkum', 'Cory', 'Iva', 'Avonburg', 'Caseyville', 'Clarksdale', 'Cobbsfork', 'Emery', 'Fishhook', 'Keomah', 'Lauer', 'Millbrook', 'Sexton', 'Starks', 'Stoy', 'Stronghurst', 'Weir'),
  F114BY804IN = c('Alford', 'Homen', 'Hosmer', 'Menfro', 'Muren', 'Rumen', 'Winfield', 'Bedford', 'Medora', 'Muscatatuck', 'Nabb', 'Rossmoyne', 'Rozetta', 'Ruma', 'Shircliff', 'Weisburg'),
  F114BY805IL = c('Marine', 'Pierron', 'Burksville', 'Rushville'),
  R114BY901IN = c('Biddle', 'Coulterville', 'Darmstadt', 'Fosterburg', 'Piasa', 'Tamalco', 'Huey'),
  R114BY902IN = c('Bethalto', 'Cowden', 'Edwardsville', 'Herrick', 'Hoyleton', 'Keller', 'Mascoutah', 'Shiloh', 'Oconee', 'Virden', 'Chauncey'),
  R114BY903IN = c('Aviston', 'Douglas', 'Elston', 'Harrison', 'Meadowbank', 'Pana', 'Waupecan', 'Wakenda', 'Wea', 'Velma', 'Assumption', 'Corydon', 'Fairmount', 'Rodman', 'Woolper')
  
)

q = unique(
  subset(missing.114B, grepl('eolian', pmkind) & T50_sand >65)$Series
)

cat(q, sep = "\', \'")

q = sort(unique(
  missing.111C$compname
))

cat(q, sep = "\', \'")





for(i in 1:length(ES114B.list)){
  ES114B.df0 <- as.data.frame(merge(names(ES114B.list[i]),ES114B.list[i])); colnames(ES114B.df0) <- c('ES','Series')
  if(i == 1){ES114B <- ES114B.df0}else{ES114B <- rbind(ES114B.df0, ES114B)}};rm(ES114B.df0) 

s.114B <- merge(ES114B, s,  by.x = 'Series', by.y = 'compname', all.y = T)
s.114B <- subset(s.114B, MLRA %in% '114B')

write.csv(s.114B, 'fy2021-refresh/mlra114B_sorts.csv', row.names = F)


stringr::str_count('cou nt ', "\\s")

missing.111A <- subset(s.111A, is.na(ES) & stringr::str_count(Series, "\\s")<1 & majcompflag %in% 'TRUE')
missing.111B <- subset(s.111B, is.na(ES) & stringr::str_count(Series, "\\s")<1 & majcompflag %in% 'TRUE')
missing.111C <- subset(s.111C, !(grepl("^F111C", Site111)|grepl("^R111C", Site111)) & stringr::str_count(compname, "\\s")<1 & majcompflag %in% 'TRUE')
missing.111D <- subset(s.111D, is.na(ES) & stringr::str_count(Series, "\\s")<1 & majcompflag %in% 'TRUE')
missing.111E <- subset(s.111E, is.na(ES) & stringr::str_count(Series, "\\s")<1 & majcompflag %in% 'TRUE')
missing.114A <- subset(s.114A, is.na(ES) & stringr::str_count(Series, "\\s")<1 & majcompflag %in% 'TRUE')
missing.114B <- subset(s.114B, is.na(ES) & stringr::str_count(Series, "\\s")<1 & majcompflag %in% 'TRUE')

