library(soilDB)
library(aqp)
library(plyr)
library(foreign)
###2021-11-12


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
mu_tgs <- merge(mu, tgs_zonal, by.x =  'lmapunitiid', by.y = 'VALUE')
mu_tgs$tgswt <- mu_tgs$MEAN*mu_tgs$COUNT
mu_tgs <- aggregate(list(tgswt=mu_tgs$tgswt, wt=mu_tgs$COUNT), by=list(muiid =mu_tgs$muiid ), FUN='sum')
mu_tgs$tgs <- mu_tgs$tgswt/mu_tgs$wt
mu <-  merge(mu, mu_tgs[,c('muiid','tgs')], by='muiid', all.x = TRUE)
mlra$MLRA.total <- mlra$Count
mlra$MLRA.value <- mlra$Value

mlra.mu <- merge(mlra[,c('MLRA.value','MLRARSYM', 'MLRA.total')], mlra.mu[,c('Value', 'Count', 'Feature_CONU2', 'northeast_gssurgo30m')], by.x='MLRA.value', by.y='Feature_CONU2')
mlra.mu <- merge(mlra.mu, mu[,c('lmapunitiid', 'muiid')], by.x='northeast_gssurgo30m', by.y='lmapunitiid')
mlra.mu$percentmlra <- mlra.mu$Count/mlra.mu$MLRA.total*100
mlra.mu.agg <- aggregate(list(pmlra = mlra.mu$percentmlra), by=list(muiid = mlra.mu$muiid, MLRARSYM = mlra.mu$MLRARSYM), FUN='sum')
mlra.mu <- merge(mlra.mu, mlra.mu.agg, by=c('muiid', 'MLRARSYM'))
mlra.mu.totals <- aggregate(list(mu.total=mlra.mu.agg$pmlra), by= list(mu = mlra.mu.agg$muiid), FUN ='sum')
mlra.mu <- merge(mlra.mu, mlra.mu.totals, by.x='muiid', by.y='mu')
mlra.mu$affinity <- mlra.mu$pmlra/mlra.mu$mu.total * 100
mlra.mu$LRU <- mlra.mu$MLRARSYM
mlra.best <- aggregate(list(best=mlra.mu$affinity), by= list(mu = mlra.mu$muiid), FUN ='max')
mlra.mu <- merge(mlra.mu, mlra.best, by.x='muiid', by.y='mu')
mlra.mu$ofbest <- mlra.mu$affinity / mlra.mu$best *100
mlra.mu.mlra <- unique(mlra.mu[,c('muiid', 'LRU')])
mlra.mu.counts <-  aggregate(list(count=mlra.mu.mlra$LRU), by= list(mu = mlra.mu.mlra$muiid), FUN ='length')
mlrabestlink <-  subset(mlra.mu, ofbest > 99)

######################################----
#usergroup <- get_group_from_table("datamapunit")
#saveRDS(usergroup, 'fy2021-refresh/usergroup20211215.RDS')
usergroup <- readRDS("fy2021-refresh/usergroup20211215.RDS")
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
compsorts <- merge(mu[,c('areasymbol','areaname','lmapunitiid','dmuiid','nationalmusym','muname','tgs')], compsorts, by.x = 'dmuiid', by.y = 'dmuiid', all.y = F )
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
                               ifelse(s$drainagecl %in% c("somewhat poorly","poorly","very poorly","subaqueous"), 50, s$Water_Table)),
                        ifelse(s$hydricrating %in% "yes",
                               ifelse(s$drainagecl %in% c("very poorly","subaqueous"),0,25),
                               ifelse(s$drainagecl %in% c("very poorly","subaqueous"),0,
                                      ifelse(s$drainagecl %in% "poorly",25,
                                             ifelse(s$drainagecl %in% "somewhat poorly",50,
                                                    ifelse(s$drainagecl %in% "moderately well",100,
                                                           250))))))



rm(compsorts, fc.hz.bedrock,  fc.hz.bedrock.min, fc.hz.carb, fc.hz50.total, fc.hz.carb.min, fc.hz.notnull, fc.hz.sum, fc.hz150.total, cm.flood, fc.hz.humic.max, fc.hz.humic, fc.hz.mineral )


summarydrcl <- aggregate(s[,c('drainagecl')], by=list(s$drainagecl, s$hydricrating,s$Water_Table), FUN = 'length')

s <- merge(mlra.mu, s, by.x = 'northeast_gssurgo30m', by.y = 'lmapunitiid')
names(s)[names(s) == 'northeast_gssurgo30m'] <- 'lmapunitiid'
s$MLRA.total <- NULL;s$MLRA.value <- NULL;s$Value <- NULL;#s$Count <- NULL



s.reserve <- s
s.lmu <- unique(subset(s, select =c(muiid, lmapunitiid)))


s <- s.reserve
mlra139.missing <- read.csv('fy2021-refresh/mlra139missingesd.csv'); colnames(mlra139.missing) <- c('lmukey','compname')
#s <- merge(s,mlra139.missing, by.x = c('lmapunitiid','compname'), by.y = c('lmukey','compname'), all.y = TRUE)
#s <- subset(s,  majcompflag %in% 'TRUE' & ofbest > 99)
s <- subset(s, majcompflag %in% 'TRUE' & ofbest > 99)

s.mlrasummary <- aggregate(s$MLRARSYM, by=list(s$MLRARSYM), FUN='length')
#soil sort for MLRA 139 ----
#0 clear sites______________________________________________________________
s$Site<-"Not"



#1 Lake and beaches______________________________________________________________
s$Site<-ifelse(s$Site %in% "Not",
               ifelse(s$muname %in% c("Water","water") & s$compname %in% c("Water","water"),"Water",
                      ifelse(s$compname %in% c("Beaches","Lake beach","Lake bluffs", "Dune land")| s$muname %in% c("Beach and Dune sand","Stony lake beaches","Dune land","Sand dunes"),"F139XY001OH","Not"))
               ,s$Site)
s$Site<-ifelse(s$Site %in% "Not",
               ifelse(s$compname %in% c("Rock outcrop"),"F139XY007OH","Not"),s$Site)
s$Site<-ifelse(s$Site %in% "Not",
               ifelse(s$compname %in% c("Riverwash"),"F139XY009OH","Not"),s$Site)
#2 Bedrock Floodplains Mucks_____________________________________________________________                
s$Site<-ifelse(s$Site %in% "Not",
               ifelse(s$muname %in% c("flooded") | s$compname %in% c("Alluvial land")|grepl("flood",s$landform_string)|s$flood %in% 'flood' | grepl('flood', s$landform_string),
                      ifelse(s$Water_Table>=50,"F139XY008OH","F139XY009OH"),
                      ifelse((s$Water_Table<0 & s$T150_OM >=20)|((grepl("Histic",s$taxsubgrp)|grepl("Histosols",s$compname)|grepl("ists",s$taxsubgrp)|grepl("ists",s$taxclname))& s$Water_Table<=0),"Muck",
                             ifelse(s$rockdepth < 150 & s$Water_Table >100,"F139XY007OH","Not")))
               ,s$Site)
#3 Sand-Till---------------_____________________________________________________________  
s$Site<-ifelse(s$Site %in% "Not",
               ifelse(
                 (((s$T150_sand >= 80 & s$T50_sand >= 70)|(s$T50_sand >= 80)|(s$T150_clay < 20 & s$T50_pH < 6 & grepl('ult', s$taxsubgrp))) & !(TRUE)), "Sandy",
                 ifelse((TRUE) & 
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
               ifelse(s$Water_Table<=25,"F139XY011OH",
                      ifelse(s$Water_Table<=50,"F139XY002OH",
                             ifelse(s$Water_Table<=100,"F139XY002OH",
                                    ifelse(s$slope_r >=15,"F139XY003OH", "F139XY003OH"))))
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
               ifelse(s$Water_Table<=25,"F139XY012OH",
                      ifelse(s$tgs < 16 & !is.na(s$tgs), "F139XY006OH",
                             ifelse(s$Water_Table<=50,"F139XY004OH",
                                    ifelse(s$Water_Table<=100,"F139XY004OH",
                                           ifelse(s$slope_r >=15,"F139XY004OH", "F139XY005OH")))))
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
               ifelse(!is.na(s$T50_pH), ifelse(s$T50_pH < 5.0, "F139XY014OH", "F139XY013OH"), ifelse(grepl('dysic',s$taxclname), "F139XY014OH", "F139XY013OH")),s$Site)


unplaced <- subset(s, is.na(Site))
s.nolmapunit<- unique(subset(s, select=-c(lmapunitiid,percentmlra,pmlra,mu.total)))
s.139.owned <- merge(s, usergroup[,c('grpname','dmuiid')], by='dmuiid')
s.139.owned <- subset(s.139.owned, grpname %in% '12-BEL Belmont, New York')

write.csv(s.139.owned,'fy2021-refresh/s.139.csv', row.names = F)


not139 <- merge(s, usergroup[,c('grpname','dmuiid')], by='dmuiid')
not139 <- subset(not139,   !grpname %in% '12-BEL Belmont, New York')

write.csv(not139, 'fy2021-refresh/not139.csv', row.names = F)



## make maps ----
library(terra)
library(sf)
soilgrid <- rast('D:/GIS/SOIL/2021/northeast_gssurgo90m.tif')

mlramap <- sf::read_sf('C:/a/Ecological_Sites/GIS/Ecoregion/RegionalReview_296/CONUS.shp')
mlra139 <- subset(mlramap, MLRARSYM %in% '139')

extof139 <- ext(mlra139)
extof139a <- extof139
for(i in 1:4){
extof139a[i] <- extof139a[i] + c(-100000, 100000, -100000, 100000)[i]}

rast139 <- crop(soilgrid, extof139a)

rast139[rast139 == 172128] <- 0
rast139[rast139 %in% c(170831, 170342)] <- 1
rast139[!rast139 %in% c(0,1)] <- -1
plot(rast139)

listlmu <- s.reserve[s.reserve$MLRARSYM %in% '139' & s.reserve$ofbest <50 & !is.na(s.reserve$taxorder), 'lmapunitiid']
rast139 <- crop(soilgrid, extof139a)
rast139[rast139 %in% listlmu] <- 1
rast139[!rast139 %in% c(0,1)] <- -1

usergrouplmuid <- merge(s.reserve, usergroup[,c('grpname','dmuiid')], by='dmuiid')
groupsummary <- unique(usergrouplmuid[usergrouplmuid$MLRARSYM %in% '139' & usergrouplmuid$ofbest > 99 & usergrouplmuid$affinity > 50& !is.na(s.reserve$taxorder),]$grpname)


rast139 <- crop(soilgrid, extof139a)
for (i in 1:length(groupsummary)){#i=3
  listlmu <- usergrouplmuid[usergrouplmuid$MLRARSYM %in% '139' & usergrouplmuid$ofbest > 99 & usergrouplmuid$affinity > 00 & usergrouplmuid$grpname %in% groupsummary[i] & !is.na(s.reserve$taxorder),'lmapunitiid']
  rast139[rast139 %in% listlmu] <- i
};rast139[!rast139 %in% 1:length(groupsummary)] <- -1
plot(rast139)
writeRaster(rast139, 'fy2021-refresh/rast139a.tif', overwrite=T)

###missing 99 ----
library(terra)
library(sf)
missing99 <- c(889672, 889676, 889672, 889676, 251810, 272002, 272006, 272010, 245294, 251839, 251839, 245295, 271358, 271359, 271364, 245074, 1846465, 1846467, 1846507, 248050, 889689, 889689, 245247, 245248, 1846530, 1846533, 2723039, 1846544, 1846547, 251979, 2723040, 2723038, 2723040, 246258, 924405, 923986, 924405, 1846596, 271434, 271435, 271986, 271986, 271436, 271438, 271440, 271440, 271440, 245329, 245329, 245334, 245334, 923985, 924021, 923987, 1846615, 1846622, 1846626, 1846628, 1846629, 1846618, 1847686, 1847692, 1847679, 271474, 245349, 924160, 2106604, 245351, 245351, 252324, 271532, 271533, 271534, 271535, 924320, 245369, 245369, 1183800, 1183800, 252409, 924032, 924404, 924404, 924419, 889724, 889724, 252459, 1847000, 271540, 271542, 271546, 271548, 275725, 272064, 271926, 272064, 272066, 2088059, 2088063, 2088064, 275725, 275726, 275727, 275728, 272136, 245366, 245366, 245367, 1847015, 1847015, 271975, 272078, 2088124, 924652, 271466, 271563, 245382, 271566)

soilgrid <- rast('D:/GIS/SOIL/2021/northeast_gssurgo90m.tif')

mlramap <- sf::read_sf('C:/a/Ecological_Sites/GIS/Ecoregion/RegionalReview_296/CONUS.shp')
mlra99 <- subset(mlramap, MLRARSYM %in% '99')

extof99 <- ext(mlra99)
extof99a <- extof99
for(i in 1:4){
  extof99a[i] <- extof99a[i] + c(-100000, 100000, -100000, 100000)[i]}


listlmu <- s.reserve[s.reserve$coiid %in% missing99, 'lmapunitiid']


rast99 <- crop(soilgrid, extof99a)

rast99[rast99 %in% listlmu] <- 1
rast99[!rast99 %in% c(0,1)] <- -1
plot(rast99)
usergrouplmuid <- merge(s.reserve, usergroup[,c('grpname','dmuiid')], by='dmuiid')
groupsummary <- unique(usergrouplmuid[usergrouplmuid$coiid %in%  missing99,]$grpname)
usergrouplmuid <- subset(usergrouplmuid, )

rast99 <- crop(soilgrid, extof99a)
for (i in 1:length(groupsummary)){#i=3
  listlmu <- usergrouplmuid[usergrouplmuid$grpname %in% groupsummary[i] & usergrouplmuid$coiid %in% missing99,'lmapunitiid']
  rast99[rast99 %in% listlmu] <- i
};rast99[!rast99 %in% 1:length(groupsummary)] <- -1
plot(rast99)
writeRaster(rast99, 'fy2021-refresh/rast99.tif', overwrite=T)

mlrasummary <- unique(usergrouplmuid[usergrouplmuid$coiid %in%  missing99 & usergrouplmuid$ofbest >99,]$MLRARSYM)
rast99 <- crop(soilgrid, extof99a)
for (i in 1:length(mlrasummary)){#i=3
  listlmu <- usergrouplmuid[usergrouplmuid$MLRARSYM %in% mlrasummary[i] & usergrouplmuid$coiid %in% missing99 & usergrouplmuid$ofbest >99,'lmapunitiid']
  rast99[rast99 %in% listlmu] <- i
};rast99[!rast99 %in% 1:length(mlrasummary)] <- -1
plot(rast99)
writeRaster(rast99, 'fy2021-refresh/rast99a.tif', overwrite=T)


## what GRR/FLI needs to add ----
library(terra)
library(sf)
soilgrid <- rast('D:/GIS/SOIL/2021/northeast_gssurgo90m.tif')

mlramap <- sf::read_sf('C:/a/Ecological_Sites/GIS/Ecoregion/RegionalReview_296/CONUS.shp')
extof9799 <- subset(mlramap, MLRARSYM %in% c('97','98','99','94A','94C'))

extof9799 <- ext(extof9799)
extof9799a <- extof9799
for(i in 1:4){
  extof9799a[i] <- extof9799a[i] + c(-100000, 100000, -100000, 100000)[i]}




listlmu <- s.reserve[s.reserve$MLRARSYM %in% c('97','98','99','94A','94C') & s.reserve$ofbest ==100, 'lmapunitiid']



usergrouplmuid <- merge(s.reserve, usergroup[,c('grpname','dmuiid')], by='dmuiid')
groupsummary <- unique(usergrouplmuid[usergrouplmuid$MLRARSYM %in% c('97','98','99','94A','94C') & usergrouplmuid$ofbest ==100 & !usergrouplmuid$grpname %in% 
                                        c('12-GRR Grand Rapids, Michigan','12-GRR Projects','12-GRR MLRA MU/DMU',
                                          '12-FLI Projects','12-FLI Flint, Michigan'),]$grpname)


rast9799 <- crop(soilgrid, extof9799a)

for (i in 1:length(groupsummary)){#i=3
  listlmu <- usergrouplmuid[usergrouplmuid$MLRARSYM %in% c('97','98','99','94A','94C') & usergrouplmuid$ofbest ==100 & usergrouplmuid$grpname %in% groupsummary[i],'lmapunitiid']
  rast9799[rast9799 %in% listlmu] <- i
};rast9799[!rast9799 %in% 1:length(groupsummary)] <- -1
plot(rast9799)

writeRaster(rast9799, 'fy2021-refresh/rast9799.tif', overwrite=T)

listGRRnew <- subset(usergrouplmuid, MLRARSYM %in% c('97','98','99','94A','94C') & ofbest ==100 & grpname %in% groupsummary)
write.csv(listGRRnew, 'fy2021-refresh/listGRRnew.csv')

## what GRR/FLI needs to subtract ----
library(terra)
library(sf)
soilgrid <- rast('D:/GIS/SOIL/2021/northeast_gssurgo90m.tif')

mlramap <- sf::read_sf('C:/a/Ecological_Sites/GIS/Ecoregion/RegionalReview_296/CONUS.shp')
extof9799 <- subset(mlramap, MLRARSYM %in% c('97','98','99','94A','94C'))

extof9799 <- ext(extof9799)
extof9799a <- extof9799
for(i in 1:4){
  extof9799a[i] <- extof9799a[i] + c(-100000, 100000, -100000, 100000)[i]}




listlmu <- s.reserve[!s.reserve$MLRARSYM %in% c('97','98','99','94A','94C') & s.reserve$ofbest ==100, 'lmapunitiid']



usergrouplmuid <- merge(s.reserve, usergroup[,c('grpname','dmuiid')], by='dmuiid')
groupsummary <- unique(usergrouplmuid[!usergrouplmuid$MLRARSYM %in% c('97','98','99','94A','94C') & usergrouplmuid$ofbest ==100 & usergrouplmuid$grpname %in% 
                                        c('12-GRR Grand Rapids, Michigan','12-GRR Projects','12-GRR MLRA MU/DMU',
                                          '12-FLI Projects','12-FLI Flint, Michigan'),]$grpname)


rast9799 <- crop(soilgrid, extof9799a)

for (i in 1:length(groupsummary)){#i=3
  listlmu <- usergrouplmuid[!usergrouplmuid$MLRARSYM %in% c('96','97','98','99','94A','94C') & usergrouplmuid$ofbest ==100 & usergrouplmuid$grpname %in% groupsummary[i],'lmapunitiid']
  rast9799[rast9799 %in% listlmu] <- i
};rast9799[!rast9799 %in% 1:length(groupsummary)] <- -1
plot(rast9799)

writeRaster(rast9799, 'fy2021-refresh/rast9799a.tif', overwrite=T)

listGRRreject <- subset(usergrouplmuid, !MLRARSYM %in% c('96','97','98','99','94A','94C') & ofbest ==100 & grpname %in% groupsummary)
write.csv(listGRRreject, 'fy2021-refresh/listGRRreject.csv')




missing140 <- c(1384984, 539438, 539439, 539440, 539441, 539442, 539443, 539444, 539445, 539446, 539463, 539585, 539464, 539465, 539466, 539467, 539468, 539469, 539308, 539309, 539310, 1948970, 1948989, 539485, 296492, 296493, 296494, 296756, 539416, 539549, 539550, 539552, 539553, 539554, 539410, 539411, 539412, 539413, 1948979, 1948980, 539414, 539415, 539416, 539557, 539558, 539559, 539560, 1883685, 1883684, 299795, 299796, 299797, 539425, 294862, 539564, 539565, 539566, 296492, 296493, 296494, 539581, 539582, 539583, 539585, 539431, 539433, 539434)





s.comp <- subset(s.reserve, !is.na(s.reserve$taxorder), select=c(MLRARSYM, compname, Count))
s.comp$compname <- stringr::str_to_lower(s.comp$compname)
s.comp$compname <- stringr::str_to_title(s.comp$compname)
s.comp <- aggregate(list(Count=s.comp$Count), by=list(MLRARSYM=s.comp$MLRARSYM, compname=s.comp$compname), FUN='sum')

c <- aggregate(list(c=s.comp$Count), by=list(MLRARSYM=s.comp$MLRARSYM), FUN='sum')
MLRA.total <- subset(MLRA.total, total > 1000000)
s.comp <-  merge(s.comp, MLRA.total, by='MLRARSYM')
s.comp$percentmlra <- s.comp$Count/s.comp$total*100
comp.total <- aggregate(list(ctotal=s.comp$percentmlra), by=list(compname=s.comp$compname), FUN='sum')
s.comp <-  merge(s.comp, comp.total, by='compname')
s.comp$affinity <- round(s.comp$percentmlra/s.comp$ctotal *100, 2)

write.csv(s.comp, 'fy2021-refresh/s.comp.affinity.csv', row.names = F)

# missing 139 ----

missing139 <- c(251531, 2238292, 2205085, 779586, 1393103, 322943, 776812, 244596, 2743334, 324700, 244608, 2746156, 774529, 788245, 788245, 788248, 324304, 322974, 367609, 324312, 788165, 816450, 816458, 1328682, 249237, 249238, 249239, 249240, 249241, 249242, 249243, 776576, 244680, 244683, 244686, 244691, 245094, 244696, 244699, 324516, 249246, 2189931, 2163155, 324517, 324518, 2238289, 2238289, 245100, 324518, 245100, 776581, 324519, 324520, 324521, 324522, 245101, 324522, 249260, 2383530, 2383534, 2414568, 245105, 2192932, 2173937, 2192933, 890479, 2383580, 2383580, 2383580, 791576, 2204950, 2204950, 2484849, 2305132, 2305132, 245113, 245113, 895392, 895401, 1239181, 245114, 324547, 324548, 245116, 245116, 249287, 249286, 249287, 816305, 788227, 816305, 251629, 788166, 2379293, 1328675, 2383564, 2383567, 2383553, 244763, 244764, 244765, 244766  )
usergrouplmuid <- merge(s, usergroup[,c('grpname','dmuiid')], by='dmuiid')
missing139 <-  subset(usergrouplmuid, coiid %in% missing139)
missing139 <- subset(missing139, ofbest == 100)
write.csv(missing139, 'fy2021-refresh/missing139.csv', row.names = F)
#mostly great lakes marsh ----
orphans <- data.frame(cbind(
  lmukey = c(192629,187357,192688, 186374, 193342, 189435, 186776, 187220, 3015048, 169588, 2633031, 187401),
esiid = c('R097XA024MI','R096XY002MI','R097XA024MI', 'R098XA002MI','R096XY002MI','R096XY002MI','F099XY010MI','F099XY010MI','F099XY010MI','F099XY010MI', 'F139XY010OH', 'F099XY010MI')

))

orphans <- merge(orphans, s.reserve[s.reserve$majcompflag %in% 'TRUE',c('lmapunitiid','dmuiid', 'muiid', 'coiid', 'compname')], by.x='lmukey', by.y = 'lmapunitiid')


