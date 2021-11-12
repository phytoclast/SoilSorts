library(soilDB)
library(aqp)
library(plyr)
library(foreign)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#establish dominant MRLA per lmukey ----
mlra <- read.csv('fy2021-refresh/mlragrid.csv')
mlra.mu <- read.csv('fy2021-refresh/combgrid.csv')
mlra <- subset(mlra, Count > 100000)
mu <- readRDS(file='fy2021-refresh/mu.RDS')
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


# mu <- get_component_correlation_data_from_NASIS_db(SS=F)
# saveRDS(mu, file='fy2021-refresh/mu.RDS')
mu <- readRDS(file='fy2021-refresh/mu.RDS')

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
spodic <- read.delim("data/Spodics.txt")


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
#flood
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



rm(compsorts, fc.hz.bedrock,  fc.hz.bedrock.min, fc.hz.carb, fc.hz50.total, fc.hz.carb.min, fc.hz.notnull, fc.hz.sum, fc.hz150.total, cm.flood, fc.hz.humic.max, fc.hz.humic, fc.hz.mineral )


summarydrcl <- aggregate(s[,c('drainagecl')], by=list(s$drainagecl, s$hydricrating,s$Water_Table), FUN = 'length')

s <- merge(mlra.mu, s, by.x = 'ncnesoil30m', by.y = 'lmapunitiid')
s.lmu <- unique(subset(s, select =c(muiid, ncnesoil30m)))
s <- unique(subset(s, select = -c(ncnesoil30m, Value, Count, percentmlra, MLRA.value, MLRA.total, mu.total)))
Appalachia <- c("124", "126", "127", "130A", "136", "139A", "139B", "140", "141", "142", "143", "144A", "144B", "145", "146", "147", "148", "149A", "149B", "153C", "153D")
Plains <- c("101", "108", "110", "111A", "111B", "111C", "111D", "111E", "114A", "90A", "93B","94AA", "94AB", "94B",  "94C",  "95B",  "96A",  "96B",  "97A", "97B", "98A1", "98A2", "98B", "99A", "99B")
Northern <- c("94AA", "94AB", "94B",  "94C",  "95B",  "96A",  "96B")
Frigid <- c("94AA", "94AB", "94B",  "94C",  "95B")


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
                      ifelse((s$Water_Table<0 & s$T150_OM >=20)|((grepl("Histic",s$taxsubgrp)|grepl("Histosols",s$compname)|grepl("ists",s$taxsubgrp)|grepl("ists",s$taxclname))& s$Water_Table<=0),"Muck",
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

s <- merge(s.lmu, s, by='muiid')

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
s.111A.missing <- subset(s.111A, ofbest >= 50 & majcompflag %in% 'TRUE' & MLRARSYM %in% '111A' & is.na(ES))

s.111A.select <- unique(subset(s.111A, ofbest >= 50 & majcompflag %in% 'TRUE' & MLRARSYM %in% '111A', select= c(ES, MLRARSYM, Series, coiid, muname)))
write.csv(s.111A.select, 'output/s.111A.select.csv', row.names = F)

#111B ----
#

ES111B.list <- list(
  F111BY101IN = c('Toledo', 'Sebring', 'Patton', 'Minster', 'Milford', 'Luray', 'Lenawee', 'Latty', 'Hoytville', 'Bono'),
  F111BY102IN = c('Toledo', 'Sebring', 'Patton', 'Minster', 'Milford', 'Luray', 'Lenawee', 'Latty', 'Hoytville', 'Bono'),
  F111BY201IN = c('Tice', 'Sloan', 'Saranac', 'Harrod', 'Glendora', 'Cohoctah', 'Ceresco'),
  F111BY202IN = c('Spencerville', 'Rossburg', 'Ross', 'Medway', 'Landes', 'Armiesburg'),
  F111BY203IN = c('Washtenaw', 'Wabasha', 'Shoals', 'Orrville', 'Newark', 'Defiance', 'Coesse', 'Algiers'),
  F111BY204IN = c('Tioga', 'Stonelick', 'Nolin', 'Lobdell', 'Lindside', 'Knoxdale', 'Gessie', 'Genesee', 'Flatrock', 'Eel', 'Chagrin'),
  F111BY302IN = c('Randolph'),
  F111BY303IN = c('Weikert', 'Ritchey', 'Newglarus', 'Milton', 'Latham', 'Heverlo', 'Eleva', 'Castalia', 'Brecksville', 'Biglick', 'Berks'),
  F111BY403IN = c('Whitaker', 'Stone', 'Sleeth', 'Riverdale', 'Matherton', 'Jimtown', 'Homer', 'Digby', 'Bogart', 'Aptakisic'),
  F111BY404IN = c('Thackery', 'Sisson', 'Scioto', 'Rush', 'Riverdale', 'Oshtemo', 'Ormas', 'Ockley', 'Muncie', 'Martinsville', 'Leoni', 'Kosciusko', 'Kalamazoo', 'Ionia', 'Haney', 'Gallman', 'Fox', 'Eldean', 'Chili', 'Casco', 'Bronson', 'Boyer', 'Belmore'),
  F111BY501IN = c('Wolcott', 'Wetzel', 'Pewamo', 'Paulding', 'Pandora', 'Mermill', 'Kokomo', 'Brookston', 'Berville', 'Barry', 'Alvada'),
  F111BY502IN = c('Tiro', 'Rimer', 'Nappanee', 'Metamora', 'Markton', 'Macomb', 'Locke', 'Lamberjack', 'Haskins', 'Elliott', 'Crosier', 'Conover', 'Capac', 'Blount', 'Aubbeenaubbee'),
  F111BY503IN = c('Williamstown', 'Wawasee', 'Strawn', 'Steamburg', 'Shawtown', 'Seward', 'Russell', 'Riddles', 'Rawson', 'Parr', 'Pacer', 'Owosso', 'Mortimer', 'Morley', 'Mississinewa', 'Miami', 'Metea', 'Marlette', 'Lykens', 'Lybrand', 'Kidder', 'Kendallville', 'Houcktown', 'Hillsdale', 'Hennepin', 'Glynwood', 'Cadmus', 'Amanda'),
  R111BY001IN = c('Roundhead', 'Palms', 'Olentangy', 'Ogden', 'McGuffey', 'Linwood', 'Kerston', 'Adrian'),
  R111BY002IN = c('Warners', 'Rollin', 'Muskego', 'Martisco', 'Edwards'),
  R111BY003IN = c('Napoleon', 'Houghton', 'Carlisle', 'Boots'),
  R111BY301IN = c('Millsdale', 'Joliet', 'Dunbridge'),
  R111BY401IN = c('Westland', 'Sebewa', 'Rensselaer', 'Pella', 'Millgrove', 'Lippincott', 'Gilford', 'Colwood'),
  R111BY402IN = c('Wasepi', 'Warsaw', 'Rodman', 'Nineveh', 'Kane', 'Crane')
 )
for(i in 1:length(ES111B.list)){
  ES111B.df0 <- as.data.frame(merge(names(ES111B.list[i]),ES111B.list[i])); colnames(ES111B.df0) <- c('ES','Series')
  if(i == 1){ES111B <- ES111B.df0}else{ES111B <- rbind(ES111B.df0, ES111B)}};rm(ES111B.df0) 

s.111B <- merge(ES111B, s,  by.x = 'Series', by.y = 'compname', all.y = T)
s.111B.missing <- subset(s.111B, ofbest >= 50 & majcompflag %in% 'TRUE' & MLRARSYM %in% '111B' & is.na(ES))

s.111B.select <- unique(subset(s.111B, ofbest >= 50 & majcompflag %in% 'TRUE' & MLRARSYM %in% '111B', select= c(ES, MLRARSYM, Series, coiid, muname)))
write.csv(s.111B.select, 'output/s.111B.select.csv', row.names = F)


#### ---- exports2 
muck <- subset(s, ofbest >= 99 & majcompflag %in% 'TRUE' & Site %in% 'M2 Mucky Depression' )
m98 <- unique(muck[muck$LRU %in% c('98A1', '98A2') ,]$ncnesoil30m)
m98B <- unique(muck[muck$LRU %in% c('98B') ,]$ncnesoil30m)
m99 <- unique(muck[muck$LRU %in% c('99A', '99B') ,]$ncnesoil30m)
m97 <- unique(muck[muck$LRU %in% c('97A') ,]$ncnesoil30m)
m97B <- unique(muck[muck$LRU %in% c('97B') ,]$ncnesoil30m)
m96 <- unique(muck[muck$LRU %in% c('96A','96B') ,]$ncnesoil30m)
m94 <- unique(muck[muck$LRU %in% c('94AA', '94AB', '94C') ,]$ncnesoil30m)
m111 <- unique(muck[muck$LRU %in% c('111', '111C', '111B') ,]$ncnesoil30m)
write.dbf(muck, "output/s.dbf") 
library(terra)
soilrast <- rast('C:/a/Ecological_Sites/GIS/Soil/ncnesoil30m.tif')
soil90 <- aggregate(soilrast, fact=3, fun='modal', na.rm=T)

soil90[soil90 %in% m97] <- 970
soil90[soil90 %in% m97B] <- 971
soil90[soil90 %in% m99] <- 990
soil90[soil90 %in% m98B] <- 981
soil90[soil90 %in% m98] <- 980
soil90[soil90 %in% m96] <- 960
soil90[soil90 %in% m94] <- 940
soil90[soil90 %in% m111] <- 1110
soil90[!soil90 %in% c(940,960,970,971,980,981,990,1110)] <- 0
writeRaster(soil90, 'output/mucks.tif', overwrite=T)


muck98 <- subset(s, ofbest >= 99 & LRU %in% c('98A1','98A2') & Site %in% 'M2 Mucky Depression' )
write.csv(muck98, 'output/muck98.csv', na='', row.names = F)
muck97 <- subset(s, ofbest >= 99 & LRU %in% c('97A') & Site %in% 'M2 Mucky Depression' )
write.csv(muck97, 'output/muck97.csv', na='', row.names = F)
muck99 <- subset(s, ofbest >= 99 & LRU %in% c('99A','99B') & Site %in% 'M2 Mucky Depression' )
write.csv(muck99, 'output/muck99.csv', na='', row.names = F)
#### ---- exports
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




