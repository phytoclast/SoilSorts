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

mlra.mu <-  unique(subset(mlra.mu, ofbest > 99, select =c(Count, affinity, ofbest, dmuiid, MLRA)))
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
lru.mu <- subset(lru.mu, lruaffinity/best >= 0, select=c(lruaffinity, dmuiid, MLRA, LRU))

mlra.mu <- merge(mlra.mu, lru.mu, by=c('dmuiid', 'MLRA'), all.x = T)

missing <- c(272449, 160652, 272463)

######################################----
#usergroup <- get_group_from_table("datamapunit")
#saveRDS(usergroup, 'fy2021-refresh/usergroup20220119.RDS')
usergroup <- readRDS("fy2021-refresh/usergroup20220119.RDS")
# comp <- get_component_data_from_NASIS_db(SS=F)
# saveRDS(comp, file='fy2021-refresh/comp.RDS')
comp <- readRDS("fy2021-refresh/comp.RDS")
mu <- readRDS("fy2021-refresh/mu.RDS")

s <- merge(mlra.mu, mu[,c('lmapunitiid','dmuiid','nationalmusym', 'muname')], by='dmuiid')
s <- merge(usergroup[,c('dmuiid','grpname')], s, by='dmuiid')

s.missing <- subset(s, lmapunitiid %in% missing)