# Kralovo pole grasslands - analysis of transect data ----
# David Zeleny, zeleny@ntu.edu.tw ----

library (vegan)
library (cluster)
library (ade4)
library (adegraphics)
library (adespatial)
library (vegetarian)

## Read data ----
setwd ('c:/Users/Zeleny/Dropbox/experimenty/KrPole train station/vegetacni snimky')

com <- read.delim ('KrPole station releves.txt', row.names = 1)
spe <- read.delim ('KrPole station species.txt')
env <- read.delim ('KrPole station env.txt')

rownames (com) <- 1:21  # replace original "releve_nr" by order along the platform transect
com.log <- log1p (com)


## Are data homogeneous or hetergeonous? ----
decorana (com.log)  # length of DCA1 = 1.8 SD - homogeneous data
ordiplot (decorana (com.log), display = 'si', type = 't') # just try if I can visualize

## NMDS analysis on species composition data and randomized transect ----
set.seed (1254468)
nmds <- metaMDS (com.log)
com.log.rand <- t(apply (com.log, MARGIN = 1, sample))
nmds.rand <- metaMDS (com.log.rand)

jpeg (filename = 'KrPole_NMDS.jpg', width = 12, height = 6, units = 'cm', res = 300, pointsize = 6)
par (mfrow = c(1,2))
ordiplot (nmds, display = 'si', type = 'n', main = 'Skutečný transekt / Real transect')
scores.nmds <- scores (nmds, display = 'si', choices = 1:2)
for (i in seq (1,20))
  lines (scores.nmds[i:(i+1),], col = 'gray', lwd = 2)
ordilabel (nmds)
par (xpd = TRUE)
mtext("A", side = 3, line = 1.3, adj = -0.2, cex = 1.4, font = 2) 

ordiplot (nmds.rand, display = 'si', type = 'n', main = 'Permutovaná data / Permuted data')
scores.nmds.rand <- scores (nmds.rand, display = 'si', choices = 1:2)
for (i in seq (1,20))
  lines (scores.nmds.rand[i:(i+1),], col = 'gray', lwd = 2)
ordilabel (nmds.rand)
mtext("B", side = 3, line = 1.3, adj = -0.2, cex = 1.4, font = 2) 
dev.off ()

## RDA with measured environmental variables ----
set.seed (123434)
rda_depth_mean <- rda (com.log ~ soil_depth_mean_cm, data = env) 
anova (rda_depth_mean)   # P = 0.069, marginally significant
RsquareAdj (rda_depth_mean)  # adjR2 = 0.020

set.seed (13441)
rda_depth_sd <- rda (com.log ~ soil_depth_sd_cm, data = env)
anova (rda_depth_sd) # P = 0.079, marginally significant
RsquareAdj (rda_depth_sd) # adjR2 = 0.019 

set.seed (1324134)
rda_ph <- rda (com.log ~ pH, data = env)
anova (rda_ph)  # P = 0.008
RsquareAdj (rda_ph) # adjR2 = 0.036

set.seed (98797345)
rda_cond <- rda (com.log ~ cond, data = env)
anova (rda_cond)  # not sign
RsquareAdj (rda_cond) # adjR2 = 0.003

## Forward selection on environmental variables ----
set.seed (1273548)
FS <- forward.sel (Y = com.log, X = env[,7:10], adjR2thresh = RsquareAdj (rda (com.log ~ as.matrix (env[,7:10])))$adj, alpha = .05)  
anova (rda (com.log ~ pH, data = env))

## dbMEM and forward selection ----
pcnm <- dbmem (env$trans_position, silent = FALSE)
pcnm.axes <- as.matrix (pcnm)
fw.pcnm <- forward.sel (com.log, pcnm.axes, adjR2thresh = RsquareAdj (rda (com.log ~ pcnm.axes))$adj, alpha = .05) #MEM1-MEM5

## Variation partitioning (env vs space) ----
VP <- varpart (com.log, env[,'pH'], pcnm.axes[,fw.pcnm$order])

set.seed (83418324)
space.env.rda <- rda (com.log, cbind (env$pH, pcnm.axes[,fw.pcnm$order]))
anova (space.env.rda)  # P < 0.001
RsquareAdj (space.env.rda) # adjR2 = 0.249

set.seed (341957232)
env.parc.space <- rda (com.log, env[,'pH'], pcnm.axes[,fw.pcnm$order])
anova (env.parc.space)  # P = 0.068 (marginally significant)
RsquareAdj (pcnm.parc.space)  # adjR2 = 0.024

set.seed (34132487)
space.parc.env <- rda (com.log, pcnm.axes[,fw.pcnm$order], env[,'pH'])  
anova (space.parc.env) # P < 0.001
RsquareAdj (space.parc.env) # adjR2 = 0.213

set.seed (42542435)
env.rda <- rda (com.log, env[,'pH'])
anova (env.rda)  # P < 0.01
RsquareAdj (env.rda) # adjR2 = 0.036

set.seed (4238935)
space.rda <- rda (com.log ~ pcnm.axes[,fw.pcnm$order])
anova (space.rda) # P < 0.001
RsquareAdj (space.rda) # adjR2 = 0.226

jpeg (filename = 'KrPole_varpart.jpg', width = 12, height = 9.5, units = 'cm', res = 300, pointsize = 6)
plot (VP, Xnames = c('Prostředí\n/ Env',"Prostor\n/ Space"), digits = 2, cex = 1.5, id.size = 1.6)
dev.off ()

## Species richness ----
rich <- apply (species, MARGIN = 1, FUN = d, q = 0)
mean (rich)  # 29.8
median (rich) # 30
range (rich)  # 24-39
