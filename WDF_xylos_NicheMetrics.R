#########################################################################################################
## Calculate SEA, TNW and Rel. Niche Index for xylos  ###################################################
## TNW, SEAc, RNI -> SIBER                     ##########################################################
## Quantifies trophic niche for each species (SEAc) and TNW for community (log) #########################
## Use RINI function to calculate species specialization per log  #######################################
## By: Philip J. Manlick, 06/14/2021 FINAL  #############################################################
#########################################################################################################

rm(list=ls())

library(dplyr)
library(tidyverse)
library(SIBER)
library(viridis)
library(performance)
library(sjPlot)

setwd("/Volumes/GoogleDrive/My Drive/Postdoc_UNM/Woodfall/analyses/FINAL")

data = read.csv("wdf_xylos.csv", strip.white = T, header = T)
data$species<-as.factor(data$species)
data$log<-as.factor(data$log)
data$type<-as.factor(data$type) 
data$density<-as.factor(data$density)

#add count data -> # of iso samples per species per log
counts = data %>% group_by(species, log) %>% count()

#basic plots

#by log
data %>% ggplot(aes(x=d13C, y=d15N)) + geom_point(aes(colour=species)) + facet_wrap(~log) + theme_bw() +
  scale_color_viridis(discrete = T) + theme(legend.position = "top")
#by species
data %>% ggplot(aes(x=d13C, y=d15N)) + geom_point(aes(colour=log)) + facet_wrap(~species) + theme_bw() +
  scale_color_viridis(discrete = T) + theme(legend.position = "top")

# #drop xylo-1-15 (LW34) and xylo-3-15 (LG16)
data <- data[-c(130, 479),]

#rename d13C and d15N to "iso1" & "iso2" for SIBER, species -> "group"
data<-rename(data, c("group"="species", "iso1"="d13C", "iso2"="d15N"))

#split df into logs
logs = split(data, data$log)

summary(logs)

#create df for each log
lg10 <- logs$`L-G-10`
lg16 <- logs$`L-G-16`
lg34 <- logs$`L-G-34`
lg7 <- logs$`L-G-7`
lw22 <- logs$`L-W-22`
lw25 <- logs$`L-W-25`
lw28 <- logs$`L-W-28`
lw34 <- logs$`L-W-34`
ly29 <- logs$`L-Y-29`
ly31 <- logs$`L-Y-31`

#############################################################################################################
#############################################################################################################
### siberKAPOW -> calculate Niche Index (Shepherd et al 2018, Ecol Letters)
### Calculates 95% SEAc for each species (i.e., isotopic niche)
### Takes union of all SEAc to calculate TNW of population (i.e., log)
### Calculate "Niche Index" -> % of TNW occupied per species (specialists vs generalists)
#############################################################################################################
#############################################################################################################

## Call functions 
source("Niche_functions.R") 

### ----------------------------------------------------------------------###
### plot niche ellipses for each log -> black = TNW, colors = species
### ----------------------------------------------------------------------###

kapow2(Y=lg10) #5spp
kapow2(Y=lg16) #5spp
kapow2(Y=lg34) #6spp
kapow2(Y=lg7) #6spp
kapow2(Y=lw22) #3spp
kapow2(Y=lw25) #3spp
kapow2(Y=lw28) #4spp
kapow2(Y=lw34) #4spp
kapow2(Y=ly29) #1spp
kapow2(Y=ly31) #1spp

####--------------------------------------------------------------------####
#### Run siberKapow for each log -> caclulates ellipses and boundary 
####--------------------------------------------------------------------####

lg10.kapow <- siberKapow(lg10, isoNames = c("iso1", "iso2"), group = "group", pEll = 0.95)
lg16.kapow <- siberKapow(lg16, isoNames = c("iso1", "iso2"), group = "group", pEll = 0.95)
lg34.kapow <- siberKapow(lg34, isoNames = c("iso1", "iso2"), group = "group", pEll = 0.95)
lg7.kapow <- siberKapow(lg7, isoNames = c("iso1", "iso2"), group = "group", pEll = 0.95)
lw22.kapow <- siberKapow(lw22, isoNames = c("iso1", "iso2"), group = "group", pEll = 0.95)
lw25.kapow <- siberKapow(lw25, isoNames = c("iso1", "iso2"), group = "group", pEll = 0.95)
lw28.kapow <- siberKapow(lw28, isoNames = c("iso1", "iso2"), group = "group", pEll = 0.95)
lw34.kapow <- siberKapow(lw34, isoNames = c("iso1", "iso2"), group = "group", pEll = 0.95)
ly29.kapow <- siberKapow(ly29, isoNames = c("iso1", "iso2"), group = "group", pEll = 0.95)
ly31.kapow <- siberKapow(ly31, isoNames = c("iso1", "iso2"), group = "group", pEll = 0.95)

####--------------------------------------------------------------------####
# Calculate SEAc, TNW, Niche Index -
####--------------------------------------------------------------------####

lg10.rini <- rini(lg10.kapow)
lg10.rini$log <- 'L-G-10' #add log number back to df

lg16.rini <- rini(lg16.kapow)
lg16.rini$log <- 'L-G-16'

lg34.rini <- rini(lg34.kapow)
lg34.rini$log <- 'L-G-34'

lg7.rini <- rini(lg7.kapow)
lg7.rini$log <- 'L-G-7'

lw22.rini <- rini(lw22.kapow)
lw22.rini$log <- 'L-W-22'

lw25.rini <- rini(lw25.kapow)
lw25.rini$log <- 'L-W-25'

lw28.rini <- rini(lw28.kapow)
lw28.rini$log <- 'L-W-28'

lw34.rini <- rini(lw34.kapow)
lw34.rini$log <- 'L-W-34'

ly29.rini <- rini(ly29.kapow)
ly29.rini$log <- 'L-Y-29'

ly31.rini <- rini(ly31.kapow)
ly31.rini$log <- 'L-Y-31'

# bind dfs
niche <- rbind(lg10.rini, lg16.rini, lg34.rini, lg7.rini,
               lw22.rini, lw25.rini, lw28.rini,
               lw34.rini, ly29.rini, ly31.rini)

# summarise orginal df to get richness/samples/etc per community
comm = data %>% group_by(group, log, density) %>% 
  summarise(richness=mean(richness), sampled=mean(sampled), type=names(which.max(table(type))))

#rename group, create merge columns
comm <-rename(comm, c("Species"="group"))
comm$merge <- as.factor(paste(comm$Species, comm$log, sep="."))
niche$merge <- as.factor(paste(niche$Species, niche$log, sep="."))

#merge dfs
dat = merge(niche, comm)
dat$prop <- dat$sampled/dat$richness #proportion of community sampled

write.csv(dat, "NicheStats_samplerichness_noutliers.csv")
save.image("/Volumes/GoogleDrive/My Drive/Postdoc_UNM/Woodfall/analyses/FINAL/Xylos_niches_dropoutliers_DATA.RData")

###########################################################################
### Niche metrics ~ species richness
###########################################################################

### ----------------------------------------------------------- ###
### TNWlog -> total isotopic niche width per log
### ----------------------------------------------------------- ###

#filter out logs with less than half of community sampled to test TNW ~ richness
#also drop L-G-31 -> unknown richness
dat.drop = dat %>% group_by(log, density, TNW, prop) %>%
  summarise(richness=mean(richness), sampled=mean(sampled), type=names(which.max(table(type))))

dat.drop$type<-as.factor(dat.drop$type)
dat.drop$density<-as.factor(dat.drop$density)

#plot TNW ~ sampled richness
dat.drop %>%  ggplot(aes(x=sampled, y=TNW)) + 
  geom_point(aes(fill=log), size=4, shape=21) + geom_smooth(method='lm') + theme_bw()

#test effect of log density -> hard vs soft
dat.drop %>% ggplot(aes(x=sampled, y=TNW)) + 
  geom_point(aes(fill=log), size=4, shape=21) + geom_smooth(method='lm') + theme_bw() +
  facet_wrap(~density, scales = "free") 

#TNW ~ prop
dat.drop %>%  ggplot(aes(x=prop, y=TNW)) + 
  geom_point(aes(fill=log), size=4, shape=21) + geom_smooth(method='lm') + theme_bw()

## Linear models -------------------------------------------- ####

# Richness tests -> number of species sampled for isotopes
tnw1 = lm(TNW ~ sampled, dat.drop) 
anova(tnw1) #Significant effect of richness
summary(tnw1) #r2 = 0.4811, p = 0.0261*
check_model(tnw1) #decent fit
plot_model(tnw1, type = "pred", show.data = T)
plot(tnw1)

### --------------------------------------------------------------------------- ###
### 95% SEAsplog stats - isotopic niche width per species per log 
### --------------------------------------------------------------------------- ###

#merge count data to account for sample size diffs
counts$merge<- as.factor(paste(counts$species, counts$log, sep="."))
dat.drop2<-merge(dat,counts)
dat.drop2<- dat.drop2 %>% select(-species)

write.csv(dat.drop2, "NicheStats_dropoutliers_withcounts.csv")

#define factors
dat.drop2$density <- as.factor(dat.drop2$density)
dat.drop2$type <- as.factor(dat.drop2$type)

#plot -  SEA ~ sampled richness
dat.drop2 %>% ggplot(aes(x=sampled, y=SEAc)) + 
  geom_point(aes(fill=Species), size=4, shape=21) + geom_smooth(method='lm') + theme_bw()

#plot - facet by density
dat.drop2 %>% ggplot(aes(x=sampled, y=SEAc)) + 
  geom_point(aes(fill=Species), size=4, shape=21) + geom_smooth(method='lm') + theme_bw() + 
  facet_wrap(~density, scales = "free") 

#plot - facet by species **some variaiton in slope by spp
dat.drop2 %>% ggplot(aes(x=sampled, y=SEAc)) + 
  geom_point(aes(fill=Species), size=4, shape=21) + geom_smooth(method='lm') + theme_bw() + 
  facet_wrap(~Species, scales = "free") 

#plot - by sample size - little relation
dat.drop2 %>% ggplot(aes(x=n, y=SEAc)) + 
  geom_point(aes(fill=Species), size=4, shape=21) + geom_smooth(method='lm') + theme_bw()

## Linear models------------------------------------------------------------ ####

#basic model
sea.1 = lm(SEAc ~ sampled, dat.drop2) 
anova(sea.1) #no significant effect of richness per se on species niche width
summary(sea.1)
check_model(sea.1) #good fit
plot(sea.1) #good fit
plot_model(sea.1, type = "pred", show.data = T)

## main model - account for species and sample sizes
sea.2 = lm(SEAc ~ sampled + Species + n, dat.drop2) 
anova(sea.2)
summary(sea.2)
check_model(sea.2) #decent fit

### -------------------------------------------------------------- ###
### Niche index stats - test niche packing hypothesis
### Niche.Index = SEAc/TNW, proportion of niche occupied per spp
### -------------------------------------------------------------- ###

#plot -  RNI ~ sampled richness
dat.drop2 %>% ggplot(aes(x=sampled, y=Niche.Index)) + 
  geom_point(aes(fill=Species), size=4, shape=21) + geom_smooth(method='lm') + theme_bw()

#plot - facet by species
dat.drop2 %>% ggplot(aes(x=sampled, y=Niche.Index)) + 
  geom_point(aes(fill=Species), size=4, shape=21) + geom_smooth(method='lm') + theme_bw() + 
  facet_wrap(~Species, scales = "free") + theme(legend.position = "none")

## Linear models ------------------------------------------------------- ####

#basic model
rni1 = lm(Niche.Index ~ sampled, dat.drop2) 
anova(rni1) #very significant effect of richness on proportion of niche occupied/spp
summary(rni1) 
check_model(rni1) #good fit

# MAIN MODEL REPORTED - accounts for sampling differences
rni2 = lm(Niche.Index ~ sampled + Species + n, dat.drop2) 
anova(rni2) #very significant effect of richness on proportion of niche occupied/spp
summary(rni2) # spp not significant, but betas very diff for some spp (see plot above)
check_model(rni2) #good fit

save.image("/Volumes/GoogleDrive/My Drive/Postdoc_UNM/Woodfall/analyses/Xylos_niches_dropoutliers_DATA.RData")


#########################################################################################################
## Figures -> 3 panel with SEA, TNW, and RNI
#########################################################################################################

library(ggthemes)
library(ggpubr)
library(ggeffects)

#A - SEAc
sea <- ggpredict(sea.2, terms ="sampled")
a=ggplot(sea, aes(x=x, y=predicted)) + theme_classic() + theme(legend.position = "none") +
  geom_point(data=dat.drop2, aes(x=sampled, y=SEAc), 
             shape=21, size = 5, stroke = .75, alpha=.75, col ="grey40", fill="grey70") + 
  labs(x="Xylophagid Species Richness", y=expression("Standard Ellipse Area (SEA; \u2030\u00b2)")) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text=element_text(size=12, color='black'),axis.title=element_text(size=14))

#B - SEAc
tnw <- ggpredict(tnw1, terms ="sampled")
b=ggplot(tnw, aes(x=x, y=predicted)) + theme_classic() + theme(legend.position = "none") +
  geom_point(data=dat.drop, aes(x=sampled, y=TNW), 
             shape=21, size = 5, stroke = .75, alpha=.75, col ="grey40", fill="grey70") + 
  labs(x="Xylophagid Species Richness", y=expression("Total Niche Width (TNW; \u2030\u00b2)")) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text=element_text(size=12, color='black'),axis.title=element_text(size=14)) + 
  geom_line(lwd=1.5, color="grey70")

#C - RNI
rni <- ggpredict(rni2, terms ="sampled")
c=ggplot(rni, aes(x=x, y=predicted)) + theme_classic() + theme(legend.position = "none") +
  geom_point(data=dat.drop2, aes(x=sampled, y=Niche.Index), 
             shape=21, size = 5, stroke = .75, alpha=.75, col ="grey40", fill="grey70") + 
  labs(x="Xylophagid Species Richness", y=expression("Relative Niche Index (RNI)")) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text=element_text(size=12, color='black'),axis.title=element_text(size=14)) + 
  geom_line(lwd=1.5, color="grey70")

ggarrange(a, b, c,ncol = 3, nrow = 1, 
          font.label = list(family = "Helvetica"))

#########################################################################################################
## (1) Data show significant increase in TNW with richness as predicted but the
## relationship is mediated by log type (L/G/Y) and density (hard/soft wood). HOWEVER - those models
## have very bad fit and don't meet assumptions of lm. Tried a lmm but "is singular" indicating
## models incorporating wood type/density are overfit. I do think there is a difference between those groups, 
## just not enough data to test it. 
## (2) SEAc (species niche width) did not respond significantly to richness, species, type, sample size. 
## Appears relatively fixed which suggests isotopic/trophic niches are relatively static and don't respond to 
## competition as would be predicted by niche packing. Rather, they appear complementary.
## (3) Relative Niche Index (RNI=SEAc/TNW) showed significant response to richness at community and species
## levels, as well interaction between species*richness. 
## IN TOTAL: This suggests niche complementarity among xylos at community-level, where additional 
## species increase TNW, but species don't significantly increase or decrease their trophic niches. 
## Niche index interaction suggests that some species respond to richness via specialization (smaller NI) 
## while others generalize (bigger NI). Would be interesting to tie response to morphology OR 
## genetic relatedness to see if this trait is conserved in any way across genera. 
#########################################################################################################

data 

# check wood against consumer values
wood = read.csv("wdf_wood.csv", strip.white = T, header = T)
wood = wood %>% filter(Species %in% c("Pinus elliotti", "Celtis laevigata", "Quercus rubra")) %>% 
  select(Sample.ID, Location, Species, name, type, d13C, d15N)

#plot softwoods - pine
wood.pe <- wood %>% filter(Species == "Pinus elliotti")
pe.xy <- data %>% filter(type == "Pinus elliotti")
pe.xy %>% ggplot(aes(x=d13C, y=d15N)) + geom_point(aes(color=species)) + 
  geom_point(aes(x=d13C, y=d15N), data=wood.pe) + facet_wrap(~log)

#plot hardwoods - sugarberry
wood.cl <- wood %>% filter(Species == "Celtis laevigata")
cl.xy <- data %>% filter(type == "Celtis laevigata")
cl.xy %>% ggplot(aes(x=d13C, y=d15N)) + geom_point(aes(color=species)) + 
  geom_point(aes(x=d13C, y=d15N), data=wood.cl) + facet_wrap(~log)

#plot hardwoods - oak
wood.qr <- wood %>% filter(Species == "Quercus rubra")
qr.xy <- data %>% filter(type == "Quercus rubra")
qr.xy %>% ggplot(aes(x=d13C, y=d15N)) + geom_point(aes(color=species)) + 
  geom_point(aes(x=d13C, y=d15N), data=wood.qr) + facet_wrap(~log)

#calculate ranges
data.range=data %>% group_by(log, type) %>% summarise(crange=max(d13C)-min(d13C), nrange=max(d15N)-min(d15N))
wood.range=wood %>% group_by(Species) %>% summarise(wcrange=max(d13C)-min(d13C), wnrange=max(d15N)-min(d15N))
wood.range= wood.range %>% rename("type"="Species")

range=merge(data.range, wood.range)
range$Cdiff<-range$crange-range$wcrange
range$Ndiff<-range$nrange-range$wnrange
range$Cdiff<-format(round(range$Cdiff, 2), nsmall = 2)
range$Ndiff<-format(round(range$Ndiff, 2), nsmall = 2)
range

save.image("/Volumes/GoogleDrive/My Drive/Postdoc_UNM/Woodfall/analyses/Xylos_niches_DATA.RData")

