#Birds and Sound Indices 86 plots R4.4.1
rm(list=ls(all=TRUE))

#- check for packages and if necessary install libraries 
libraries <- c("SpadeR","ecodist","vegan","snowfall")       # libraries required for the session
inst      <- libraries %in% installed.packages()         # Install CRAN packages (if not already installed)
if(length(libraries[!inst]) > 0) install.packages(libraries[!inst])
lapply(libraries, require, character.only=TRUE)     # Load packages into session
rm(list=ls())

library(SpadeR)
library(ecodist)
library(vegan)
library(snowfall)
library(effects)

source("div_pair.R")

data_org<- read.csv2("BirdCommunityNamePhylogeny_Indices.csv", row.names=1,header=T)
rownames(data_org)
names(data_org)
birds<-data_org[,14:347]
tab<-birds
dim(tab)

tab[is.na(tab)] <- 0

# hier bei cpus die anzahl deiner rechenkerne eingeben
tt <- div_pair(com = tab, nboot = 3, parallel = T, cpus = 12)
tt_org <- tt

# an folgendem schritt w�hlen, welche artdistanz in die analyse geht
names(tt) # jac, sor - > q = 0, hor -> q = 1, mor_hor -> q = 2

#Distances observed species along Hill
dis_com_q1obs <- matrix(tt[,"hor_obs"], ncol = nrow(tab), byrow = T)
dis_com_q1obs <- (1 - dis_com_q1obs)# in dissimilarity convertieren
row.names(dis_com_q1obs)<-row.names(tab)



###############################################################################
#q=1
q1obs.NMDS <-
  metaMDS(dis_com_q1obs,
          #distance = "bray",
          k = 2,
          maxit = 999, 
          trymax = 100)

clrs <- c("deepskyblue4","darkslategray","deepskyblue1","chartreuse4","darkseagreen1","darksalmon","darkseagreen",
			"red","yellow","brown")


plot(q1obs.NMDS,display = "sites",type="p", main="Bird community, q=1") # 
ordihull(q1obs.NMDS, data_org$Category10, display = "sites", draw = c("polygon"), col=clrs )
ordispider(q1obs.NMDS, data_org$Category10,  
	 spiders = c("median"),  
         label = TRUE, col = clrs)

data_org$Birds_expert_Axis1<-q1obs.NMDS$points[,1]
data_org$Birds_expert_Axis2<-q1obs.NMDS$points[,2]

#Model of Vocalizing vertebrate total Richness# Table 1
modCom_SoundIndicesAx1<-lm(Birds_expert_Axis1 ~  TemporalEntropy + AcousticComplexity + EntropyOfVarianceSpectrum 
											+ SoundscapeSaturation + EventsPerSecond, data=data_org)
summary(modCom_SoundIndicesAx1)
plot(predict(modCom_SoundIndicesAx1), data_org$Birds_expert_Axis1, xlab="Predicted values",ylab="Observed values")
abline(lm(data_org$Birds_expert_Axis1 ~ predict(modCom_SoundIndicesAx1)))
text(-0.4,0.2,"R²=0.60")
ComAx1<- modCom_SoundIndicesAx1
											
summary(ComAx1)
Coefficients:
                          Estimate Std. Error t value Pr(>|t|)    
(Intercept)                2.41759    0.74430   3.248 0.001707 ** 
TemporalEntropy            2.99227    0.77653   3.853 0.000236 ***
AcousticComplexity        -8.18147    1.53297  -5.337 8.80e-07 ***
EntropyOfVarianceSpectrum  1.84451    0.46845   3.938 0.000176 ***
SoundscapeSaturation       0.72153    0.22387   3.223 0.001845 ** 
EventsPerSecond           -0.35386    0.06863  -5.156 1.82e-06 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.1717 on 79 degrees of freedom
Multiple R-squared:  0.6268,    Adjusted R-squared:  0.6031 
F-statistic: 26.53 on 5 and 79 DF,  p-value: 1.231e-15

#Calculate partial effects without residuals:
eff.ComAx1<- allEffects(ComAx1, residuals=FALSE)
plot(eff.ComAx1)





modCom_SoundIndicesAx2<-lm(Birds_expert_Axis2 ~  TemporalEntropy + AcousticComplexity + EntropyOfVarianceSpectrum 
											+ SoundscapeSaturation + EventsPerSecond, data=data_org)
summary(modCom_SoundIndicesAx2)
ComAx2<- modCom_SoundIndicesAx2
eff.ComAx2<- allEffects(ComAx2, residuals=FALSE)
plot(eff.ComAx2)

# Creating the distance by sound----------------------------------------------
names(data_org)
#TemporalEntropy 12# AcousticComplexity 11 # EntropyOfVarianceSpectrum 10 #	SoundscapeSaturation 9 # EventsPerSecond 13

TemporalEntropy <- as.data.frame(data_org[,12])
rownames(TemporalEntropy )<-rownames(data_org)
TemporalEntropy_s<-decostand(TemporalEntropy , "standardize")
#Standardize  to make estimates comparabel later
TemporalEntropy_s_dist <- dist(TemporalEntropy_s, method = "euclidean")
tt <- data.frame(as.matrix(TemporalEntropy_s_dist))
tt <- tt[(row.names(tt)),]
tt <- tt[,(colnames(tt))] 
dis_TemporalEntropy_s <- as.dist(tt, upper = T, diag = T)
all(row.names(as.matrix(dis_TemporalEntropy_s)) == row.names(as.matrix(dis_com_q1obs)))
####################################################################################
AcousticComplexity<- as.data.frame(data_org[,11])
rownames(AcousticComplexity)<-rownames(data_org)
AcousticComplexity_s<-decostand(AcousticComplexity, "standardize")
#Standardize  to make estimates comparabel later
AcousticComplexity_s_dist <- dist(AcousticComplexity_s, method = "euclidean")
tt <- data.frame(as.matrix(AcousticComplexity_s_dist))
tt <- tt[(row.names(tt)),]
tt <- tt[,(colnames(tt))] 
dis_AcousticComplexity_s <- as.dist(tt, upper = T, diag = T)
all(row.names(as.matrix(dis_AcousticComplexity_s)) == row.names(as.matrix(dis_com_q1obs)))
####################################################################################
EntropyOfVarianceSpectrum <- as.data.frame(data_org[,10])
rownames(EntropyOfVarianceSpectrum)<-rownames(data_org)
EntropyOfVarianceSpectrum_s<-decostand(EntropyOfVarianceSpectrum, "standardize")
#Standardize  to make estimates comparabel later
EntropyOfVarianceSpectrum_s_dist <- dist(EntropyOfVarianceSpectrum_s, method = "euclidean")
tt <- data.frame(as.matrix(EntropyOfVarianceSpectrum_s_dist))
tt <- tt[(row.names(tt)),]
tt <- tt[,(colnames(tt))] 
dis_EntropyOfVarianceSpectrum_s <- as.dist(tt, upper = T, diag = T)
all(row.names(as.matrix(dis_EntropyOfVarianceSpectrum_s)) == row.names(as.matrix(dis_com_q1obs)))

####################################################################################
SoundscapeSaturation<- as.data.frame(data_org[,9])
rownames(SoundscapeSaturation)<-rownames(data_org)
SoundscapeSaturation_s<-decostand(SoundscapeSaturation, "standardize")
#Standardize  to make estimates comparabel later
SoundscapeSaturation_s_dist <- dist(SoundscapeSaturation_s, method = "euclidean")
tt <- data.frame(as.matrix(SoundscapeSaturation_s_dist))
tt <- tt[(row.names(tt)),]
tt <- tt[,(colnames(tt))] 
dis_SoundscapeSaturation_s <- as.dist(tt, upper = T, diag = T)
all(row.names(as.matrix(dis_SoundscapeSaturation_s)) == row.names(as.matrix(dis_com_q1obs)))

####################################################################################
EventsPerSecond <- as.data.frame(data_org[,13])
rownames(EventsPerSecond)<-rownames(data_org)
EventsPerSecond_s<-decostand(EventsPerSecond, "standardize")
#Standardize  to make estimates comparabel later
EventsPerSecond_s_dist <- dist(EventsPerSecond_s, method = "euclidean")
tt <- data.frame(as.matrix(EventsPerSecond_s_dist))
tt <- tt[(row.names(tt)),]
tt <- tt[,(colnames(tt))] 
dis_EventsPerSecond_s <- as.dist(tt, upper = T, diag = T)
all(row.names(as.matrix(dis_EventsPerSecond_s)) == row.names(as.matrix(dis_com_q1obs)))


dis_com_q1obs <- as.dist(dis_com_q1obs, upper = T, diag = T)


MRM(dis_com_q1obs ~ dis_EventsPerSecond_s + dis_SoundscapeSaturation_s + dis_EntropyOfVarianceSpectrum_s 
		+ dis_AcousticComplexity_s+ dis_TemporalEntropy_s, method = "linear")









