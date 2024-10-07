

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


##############################################################################################

load(file= "data/distances_tax83.rda")
load(file= "data/distances_phy83.rda")
load(file= "data/distances_func83.rda")

#### load environment 

com_env <- read.csv2("data/EnvironAll86plots.csv", header=T)
names(com_env)

com_env$Plot == rownames(dis_birds_q0est)

com_env$Category4[com_env$Category10 == "CA_M"] <- "AG_M"
com_env$Category4[com_env$Category10 == "PA_M"] <- "AG_M"
com_env$Category4[com_env$Category10 == "OG_M"] <- "OG_M"
#######################################################################################

clr_data <- data.frame(clr=c("#ffd500","#ffd500","#92cb11","#92cb11","#359619","#359619","#18622a","#B1E3F0","#4974a5","#7a49a5"),
                       Category10 = c("Cacao","Pasture","C_Reg_I","P_Reg_I","C_Reg_II","P_Reg_II","Oldgrowth","PA_M","CA_M","OG_M"))

com_env$clrs <- clr_data[match(com_env$Category10, clr_data$Category10), "clr"]


clr_data <- data.frame(clr=c("#ffd500","#92cb11","#359619","#18622a","#4974a5","#7a49a5"),
                       Category4 = c("Agri","Reg1","Reg2","Old","AG_M","OG_M"))

com_env$clrs2 <- clr_data[match(com_env$Category4, clr_data$Category4), "clr"]
######################################################################################



png(filename="figures/nmds_ellipse4_sc83.png", width = 170, height = 200, units = "mm",
    res=1000)

layout(matrix(c(0,1,2,
                0,3,4,
                0,5,6,
                0,0,0),ncol=3,byrow=T),width=c(0.2,5,5),heights=c(5,5,5,0.3))
par(mar = c(2,2,2,0.5))


q0est.NMDS <-
  metaMDS(dis_birds_q0est, k = 2,maxit = 999,trymax = 100)

#names(table(com_env$Category10))
#clrs <- c("#92cb11","#359619","#4974a5","#ffd500","#7a49a5","#18622a","#92cb11","#359619","#B1E3F0","#ffd500")

names(table(com_env$Category4))
clrs <- c("#4974a5","#ffd500","#7a49a5","#18622a","#92cb11","#359619")




plot(q0est.NMDS,display = "sites",type="p",cex=0.8, main="Estimated taxonomic q=0")

ds <- data.frame(scores(q0est.NMDS, "sites"))
ds$Plot <- rownames(ds)
ds <- merge(com_env, ds, by="Plot")
points(ds$NMDS1, ds$NMDS2, pch=16, cex=1, col=ds$clrs2)


ordiellipse(q0est.NMDS, com_env$Category4, display = "sites", draw = c("polygon"), col=clrs, 
            border=clrs, label = TRUE)

mtext("NMDS 2", 2, line = 2, cex=0.7)
#------------------------------------------------------------------------------------------

q2est.NMDS <-
  metaMDS(dis_birds_q2est, k = 2,maxit = 999,trymax = 100)


plot(q2est.NMDS,display = "sites",type="p",cex=0.8, main="Estimated taxonomic q=2")

ds <- data.frame(scores(q2est.NMDS, "sites"))
ds$Plot <- rownames(ds)
ds <- merge(com_env, ds, by="Plot")
points(ds$NMDS1, ds$NMDS2, pch=16, cex=1, col=ds$clrs2)


ordiellipse(q2est.NMDS, com_env$Category4, display = "sites", draw = c("polygon"), col=clrs, 
            border=clrs, label = TRUE)


#------------------------------------------------------------------------------------------

q0est.NMDS <-
  metaMDS(dis_birds_q0est_pd, k = 2,maxit = 999,trymax = 100)


plot(q0est.NMDS,display = "sites",type="p",cex=0.8, main="Estimated phylogenetic q=0")

ds <- data.frame(scores(q0est.NMDS, "sites"))
ds$Plot <- rownames(ds)
ds <- merge(com_env, ds, by="Plot")
points(ds$NMDS1, ds$NMDS2, pch=16, cex=1, col=ds$clrs2)


ordiellipse(q0est.NMDS, com_env$Category4, display = "sites", draw = c("polygon"), col=clrs, 
            border=clrs, label = TRUE)

mtext("NMDS 2", 2, line = 2, cex=0.7)
#------------------------------------------------------------------------------------------

q2est.NMDS <-
  metaMDS(dis_birds_q2est_pd, k = 2,maxit = 999,trymax = 100)


plot(q2est.NMDS,display = "sites",type="p",cex=0.8, main="Estimated phylogenetic q=2")

ds <- data.frame(scores(q2est.NMDS, "sites"))
ds$Plot <- rownames(ds)
ds <- merge(com_env, ds, by="Plot")
points(ds$NMDS1, ds$NMDS2, pch=16, cex=1, col=ds$clrs2)


ordiellipse(q2est.NMDS, com_env$Category4, display = "sites", draw = c("polygon"), col=clrs, 
            border=clrs, label = TRUE)



#------------------------------------------------------------------------------------------

q0est.NMDS <-
  metaMDS(dis_birds_q0est_fd, k = 2,maxit = 999,trymax = 100)


plot(q0est.NMDS,display = "sites",type="p",cex=0.8, main="Estimated functional q=0")

ds <- data.frame(scores(q0est.NMDS, "sites"))
ds$Plot <- rownames(ds)
ds <- merge(com_env, ds, by="Plot")
points(ds$NMDS1, ds$NMDS2, pch=16, cex=1, col=ds$clrs2)


ordiellipse(q0est.NMDS, com_env$Category4, display = "sites", draw = c("polygon"), col=clrs, 
            border=clrs, label = TRUE)

mtext("NMDS 2", 2, line = 2, cex=0.7)
mtext("NMDS 1", 1, line = 2, cex=0.7)
#------------------------------------------------------------------------------------------

q2est.NMDS <-
  metaMDS(dis_birds_q2est_fd, k = 2,maxit = 999,trymax = 100)


plot(q2est.NMDS,display = "sites",type="p",cex=0.8, main="Estimated functional q=2")

ds <- data.frame(scores(q2est.NMDS, "sites"))
ds$Plot <- rownames(ds)
ds <- merge(com_env, ds, by="Plot")
points(ds$NMDS1, ds$NMDS2, pch=16, cex=1, col=ds$clrs2)


ordiellipse(q2est.NMDS, com_env$Category4, display = "sites", draw = c("polygon"), col=clrs, 
            border=clrs, label = TRUE)


mtext("NMDS 1", 1, line = 2, cex=0.7)

dev.off()
system("open figures/nmds_ellipse4_sc83.png")



