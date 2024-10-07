rm(list=ls())

library(vegan)



## ================================================================================================ ##
## Use abundance data / community composition "BirdCommunityNamePhylogeny_Indices.csv" from Joerg.
## There are 85 assemblages (rows) and 334 birds (columns), identified by experts.
## The table contains also 8 different sound indices and plot categories (Pasture/Cacao -> Old Growth)

## Compute and plot nonmetric multidimensional scaling (nmds)


# Load community composition data containing plot categories
data_com <- read.csv2("data/BirdCommunityNamePhylogeny_Indices.csv", row.names=1,header=T)


# Load taxonomic dissimilarity matrices
load(file = "taxonomic_diversity/distances_com_exp_tax.rda")

# Load phylogenetic dissimilarity matrices
load(file = "phylogenetic_diversity/distances_com_exp_phy.rda")

# Load functional dissimilarity matrices
load(file = "functional_diversity/distances_com_exp_func.rda")


# Define color palettes
palette_Forest <- c("#FFFF00", "greenyellow", "forestgreen", "#006400")
palette_Agriculture <- c("orange", "sienna")

names(table(data_com$Category10))
# > "A_Caca"  "A_Old"   "A_Past"  "F_Caca"  "F_CReg1" "F_CReg2" "F_Old"   "F_Past"  "F_PReg1" "F_PReg2"

# Manual color mapping
color_mapping <- c(
  "A_Caca" = "orange",
  "A_Old" = "sienna",
  "A_Past" = "orange",
  "F_Caca" = "yellow",
  "F_CReg1" = "greenyellow",
  "F_CReg2" = "chartreuse3",
  "F_Old" = "darkgreen",
  "F_Past" = "yellow",
  "F_PReg1" = "greenyellow",
  "F_PReg2" = "chartreuse3"
)

# Assign colors based on the mapping and add them to the table
data_com$Color <- color_mapping[data_com$Category10]


######################## Perform Nonmetric Multidimensional Scaling (NMDS) ###########################

png(filename="plots/nmds.png", width = 170, height = 200, units = "mm",
    res=1000)

layout(matrix(c(1, 2, 3, 
                4, 5, 6, 
                7, 8, 9), ncol=3, byrow=TRUE), widths=c(1, 1, 1), heights=c(1, 1, 1))
par(mar = c(2, 2, 2, 1), oma = c(4, 4, 4, 1))

# Helper function to plot NMDS, add hulls and spiders
plot_nmds <- function(nmds) {
  ordiplot(nmds, display = "sites", type="n")
  points(nmds, display = "sites", col=data_com$Color, pch=19, cex=0.8)
  ordihull(nmds, data_com$Category10, display = "sites", draw = "polygon", col=color_mapping)
  ordispider(nmds, data_com$Category10, spiders = "median", col=color_mapping, label = TRUE)
}

# ---------------- TD q0 ----------------
tdq0est.NMDS <- metaMDS(dis_com_exp_q0est, k = 2, maxit = 999, trymax = 100)
print(tdq0est.NMDS$stress)
data_com$TD_q0_est_Axis1 <- tdq0est.NMDS$points[,1]
data_com$TD_q0_est_Axis2 <- tdq0est.NMDS$points[,2]
plot_nmds(tdq0est.NMDS)

# ---------------- TD q1 ----------------
tdq1est.NMDS <- metaMDS(dis_com_exp_q1est, k = 2, maxit = 999, trymax = 100)
print(tdq1est.NMDS$stress)
data_com$TD_q1_est_Axis1 <- tdq1est.NMDS$points[,1]
data_com$TD_q1_est_Axis2 <- tdq1est.NMDS$points[,2]
plot_nmds(tdq1est.NMDS)

# ---------------- TD q2 ----------------
tdq2est.NMDS <- metaMDS(dis_com_exp_q2est, k = 2, maxit = 999, trymax = 100)
print(tdq2est.NMDS$stress)
data_com$TD_q2_est_Axis1 <- tdq2est.NMDS$points[,1]
data_com$TD_q2_est_Axis2 <- tdq2est.NMDS$points[,2]
plot_nmds(tdq2est.NMDS)

# ---------------- FD q0 ----------------
fdq0est.NMDS <- metaMDS(dis_com_exp_q0est_fd, k = 2, maxit = 999, trymax = 100)
print(fdq0est.NMDS$stress)
data_com$FD_q0_est_Axis1 <- fdq0est.NMDS$points[,1]
data_com$FD_q0_est_Axis2 <- fdq0est.NMDS$points[,2]
plot_nmds(fdq0est.NMDS)

# ---------------- FD q1 ----------------
fdq1est.NMDS <- metaMDS(dis_com_exp_q1est_fd, k = 2, maxit = 999, trymax = 100)
print(fdq1est.NMDS$stress)
data_com$FD_q1_est_Axis1 <- fdq1est.NMDS$points[,1]
data_com$FD_q1_est_Axis2 <- fdq1est.NMDS$points[,2]
plot_nmds(fdq1est.NMDS)

# ---------------- FD q2 ----------------
fdq2est.NMDS <- metaMDS(dis_com_exp_q2est_fd, k = 2, maxit = 999, trymax = 100)
print(fdq2est.NMDS$stress)
data_com$FD_q2_est_Axis1 <- fdq2est.NMDS$points[,1]
data_com$FD_q2_est_Axis2 <- fdq2est.NMDS$points[,2]
plot_nmds(fdq2est.NMDS)

# ---------------- PD q0 ----------------
pdq0est.NMDS <- metaMDS(dis_com_exp_q0est_pd, k = 2, maxit = 999, trymax = 100)
print(pdq0est.NMDS$stress)
data_com$PD_q0_est_Axis1 <- pdq0est.NMDS$points[,1]
data_com$PD_q0_est_Axis2 <- pdq0est.NMDS$points[,2]
plot_nmds(pdq0est.NMDS)

# ---------------- PD q1 ----------------
pdq1est.NMDS <- metaMDS(dis_com_exp_q1est_pd, k = 2, maxit = 999, trymax = 100)
print(pdq1est.NMDS$stress)
data_com$PD_q1_est_Axis1 <- pdq1est.NMDS$points[,1]
data_com$PD_q1_est_Axis2 <- pdq1est.NMDS$points[,2]
plot_nmds(pdq1est.NMDS)

# ---------------- PD q2 ----------------
pdq2est.NMDS <- metaMDS(dis_com_exp_q2est_pd, k = 2, maxit = 999, trymax = 100)
print(pdq2est.NMDS$stress)
data_com$PD_q2_est_Axis1 <- pdq2est.NMDS$points[,1]
data_com$PD_q2_est_Axis2 <- pdq2est.NMDS$points[,2]
plot_nmds(pdq2est.NMDS)


# Add labels for rows and columns
mtext("q0", side=3, line=0.5, outer=TRUE, at=0.167, cex=1.2)
mtext("q1", side=3, line=0.5, outer=TRUE, at=0.5, cex=1.2)
mtext("q2", side=3, line=0.5, outer=TRUE, at=0.833, cex=1.2)

mtext("Taxonomic Diversity", side=2, line=2, outer=TRUE, at=0.83, cex=1.2)
mtext("Functional Diversity", side=2, line=2, outer=TRUE, at=0.5, cex=1.2)
mtext("Phylogenetic Diversity", side=2, line=2, outer=TRUE, at=0.17, cex=1.2)

dev.off()

write.csv(data_com, file = "data/BirdCommunityNamePhylogeny_Indices_nmds.csv")
