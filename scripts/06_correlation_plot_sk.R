rm(list=ls(all=TRUE))
Sys.setenv(LANG = "en")

# Working directory
setwd("E:/Manuscripts/1_3D_ecuador")


# load nmds values
nmds <- read.csv("data/plots_categories_indices_nmds.csv")


set.seed(123)


nmds_dataframe <- data.frame(
  TD_q0_Axis1 = nmds$TD_q0_Axis1,
  TD_q1_Axis1 = nmds$TD_q1_Axis1,
  TD_q2_Axis1 = nmds$TD_q2_Axis1,
  FD_q0_Axis1 = nmds$FD_q0_Axis1,
  FD_q1_Axis1 = nmds$FD_q1_Axis1,
  FD_q2_Axis1 = nmds$FD_q2_Axis1,
  PD_q0_Axis1 = nmds$PD_q0_Axis1,
  PD_q1_Axis1 = nmds$PD_q1_Axis1,
  PD_q2_Axis1 = nmds$PD_q2_Axis1
)

library(ggplot2)
library(reshape2)

# Calculate correlations
cor_matrix <- cor(nmds_dataframe)

# Melt correlation matrix
melted_cor <- melt(cor_matrix)


# Create heatmap
ggplot(melted_cor, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0.75, limit = c(0.5,1), space = "Lab",
                       name="Correlation") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA)
  ) +
  coord_fixed() +
  labs(title = "NMDS Correlation Heatmap",
       x = "NMDS1",
       y = "NMDS1")

# Save heatmap to file
ggsave("plots/nmds_correlation_heatmap.png", width = 8, height = 6, dpi = 300)
