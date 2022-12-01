
library(ggplot2)

load('~/scratch/principal_component_GWAS/sladek_gip1/gip1_dataframe.RData')


#====================================
#genetic correlation heatmap
library(reshape2)
genetic_correlation = temp_correlation
rownames(genetic_correlation) = colnames(genetic_correlation)
genetic_correlation_df = melt(genetic_correlation)
colnames(genetic_correlation_df) = c("phenotype_1", "phenotype_2", "genetic_correlation")


heatplot_correlation = ggplot(genetic_correlation_df, aes(phenotype_1, ordered(phenotype_2, levels = rev(sort(unique(phenotype_2)))), fill= genetic_correlation)) + 
  geom_tile() +
  ggtitle("Heatmap of Genetic Correlation") +
  xlab("Phenotype 1") +
  ylab("Phenotype 2") +
  geom_text(aes(label = round(genetic_correlation, 2)), colour="red") + 
  theme(axis.text.x=element_text(angle=90,hjust=1))

ggsave(file = "~/scratch/principal_component_GWAS/sladek_plots/heatplot_correlation.png", plot = heatplot_correlation, width = 7, height = 5)





#second set of plot of gip2 was based on matrix of n's not a single n

#====================
#PC loading bar plots

library(reshape2)
rownames(eigenvectors) = colnames(genetic_covariance)
colnames(eigenvectors) = c(paste0("GIP-", 1:22))
eigenvectors_df = melt(eigenvectors)
eigenvectors_df$colour = ifelse(eigenvectors_df$value < 0,"negative","positive")

gip1_plot = ggplot(data = eigenvectors_df[which(eigenvectors_df$Var2 == "GIP-1"),], aes(y = value, x = Var1, fill = colour)) +
  geom_col() +
  coord_flip()+ theme(legend.position = "none") + scale_fill_manual(values=c(positive="firebrick1",negative="steelblue")) + 
  scale_y_continuous(breaks = round(seq(min(eigenvectors_df$value), max(eigenvectors_df$value), by = 0.5),1))+ 
  xlab("Phenotype") +
  ylab("PC Loading")


gip2_plot = ggplot(data = eigenvectors_df[which(eigenvectors_df$Var2 == "GIP-2"),], aes(y = value, x = Var1, fill = colour)) +
  geom_col() +
  coord_flip()+ theme(legend.position = "none") + scale_fill_manual(values=c(positive="firebrick1",negative="steelblue")) + 
  scale_y_continuous(breaks = round(seq(min(eigenvectors_df$value), max(eigenvectors_df$value), by = 0.5),1))+ 
  xlab("Phenotype") +
  ylab("PC Loading")

gip3_plot = ggplot(data = eigenvectors_df[which(eigenvectors_df$Var2 == "GIP-3"),], aes(y = value, x = Var1, fill = colour)) +
  geom_col() +
  coord_flip()+ theme(legend.position = "none") + scale_fill_manual(values=c(positive="firebrick1",negative="steelblue")) + 
  scale_y_continuous(breaks = round(seq(min(eigenvectors_df$value), max(eigenvectors_df$value), by = 0.5),1))+ 
  xlab("Phenotype") +
  ylab("PC Loading")

gip4_plot = ggplot(data = eigenvectors_df[which(eigenvectors_df$Var2 == "GIP-4"),], aes(y = value, x = Var1, fill = colour)) +
  geom_col() +
  coord_flip()+ theme(legend.position = "none") + scale_fill_manual(values=c(positive="firebrick1",negative="steelblue")) + 
  scale_y_continuous(breaks = round(seq(min(eigenvectors_df$value), max(eigenvectors_df$value), by = 0.5),1))+ 
  xlab("Phenotype") +
  ylab("PC Loading")

gip5_plot = ggplot(data = eigenvectors_df[which(eigenvectors_df$Var2 == "GIP-5"),], aes(y = value, x = Var1, fill = colour)) +
  geom_col() +
  coord_flip()+ theme(legend.position = "none") + scale_fill_manual(values=c(positive="firebrick1",negative="steelblue")) + 
  scale_y_continuous(breaks = round(seq(min(eigenvectors_df$value), max(eigenvectors_df$value), by = 0.5),1))+ 
  xlab("Phenotype") +
  ylab("PC Loading")

gip6_plot = ggplot(data = eigenvectors_df[which(eigenvectors_df$Var2 == "GIP-6"),], aes(y = value, x = Var1, fill = colour)) +
  geom_col() +
  coord_flip()+ theme(legend.position = "none") + scale_fill_manual(values=c(positive="firebrick1",negative="steelblue")) + 
  scale_y_continuous(breaks = round(seq(min(eigenvectors_df$value), max(eigenvectors_df$value), by = 0.5),1))+ 
  xlab("Phenotype") +
  ylab("PC Loading")

gip7_plot = ggplot(data = eigenvectors_df[which(eigenvectors_df$Var2 == "GIP-7"),], aes(y = value, x = Var1, fill = colour)) +
  geom_col() +
  coord_flip()+ theme(legend.position = "none") + scale_fill_manual(values=c(positive="firebrick1",negative="steelblue")) + 
  scale_y_continuous(breaks = round(seq(min(eigenvectors_df$value), max(eigenvectors_df$value), by = 0.5),1))+ 
  xlab("Phenotype") +
  ylab("PC Loading")

gip8_plot = ggplot(data = eigenvectors_df[which(eigenvectors_df$Var2 == "GIP-8"),], aes(y = value, x = Var1, fill = colour)) +
  geom_col() +
  coord_flip()+ theme(legend.position = "none") + scale_fill_manual(values=c(positive="firebrick1",negative="steelblue")) + 
  scale_y_continuous(breaks = round(seq(min(eigenvectors_df$value), max(eigenvectors_df$value), by = 0.5),1))+ 
  xlab("Phenotype") +
  ylab("PC Loading")

library(ggpubr)
bar_plot = ggarrange(gip1_plot, gip2_plot, gip3_plot, gip4_plot, gip5_plot, gip6_plot, gip7_plot, gip8_plot,
                     labels = paste0("GIP-", 1:8),
                     ncol = 8, nrow = 1)


ggsave(file = "~/scratch/principal_component_GWAS/sladek_plots/bar_plot.png", plot = bar_plot, width = 20, height = 3)




#==========================================================

