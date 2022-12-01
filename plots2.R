


#===============================================
load('~/scratch/principal_component_GWAS/sladek_gip1/gip1_dataframe.RData')
load('~/scratch/principal_component_GWAS/sladek_gip1/gip1_positions.RData')



colnames(SNPs) = c("rsid", "chr", "start", "stop")

gip1_positions <- merge(SNPs, gip1_dataframe,by="rsid")


gip1_positions = gip1_positions[which(gip1_positions$chr %in% c(1:22)),]

gip1_positions = gip1_positions[order(factor(gip1_positions$chr, levels = c(paste0(1:22)))),]


gip1_positions = gip1_positions[which(gip1_positions$start >= 0),]

save(gip1_positions, file = '~/scratch/principal_component_GWAS/sladek_gip1/gip1_positions_filtered.RData')


print(5)

#======
#plot
#=====

load('~/scratch/principal_component_GWAS/sladek_gip1/gip1_positions_filtered.RData')


# gip1_positions = na.omit(gip1_positions)
gip1_positions = gip1_positions[is.finite(gip1_positions$neg_log10_p_val),]



#addition
gip1_positions$chr = as.numeric(gip1_positions$chr)
gip1_positions = gip1_positions[order(gip1_positions$chr),]
#addition



library(ggrepel)
library(ggtext)
library(ggplot2)


# ggplot(gip1_dataframe, aes(1:nrow(gip1_dataframe), neg_log10_p_val)) + geom_point()
plot = ggplot(gip1_positions,
              aes(start, neg_log10_p_val, label = rsid)) +
  geom_point(size= 0.0005) +
  facet_grid(~chr,
             scales = "free_x",
             space = "free_x",
             switch = "x") +
  theme(axis.title = element_blank()) + theme(axis.title.x=element_blank(),
                                              axis.text.x=element_blank(),
                                              axis.ticks.x=element_blank()) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_text_repel(aes(label=ifelse(neg_log10_p_val>7.30102999566,as.character(rsid),'')),hjust=0,vjust=0, size=3) + 
  coord_cartesian(clip = "off") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
  
  
  
  
  
  
  
# elog10_textl>2,as.charater(rsid),'')),hjust=0,vjust=0)

print(6)

ggsave(file = "~/scratch/principal_component_GWAS/sladek_plots/manhattan_gip1.png", plot = plot, width = 12, height = 4)



top_snps = gip1_positions[order(gip1_positions$p_val),]

print("done 7")
