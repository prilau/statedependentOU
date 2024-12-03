library(ape)
library(tidyverse)
library(viridis)
library(ggtree)

tree <- read.tree("data/3_empirical/mammal_perMY_n2955.tre")
sp_order <- list()
for (i in 1:length(tree$tip.label)){
  tmp_order <- df_hpa$order[which(df_hpa$binomial_name == tree$tip.label[i])]
  sp_order[[tmp_order]] <- append(sp_order[[tmp_order]], tree$tip.label[i])
}

tree <- tidytree::groupOTU(tree, sp_order)

colors <- rocket(length(unique(df_hpa$order)))
names(colors) <- unique(df_hpa$order)

p <- ggtree(tree, layout = "circular", aes(color=group), linewidth = 0.1) +
  ggtree::geom_tippoint(ggplot2::aes(colour = group), size = 0.1, shape = 1, alpha = 0.00) +
  #geom_tiplab(aes(color=group), size = 0.4, offset = 0.5) +
  scale_color_manual(values=colors) +
  guides(color = guide_legend(title="Order",
    override.aes=list(size = 4, shape = 20, alpha = 1.0))) +
  theme(legend.title = element_text(hjust = 0.5, size = 14)) 
p

ggsave(paste0("figures/3_empirical/phylo_by_order.pdf"), p, width = 8, height = 6)








