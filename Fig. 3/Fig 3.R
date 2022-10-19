library(ggplot2)
library(ggpmisc)
library(ggpubr)
otu_network_result <- read.table("network_result.txt", sep = "\t", header = T)
p<-ggplot(otu_network_result,aes(x =Leaf*100, y = Index.value ,color = Group))   +
  geom_smooth(method = 'lm', formula = y ~ x, se = T) + 
  stat_poly_eq(aes(label = paste(..adj.rr.label.., sep = '~~~~')), formula = y ~ x, parse = T) + 
  scale_linetype_manual() + 
  geom_point(size=2) +
  facet_grid(index~Group,scales= "free" ) +
  stat_fit_glance(method = "lm",  label.x = "centre",label.y = "top",aes(label = paste("italic(P)*\"-value = \"*", signif(..p.value.., digits = 4), sep = "")),parse = TRUE ) +
  xlab("Proportion of leaf (%)") 
p