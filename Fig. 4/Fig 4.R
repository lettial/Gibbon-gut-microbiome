####co-relationship #####
#整合数据并计算相关性
#fig. 4b,c
library(Hmisc)
core_otu_table <- read.table("core_otu.txt",sep = "\t",row.names = 1,header = T)
metadata <-read.table(file = "metadata",sep = "\t",row.names = 1,header = T)
core_otu_table <- data.frame(t(core_otu_table))
core_otu_table$sample_id <- row.names(core_otu_table)
metadata$sample_id <- row.names(metadata)
core_otu_diet <- merge(metadata,core_otu_table,by="sample_id")
row.names(core_otu_diet) <- core_otu_diet$sample_id
metadata_a1 <- subset(metadata,Subject == "A1")
core_otu_diet <- core_otu_diet[,-c(1:2)] 
core_a1 <- core_otu_diet[row.names(metadata_a1),]
data=as.matrix(core_a1)
corr=rcorr(data, type="spearman")
rcorr=corr$r
pcorr=corr$P
ecocor=as.matrix(rcorr[6:304, 1:5])
ecocop=as.matrix(pcorr[6:304, 1:5])
ecocor <- data.frame(ecocor)
ecocop <- data.frame(ecocop)

#adjust p-value
ecocop_adjust <- data.frame(apply(ecocop, 2, function(x) p.adjust(x, method = "fdr")))
write.table(ecocor, file = "a1_ecocor.txt", sep = "\t",quote = F)
write.table(ecocop, file = "a1_ecocop.txt", sep = "\t",quote = F)
write.table(ecocop_adjust, file = "a1_ecocop_adjust.txt", sep = "\t",quote = F)
ecocop_leaf <- subset(ecocop_adjust , Leaf < 0.05)
ecocop_fruit <- subset(ecocop_adjust , Fruit < 0.05)
ecocor_leaf <- subset(ecocor , Leaf > 0)
ecocor_fruit <- subset(ecocor , Fruit > 0)
r_leaf_select <- ecocor_leaf[c(intersect(rownames(ecocop_leaf), rownames(ecocor_leaf))),]
p_leaf_select <-data.frame(ecocop_leaf[c(intersect(rownames(ecocop_leaf), rownames(ecocor_leaf))),])
a1_leaf_select <- cbind(r_leaf_select,p_leaf_select)
r_fruit_select <- ecocor_fruit[c(intersect(rownames(ecocop_fruit), rownames(ecocor_fruit))),]
p_fruit_select <-ecocop_fruit[c(intersect(rownames(ecocop_fruit), rownames(ecocor_fruit))),]
a1_fruit_select <- cbind(r_fruit_select,p_fruit_select)
a1_leaf <- intersect(rownames(ecocor_leaf), rownames(ecocop_leaf))
a1_fruit <- intersect(rownames(ecocor_fruit), rownames(ecocop_fruit))
write.table(a1_leaf_select, file = "a1_leaf_select.txt", sep = "\t")
write.table(a1_fruit_select, file = "a1_fruit_select.txt", sep = "\t")

#Repeat the above code, replacing different individuals
all_leaf_intersect <- Reduce(intersect, list(c(a1_leaf),c(a2_leaf),c(a3_leaf),c(a4_leaf),c(b1_leaf),c(b2_leaf)))
all_fruit_intersect <- Reduce(intersect, list(c(a1_fruit),c(a2_fruit),c(a3_fruit),c(a4_fruit),c(b1_fruit),c(b2_fruit)))
A_all_leaf_intersect <- Reduce(intersect, list(c(a1_leaf),c(a2_leaf),c(a3_leaf),c(a4_leaf)))
A_all_fruit_intersect <- Reduce(intersect, list(c(a1_fruit),c(a2_fruit),c(a3_fruit),c(a4_fruit)))
B_all_leaf_intersect <- Reduce(intersect, list(c(b1_leaf),c(b2_leaf)))
B_all_fruit_intersect <- Reduce(intersect, list(c(b1_fruit),c(b2_fruit)))
all_leaf_intersect_data <- data.frame(t(core_otu_table)[c(all_leaf_intersect),])

View(nocopy_core_otu_data)
library(ggplot2)
library(ggpmisc)
library(ggthemes)
library(reshape2)
library(ggpubr)
all_leaf_intersect_data<- data.frame(t(all_leaf_intersect_data))
all_leaf_intersect_data$sample_id <- row.names(all_leaf_intersect_data)
all_leaf_intersect_data_3col <- melt(all_leaf_intersect_data, id.vars = "sample_id", variable.name = "Zotu", value.name = "abundance")
all_leaf_intersect_data_2 <- merge(metadata_2,all_leaf_intersect_data_3col,by = "sample_id")

p<-ggplot(all_leaf_intersect_data_2,aes(x = Leaf * 100, y = as.numeric(abundance)/82600*100,color= Subject)) +
  geom_smooth(method = 'lm', formula = y ~ x, se = T,aes(color= Subject)) + 
  stat_poly_eq(aes(label = paste(..adj.rr.label.., sep = '~~~~')), formula = y ~ x, parse = T) + 
  scale_linetype_manual() + 
  geom_point(size=1,aes(color= Subject)) +
  facet_wrap(Zotu~Subject,strip.position = c("top"),scales = "free_y",ncol = 6)+ 
  theme(strip.background = element_blank(), strip.placement = "outside",text = element_text(size = 5))+
  theme_classic()+
  scale_fill_manual(values = subject_color) +
  scale_color_manual(values = subject_color) +
  stat_fit_glance(method = "lm",  label.x = "centre",label.y = "top",aes(label = paste("italic(P)*\"-value = \"*", signif(..p.value.., digits = 4), sep = "")),parse = TRUE ) +
  xlab("Proportion of leaf (%)") + ylab("Relative abundance of Zotu (%)")
p

