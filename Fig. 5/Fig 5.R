#significant difference ko
library(tidyverse)
data <- read.table('kegg.txt',sep = "\t",row.names = 1,header = T)
#NK result
kegg_nk <- data[row.names(metadata_nk),]
group <- metadata_nk
data <- data.frame(t(kegg_nk))
#Remove KO with too low abundance
count <- data %>% filter(apply(data,1,mean) > 300)
data <- t(count)
data1 <- data.frame(data,group$leaf_proportion)
colnames(data1) <- c(colnames(data),"Group")
data1$Group <- as.factor(data1$Group)
diff <- data1 %>%
  select_if(is.numeric) %>%
  map_df(~ broom::tidy(wilcox.test(. ~ Group,data = data1, alternative = "two.sided")), .id = 'var')
#p.adjust
diff$p.adj <- p.adjust(diff$p.value,"fdr")
diff_b <- diff %>% filter(p.adj < 0.05)

#BC result
kegg_bc <- metagenomic_ko[row.names(metadata_bc),]
group <- metadata_bc
data <- data.frame(t(kegg_bc))
#Remove KO with too low abundance
count <- data %>% filter(apply(data,1,mean) > 300)
data <- t(count)
data1 <- data.frame(data,group$leaf_proportion)
colnames(data1) <- c(colnames(data),"Group")
data1$Group <- as.factor(data1$Group)
diff <- data1 %>%
  select_if(is.numeric) %>%
  map_df(~ broom::tidy(wilcox.test(. ~ Group,data = data1, alternative = "two.sided")), .id = 'var')
#p.adjust
diff$p.adj <- p.adjust(diff$p.value,"fdr")
diff_bc <- diff %>% filter(p.adj < 0.05)
