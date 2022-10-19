library(reshape2)
library(ggplot2)
#### fig1 taxonomy barplox ####
phylum <- read.table('phylum.txt',sep = "\t",row.names = 1,header = T)
family <- read.table('family.txt',sep = "\t",row.names = 1,header = T)
metadata <- read.table(metadata.txt,sep = "\t",row.names = 1,header = T)
phylum$sample_id <- rownames(phylum)
phylum_3col <- melt(phylum, id.vars = "sample_id", variable.name = "Phylum", value.name = "Relative_Abundance")
family$sample_id <- rownames(family)
family <- melt(family, id.vars = "sample_id", variable.name = "Family", value.name = "Relative_Abundance")
metadata$sample_id <- rownames(metadata)
phylum_data <- merge(phylum_3col, metadata, by = "sample_id")
gg <- ggplot(phylum_data) + 
  geom_bar(aes(x=sample_id, y = Relative_Abundance,fill=Phylum,
                   color= Phylum),position = "fill",stat="identity") +
  theme_bw()
gg

gg <- ggplot(family_data) + 
  geom_bar(aes(x=sample_id, y = RA,fill=Phylum,
               color= Family),position = "fill",stat="identity") +
  theme_bw()
gg

