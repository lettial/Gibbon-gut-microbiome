library(ggplot2)
library(ggdendro)
library(remotes)
library(tidyverse)
library(ggtree)
library(ggstance)
library(ggnewscale)
library(ggpubr)
library(ggpmisc)
#### fig 2a,b ####
diet_nk <- read.table('nk_diet',header = T,sep = "\t",row.names = 1)
diet_bc <- read.table('bc_diet',header = T,sep = "\t",row.names = 1)
diet_nk_dist <- vegdist(diet_nk, method = "bray")
hc_nk <- hclust(diet_nk_dist, method = "average")
clus <- cutree(hc_nk, 2)
d = data.frame(label = names(clus), 
               member = factor(clus))
dd = merge(d,diet_nk,by = "row.names",all = F)
row.names(dd) = dd$Row.names 
dd$Row.names = NULL
hc_nk_backup <- hc_nk
p  = ggtree(hc_nk) %<+% dd + 
  geom_tippoint(size=.01, shape= 21, aes( x=x+.01,color= member,fill= member))+
  geom_tiplab(aes(x=x*1.2,color = member), size=.7) 
p  = ggtree(hc_nk) 
p
diet_nk$sample_id <- rownames(diet_nk)
diet_nk_3col <- melt(diet_nk, id.vars = "sample_id", variable.name = "Diet", value.name = "RA")
p <- p + new_scale_fill()
p3 <- facet_plot(p, panel = 'Stacked Barplot', 
                 data = diet_a_3col, geom = geom_barh,
                 mapping = aes(x = abundance, fill = Sample),
                 color = NA,stat='identity' ) 
p3
#Repeat the above code to get BC result

### fig. 2c ####

otu_table <- read.table("otutable.txt",sep = "\t",row.names = 1,header = T)
otu_table <- data.frame(t(otu_table))
otu_table_bc <- otu_table[row.names(diet_bc),]
otu_table_bc_dist <- vegdist(otu_table_bc,method = "bray")
#Calculate NK result
### clr-transfornmation of the counts ###
#for more information on compositional microbiome analysis, see Gloor et al. 2017 and Gloor's github tutorial (https://github.com/ggloor/CoDa_microbiome_tutorial/wiki/Part-1:-Exploratory-Compositional-PCA-biplot)
pseudocount <- 0.65
clr <- t(apply(t(otu_table_nk)+pseudocount, 2, compositions::clr))
dist(clr)
### Perform the PCA ###
pca <- prcomp(otu_table_bc_dist)
# Calculate the variance explained by PC1 and PC2
d.mvar <- sum(pca$sdev^2) # total variance
PC1_bc <- paste("PC1 ","(",round(sum(pca$sdev[1]^2)/d.mvar, 3)*100,"%",")",sep="")

#Add metadata
df<-data.frame(pca$x)
df$sample_id<-rownames(df)
df_bc<-merge(diet_bc,df,by="sample_id",all.x=TRUE,all.y=FALSE)
#Calculate BC result
otu_table_bc <- otu_table[row.names(diet_bc),]
pseudocount <- 0.65
clr <- t(apply(t(otu_table_nk)+pseudocount, 2, compositions::clr))
dist(clr)
pca <- prcomp(clr)
d.mvar <- sum(pca$sdev^2) # total variance
PC1_bc <- paste("PC1 ","(",round(sum(pca$sdev[1]^2)/d.mvar, 3)*100,"%",")",sep="")
df<-data.frame(pca$x)
df$sample_id<-rownames(df)
df_nk<-merge(diet_nk,df,by="sample_id",all.x=TRUE,all.y=FALSE)
#merge NK and BC result
dat1<-data.frame(rbind(df_nk[,1:7],df_b[,1:7]),
                 Group=rep(c("NK","BC"),c(159,115)))

pcoa_bray_gg <- ggplot(df_bc,aes(x=leaf*100,y=PC1)) 
pcoa_bray_gg + geom_point()   +
  theme(panel.grid.minor = element_blank()) + 
  theme(axis.title= element_text(size=12, color="black", face="bold", vjust=0.5, 
                                 hjust=0.5))  + labs(x="Proportion of leaf (%)",y="PC1") +
  geom_smooth(method = 'lm', formula = y ~ x, se = F) +
  stat_poly_eq(aes(label = paste(..adj.rr.label.., sep = '~~~~')), formula = y ~ x, parse = T) +
  stat_fit_glance(method = "lm",  label.x = "centre",label.y = "top",aes(label = paste("italic(P)*\"-value = \"*", signif(..p.value.., digits = 4), sep = "")),parse = TRUE )
  


####fig 2e f ######
### Loading scores of taxa on PC1 ###
taxa <- read.table('taxonomy.txt',sep = "\t",row.names = 1,header = T)
pca$rotation[1:5,1:5]
t1<-data.frame(pca$rotation)
t1$taxa_id<-rownames(t1)
t1[1:5,1:5]
barplot(sort(t1$PC1,decreasing=TRUE),ylab="Loading PC1")
taxa1<-merge(taxonomy,t1[c("taxa_id","PC1")],by="taxa_id")
taxa2<-taxa1[with(taxa1,order(-PC1)), ] 
head(taxa2,n=20)

#### fig 2ghij #####
#alpha diversity 
library(microbiome)
alpha <- alpha(data.frame(t(otu_table)), index = "all")
nocopy_alpha$sample_id <- row.names(nocopy_alpha)
alpha_data <- merge(nocopy_alpha,metadata_2,by="sample_id")
p_1<-ggplot(alpha_data,aes(x = Leaf*100, y = diversity_shannon ,group = Group, color = Group))  +
  geom_smooth(method = 'lm', formula = y ~ x, se = T, size = 2) + 
  stat_poly_eq(aes(label = paste(..adj.rr.label.., sep = '~~~~')), formula = y ~ x, parse = T) + 
  scale_linetype_manual() + 
  geom_point(size=2) + 
   xlab("Proportion of leaf (%)") + ylab("Shannon diversity") +
  stat_fit_glance(method = "lm",  label.x = "centre",label.y = "top",aes(label = paste("italic(P)*\"-value = \"*", signif(..p.value.., digits = 4), sep = "")),parse = TRUE )
p_1
p_2<-ggplot(alpha_data,aes(x = Fruit*100, y = diversity_shannon ,group = Group, color = Group))  +
  geom_smooth(method = 'lm', formula = y ~ x, se = T, size = 2) + 
  stat_poly_eq(aes(label = paste(..adj.rr.label.., sep = '~~~~')), formula = y ~ x, parse = T) + 
  scale_linetype_manual() + 
  geom_point(size=2) + 
  xlab("Proportion of fruit (%)") + ylab("Shannon diversity") +
  stat_fit_glance(method = "lm",  label.x = "centre",label.y = "top",aes(label = paste("italic(P)*\"-value = \"*", signif(..p.value.., digits = 4), sep = "")),parse = TRUE )
p_2
p_3<-ggplot(alpha_data,aes(x = Leaf*100, y = evenness_pielou ,group = Group, color = Group))  +
  geom_smooth(method = 'lm', formula = y ~ x, se = T, size = 2) + 
  stat_poly_eq(aes(label = paste(..adj.rr.label.., sep = '~~~~')), formula = y ~ x, parse = T) + 
  scale_linetype_manual() + 
  geom_point(size=2) + 
  xlab("Proportion of leaf (%)") + ylab("Communitity Evenness") +
  stat_fit_glance(method = "lm",  label.x = "centre",label.y = "top",aes(label = paste("italic(P)*\"-value = \"*", signif(..p.value.., digits = 4), sep = "")),parse = TRUE )
p_3
p_4<-ggplot(alpha_data,aes(x = Fruit*100, y = evenness_pielou ,group = Group, color = Group))  +
  geom_smooth(method = 'lm', formula = y ~ x, se = T, size = 2) + 
  stat_poly_eq(aes(label = paste(..adj.rr.label.., sep = '~~~~')), formula = y ~ x, parse = T) + 
  scale_linetype_manual() + 
  geom_point(size=2) + 
  xlab("Proportion of fruit (%)") + ylab("Communitity Evenness") +
  stat_fit_glance(method = "lm",  label.x = "centre",label.y = "top",aes(label = paste("italic(P)*\"-value = \"*", signif(..p.value.., digits = 4), sep = "")),parse = TRUE )
p_4
ggarrange(p_1,p_2,p_3,p_4,ncol = 2,nrow = 2)
