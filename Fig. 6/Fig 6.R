#fig 7a #
library(ggheatmap)
cazy_select <- read.table("cazy_select_result.txt",sep = "\t",row.names = 1,header = T)
cazy_select <- data.frame(t(cazy_select))
nk_cazy_select <- cazy_select[row.names(metadata_nk),]
bc_cazy_select <- cazy_select[row.names(metadata_bc),]
group <- read.table("group.txt",sep = "\t",row.names = 1, header = T)
group$Sample <- rownames(group)
group_nk <- group[row.names(metadata_nk),]
group_bc <- group[row.names(metadata_bc),]
write.table(group_a, "group_nk.txt", sep = "\t", quote = F)
group_nk <- read.table("group_nk.txt", header = T, sep = "\t", row.names = 1)
group_nk_2 <- group_a
group_nk_2$Sample <- rownames(group_nk_2)
group_nk_2 <- group_nk_2[order(group_nk_2$group),]
names(group_nk) <- "Sampletype"
cazy_group <- read.table("cazy_group.txt",sep = "\t",row.names = 1, header = T)
names(cazy_group) <- "GHtype"
cazy_group_2 <- cazy_group
cazy_group_2$GH <- rownames(cazy_group_2)
Samplecol <- c("#EE0000FF","#008B45FF" )
names(Samplecol) <- c("F","L")
GHcol <- c("#845EC2","#D65DB1","#FF6F91","#FF9671","#FFC75F","#F9F871")
GHcol_2 <- c("79b99e","#326e8c","#845EC2","#f19c5c","#cd3b3e","#f7efcf")
names(GHcol_2) <- c("Cellulase","Hemicellulase","Pectinase","Debranching enzymes", "Amylases", "Oligosaccharide")
col <- list(Sampletype=Samplecol, GHtype=GHcol_2)
p<- ggheatmap(log(as.matrix(t(a_cazy_select)+1)),cluster_rows = T,cluster_cols = F,scale = "row",
              text_show_rows = rownames(t(nk_cazy_select)),
              annotation_rows = cazy_group,
              color=colorRampPalette(c( "#04027B","white","#BF201D"))(100),
              annotation_cols = group_nk,
              annotation_color = col,
              levels_cols = group_nk_2$Sample,
              levels_rows = cazy_group_2$GH
)

p


p%>%
  ggheatmap_theme(1,
                  theme =list(
                    theme(axis.text.x = element_text(angle = 90,face = "bold",size = 10),
                          axis.text.y = element_text(face = "bold")),
                    theme(legend.title = element_text(face = "bold")),
                    theme(legend.title = element_text(face = "bold")),
                    theme(legend.title = element_text(face = "bold")),
                    theme(legend.title = element_text(face = "bold"))
                  ))
