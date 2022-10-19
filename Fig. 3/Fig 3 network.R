# Load essential packages
library(SpiecEasi)
library(igraph)
library(Matrix)
source("info.centrality.R") # download from https://github.com/Mengting-Maggie-Yuan/warming-network-complexity-stability

# Customize function 
Stand_Error <- function(x){
  sd(x)/length(x)
}

# Other function (https://github.com/Mengting-Maggie-Yuan/warming-network-complexity-stability), including rand.remov.once and rmsimu.

# read file
otu_table_t <- read.table('otutable_above_50.txt', header=T, sep="\t", row.names = 1)
otu_table <- data.frame(t(otu_table_t))
rm.p.list=seq(0.05,0.2,by=0.05)
# Defined subset
A1L_F <- c("A030","A035","A040","A043","A053","A060","A065","A068","A084","A088","A092","D011","D014","D022","D023","D042","D049","D059","D060","D062","D069","D075","E008","E014","E018","E027","E028") # changed according to your case
rarefy_times <- 999 # permutation times
min_sample_num <- 19 # changed according to your case
A1L_F_otu_table <- otu_table[A1L_F,]
A1L_F_otu_table_t <- data.frame(t(A1L_F_otu_table))
A1L_F_otu_table_t_t <- subset(A1L_F_otu_table_t, (Reduce('+', as.data.frame(A1L_F_otu_table_t > 0)))>13)
A1L_F_otu_table_t_t_t <- data.frame(t(A1L_F_otu_table_t_t))

# predefine list name
speic_list_A1L <- list()
AD_A1L <- list()
AP_A1L <- list()
CO_A1L <- list()
NE_A1L <- list()
TL_A1L <- list()
ED_A1L <- list()
TA_A1L <- list()
NC_A1L <- list()
CB_A1L <- list()
CD_A1L <- list()
MD_A1L <- list()
VU_A1L <- list()
WS_A1L <- list()
UWS_A1L <- list ()

# permutations, for more information, see https://github.com/Mengting-Maggie-Yuan/warming-network-complexity-stability
for (i in 1:rarefy_times){
  otu_table_rarefy <- A1L_F_otu_table_t_t_t[sample(1:nrow(A1L_F_otu_table_t_t_t), min_sample_num),]
  otu_table_rarefy_f <- otu_table_rarefy[,colSums(otu_table_rarefy) > 4]
  comm<-otu_table_rarefy_f
  sp.ra<-colMeans(comm)/82600
  cormatrix=matrix(0,ncol(comm),ncol(comm))

  for (a in 1:ncol(comm)){
    for (j in a:ncol(comm)){
      speciesi<-sapply(1:nrow(comm),function(k){
        ifelse(comm[k,a]>0,comm[k,a],ifelse(comm[k,j]>0,0.01,NA))
      })
      speciesj<-sapply(1:nrow(comm),function(k){
        ifelse(comm[k,j]>0,comm[k,j],ifelse(comm[k,a]>0,0.01,NA))
      })
      corij<-cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
      cormatrix[a,j]<-cormatrix[j,a]<-corij
    
    }}

  row.names(cormatrix)<-colnames(cormatrix)<-colnames(comm)

  cormatrix2<-cormatrix*(abs(cormatrix)>=0.80)  #only keep links above the cutoff point
  cormatrix2[is.na(cormatrix2)]<-0
  diag(cormatrix2)<-0
  network.raw<-cormatrix2[colSums(abs(cormatrix2))>0,colSums(abs(cormatrix2))>0]
  sp.ra2<-sp.ra[colSums(abs(cormatrix2))>0]
  Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
  Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=F,nperm=100)
  speic <- spiec.easi(data = as.matrix(otu_table_rarefy_f), method = "mb", lambda.min.ratio = 0.001, nlambda = 20, 
                      pulsar.params = list(rep.num=50, ncores=20))
  # Store the result
  speic_list_A1L[[i]] <- otu_table_rarefy_f 
  speic_ig <- adj2igraph(getRefit(speic), vertex.attr=list(label=colnames(otu_table_rarefy_f)))
  iso_node_id = which(degree(speic_ig)==0)
  g2 = delete.vertices(speic_ig, iso_node_id)
  node.vul<-info.centrality.vertex(g2)
  VU_A1L[i] <- max(node.vul)
  WS_A1L[i] <- Weighted.simu[10,1]
  UWS_A1L[i] <- Unweighted.simu[10,1]
  AD_A1L[i] <- mean(igraph::degree(speic_ig))
  AP_A1L[i] <- average.path.length(speic_ig)
  CO_A1L[i] <- edge_connectivity(speic_ig)
  NE_A1L[i] <- length(E(speic_ig))
  TL_A1L[i] <- length(V(speic_ig))
  ED_A1L[i] <- edge_density(speic_ig, loops=F)
  TA_A1L[i] <- transitivity(speic_ig)
  NC_A1L[i] <- no.clusters(speic_ig)
  CB_A1L[i] <- centralization.betweenness(speic_ig)$centralization
  CD_A1L[i] <- centralization.degree(speic_ig)$centralization
  wtc <- cluster_walktrap(speic_ig)
  MD_A1L[i] <- modularity(wtc)
  print(paste("A1L",i,sep="-"))
  }
AD_A1L_mean <- mean(unlist(AD_A1L))
AD_A1L_sd <- Stand_Error(unlist(AD_A1L))
AP_A1L_mean <- mean(unlist(AP_A1L))
AP_A1L_sd <- Stand_Error(unlist(AP_A1L))
CO_A1L_mean <- mean(unlist(CO_A1L))
CO_A1L_sd <- Stand_Error(unlist(CO_A1L))
NE_A1L_mean <- mean(unlist(NE_A1L))
NE_A1L_sd <- Stand_Error(unlist(NE_A1L))
TL_A1L_mean <- mean(unlist(TL_A1L))
TL_A1L_sd <- Stand_Error(unlist(TL_A1L))
ED_A1L_mean <- mean(unlist(ED_A1L))
ED_A1L_sd <- Stand_Error(unlist(ED_A1L))
TA_A1L_mean <- mean(unlist(TA_A1L))
TA_A1L_sd <- Stand_Error(unlist(TA_A1L))
NC_A1L_mean <- mean(unlist(NC_A1L))
NC_A1L_sd <- Stand_Error(unlist(NC_A1L))
CB_A1L_mean <- mean(unlist(CB_A1L))
CB_A1L_sd <- Stand_Error(unlist(CB_A1L))
CD_A1L_mean <- mean(unlist(CD_A1L))
CD_A1L_sd <- Stand_Error(unlist(CD_A1L))
MD_A1L_mean <- mean(unlist(MD_A1L))
MD_A1L_sd <- Stand_Error(unlist(MD_A1L))
VU_A1L_mean <- mean(unlist(VU_A1L))
VU_A1L_sd <- Stand_Error(unlist(VU_A1L))
WS_A1L_mean <- mean(unlist(WS_A1L))
WS_A1L_sd <- Stand_Error(unlist(WS_A1L))
UWS_A1L_mean <- mean(unlist(UWS_A1L))
UWS_A1L_sd <- Stand_Error(unlist(UWS_A1L))

save(speic_list_A1L,AD_A1L,AP_A1L,CO_A1L,NE_A1L,TL_A1L,ED_A1L,TA_A1L,NC_A1L,CB_A1L,CD_A1L,MD_A1L,VU_A1L,WS_A1L,UWS_A1L,AD_A1L_mean,AD_A1L_sd,AP_A1L_mean,AP_A1L_sd,CO_A1L_mean,CO_A1L_sd,NE_A1L_mean,NE_A1L_sd,TL_A1L_mean,TL_A1L_sd,ED_A1L_mean,ED_A1L_sd,TA_A1L_mean,TA_A1L_sd,NC_A1L_mean,NC_A1L_sd,CB_A1L_mean,CB_A1L_sd,CD_A1L_mean,CD_A1L_sd,MD_A1L_mean,MD_A1L_sd,VU_A1L_mean,VU_A1L_sd,WS_A1L_mean,WS_A1L_sd,UWS_A1L_mean,UWS_A1L_sd,file="A1F.Rdata")


