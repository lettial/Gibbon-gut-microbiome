#  Load essential packages
library(randomForest)

# Customize function
Stand_Error <- function(x){
  sd(x)/length(x)
}

maefun <- function(pred, obs) mean(abs(pred - obs))  

msefun <- function(pred, obs) mean(abs(pred - obs)^2)  

nmsefun <- function(pred, obs) mean((pred - obs)^2)/mean((mean(obs) - obs)^2)

# read the file
otu_table <- read.table("core_otu.txt",sep = "\t",row.names = 1,header = T)
diet <- read.table('diet.txt',header = T,sep = "\t",row.names = 1)
taxon_select <- as.character("Proteobacteria") #changed according to your case
subject_select <- as.character(A)
diet_select <- subset(diet, Subject == subject_select)
otu_table_select <- otu_table[rownames(diet_select),]
diet_select <- diet_select[,c(6:10)]
set.seed(123)
mae_list <- list()
mse_list <- list()
nmse_list <- list()
mae_null_list <- list()
mse_null_list <- list()
nmse_null_list <- list()
predict_data_frame <- data.frame()

# permutation
for (i in 1:1000){
  ind=sample(2,nrow(diet_select),replace=TRUE, prob=c(0.7,0.3)) # 70% training, 30% validation
  genus_rm = randomForest(diet_select[ind == 1,], otu_table_select[ind== 1,][,taxon_select], importance=TRUE, proximity=TRUE, ntree = 1000)
  # Null model
  table_null <- data.frame(cbind(rownames(otu_table_select), otu_table_select[,taxon_select]))
  table_null$taxon_select_2 <- sample(table_null$X2)
  rownames(table_null) <- table_null$X1
  genus_rm_null = randomForest(diet_select[ind == 1,], as.numeric(table_null[ind== 1,]$taxon_select_2), importance=TRUE, proximity=TRUE, ntree = 1000)
  train = diet_select[ind==1,]
  test = diet_select[ind==2,]
  plant_predict <- predict(genus_rm, test)
  plant_predict_null <- predict(genus_rm_null, test)
  mae_list[i] = maefun(plant_predict, otu_table_select[ind== 2,][,taxon_select])  
  mse_list[i] = msefun(plant_predict, otu_table_select[ind== 2,][,taxon_select]) 
  nmse_list[i] = nmsefun(plant_predict, otu_table_select[ind== 2,][,taxon_select])
  mae_null_list[i] = maefun(plant_predict_null, otu_table_select[ind== 2,][,taxon_select]) 
  mse_null_list[i] = msefun(plant_predict_null, otu_table_select[ind== 2,][,taxon_select])
  nmse_null_list[i] = nmsefun(plant_predict_null, otu_table_select[ind== 2,][,taxon_select])
  predict_data <- data.frame(plant_predict)
  for (sample_name in rownames(predict_data)){
    predict_data_frame[i, sample_name] <- predict_data[sample_name,]
}   
  print(paste(subject_select,i,sep="-"))
}
  Mean_Table <- data.frame(colMeans(predict_data_frame,na.rm=T))
  Mean_Table_order <- data.frame(Mean_Table[order(rownames(Mean_Table)),])
  colnames(Mean_Table_order) <- taxon_select
  write.table(Mean_Table_order, "Mean.txt", sep="\t", quote=F, row.names=F)
  
  # make sure the value is valid.
  if (T %in% sapply(mae_list, is.infinite)){
    is.na(mae_list) <- sapply(mae_list, is.infinite)
}

  if (T %in% sapply(mse_list, is.infinite)){
    is.na(mse_list) <- sapply(mse_list, is.infinite)
}

  if (T %in% sapply(nmse_list, is.infinite)){
    is.na(nmse_list) <- sapply(nmse_list, is.infinite)
}

  if (T %in% sapply(mae_null_list, is.infinite)){
    is.na(mae_null_list) <- sapply(mae_null_list, is.infinite)
}

  if (T %in% sapply(mse_null_list, is.infinite)){
    is.na(mse_null_list) <- sapply(mse_null_list, is.infinite)
}

  if (T %in% sapply(nmse_null_list, is.infinite)){
    is.na(nmse_null_list) <- sapply(nmse_null_list, is.infinite)
}
  print(unlist(mae_list))
  print(unlist(mae_null_list))
  cat(subject_select, taxon_select, mean(unlist(mae_list), na.rm=T), sd(unlist(mae_list), na.rm=T), Stand_Error(unlist(mae_list)), mean(unlist(mse_list), na.rm=T), sd(unlist(mse_list), na.rm=T), Stand_Error(unlist(mse_list)), mean(unlist(nmse_list), na.rm=T), sd(unlist(nmse_list), na.rm=T), Stand_Error(unlist(nmse_list)), subject_select, taxon_select, "Null", mean(unlist(mae_null_list), na.rm=T), sd(unlist(mae_null_list), na.rm=T), Stand_Error(unlist(mae_null_list)), mean(unlist(mse_null_list), na.rm=T), sd(unlist(mse_null_list), na.rm=T), Stand_Error(unlist(mse_null_list)), mean(unlist(nmse_null_list), na.rm=T), sd(unlist(nmse_null_list), na.rm=T), Stand_Error(unlist(nmse_null_list)), subject_select, taxon_select,wilcox.test(unlist(mae_list), unlist(mae_null_list))$p.value,wilcox.test(unlist(mse_list), unlist(mse_null_list))$p.value, wilcox.test(unlist(nmse_list), unlist(nmse_null_list))$p.value,"\n", file="statistical_result.txt", append=T)


save(mae_list,mse_list,nmse_list,mae_null_list,mse_null_list,nmse_null_list,predict_data_frame,file="final_result.Rdata", sep="")

