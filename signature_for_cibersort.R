library(Seurat)
scthca = readRDS('THCA.rds')

exp_matrix = scthca@assays$RNA@counts#expression matrix
scthca@meta.data$celltype = Idents(scthca)
metadata = scthca@meta.data#metadata
identical(row.names(scthca@meta.data),colnames(exp_matrix))#TRUE

set.seed(20230401)
##sampling
cellnames = as.character(levels(scthca$celltype))
sample_name = c()
nsample = 200

for (i in 1:length(cellnames)) {
  newmetadata = metadata[metadata$celltype == cellnames[i],]
  sample_name = c(sample_name,
                  row.names(newmetadata[sample(1:nrow(newmetadata),nsample,replace = F),]))
  sample_name = gsub("\\.[0-9]","",sample_name)
  rm(newmetadata)
}

cell_names = metadata[sample_name,]$celltype
exp_matrix = exp_matrix[,sample_name]
colnames(exp_matrix) = cell_names
exp_matrix = as.data.frame(exp_matrix)
exp_matrix = exp_matrix[rowSums(exp_matrix)>0,]

#save the signature matrix
write.table(exp_matrix,file = "exp_matrix.txt",sep = "\t",col.names = T,row.names = T)

##then we used CIBERSORTx website.(https://cibersortx.stanford.edu/) 
