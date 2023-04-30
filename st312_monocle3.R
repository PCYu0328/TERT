library(Seurat)
library(monocle3)
library(ggplot2)
library(patchwork)
#get counts matrix from Seurat Object#
data <- GetAssayData(BTC_1, assay = 'SCT', slot = 'counts')
cell_metadata <- BTC_1@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds)

#clustering#
cds = cluster_cells(cds, resolution=1e-5)

#Interactive selection of thyroid spots#
cds_subset <- choose_cells(cds)
plot_cells(cds_subset)
cds_subset = cluster_cells(cds_subset, resolution=1e-5)
cds_subset <- learn_graph(cds_subset)

#indicate the well-differentiated cells as the root#
cds_subset = order_cells(cds_subset)

plot_cells(cds_subset,
           color_cells_by = "pseudotime",
           trajectory_graph_color = "#1CCCFF",
           trajectory_graph_segment_size = 3,
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1.5)

