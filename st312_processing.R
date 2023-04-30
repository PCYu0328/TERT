library(Seurat)
library(dplyr)
library(ggplot2)
library(hdf5r)
library(patchwork)
library(glmGamPoi)

#load spatial files#
BTC <- Load10X_Spatial(data.dir = "F:/yupengcheng_spaceranger_results.tar/312BTC_spatial",slice = "BTC")

#we used cloupe to recognize the blank region#
uninc <- read.csv('zone-1.csv',header = T)
uniID <- uninc[,1]
BTC_1 <- subset(BTC,cells = uniID,invert=T)

BTC_1 <- SCTransform(BTC_1, method = "glmGamPoi",assay ="Spatial",verbose = F)

DefaultAssay(BTC_1) <- "SCT"

BTC_1 <- RunPCA(BTC_1)
BTC_1 <- FindNeighbors(BTC_1,dims = 1:30)
BTC_1 <- FindClusters(BTC_1)

p1_1 <- SpatialDimPlot(BTC_1,label = T,label.size = 3)
p1_1

BTC_1 <- RunTSNE(BTC_1,dims = 1:30,label=T)
p2_1 <- DimPlot(BTC_1, reduction = "tsne",pt.size = 1.5,label = T,label.box = T,repel = T) + NoLegend()
p2_1

#find marker genes and annotate cell types#
DefaultAssay(BTC_1) <- "Spatial"
markers <- FindAllMarkers(object = BTC_1, test.use="wilcox" ,
                          only.pos = F,
                          logfc.threshold = 0.25)   
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
## Identify the cell type based on the DEGs of each cluster
#thyroid epithelial cells (Epcam, Tg, Tpo)
#salivary gland cells (Pigr, Msln, Agr2)
#adipocyte cells (Adipoq, Cidec, Cfd) 
#immune cells (Ptprcap, Cd79b, Ccl5, Cxcl9)
#muscle cells (Actn2, Trdn, Mypz1) 
#erythroid cells (Hbb-bt, Hba-a1, Hba-a2) 
#esophagus epithelial cells (Mt4, Lgals7, Krt78, Fam25c) 


#analyze thyroid spots#
thyro2 <- subset(BTC_1,ident = "Thyroid epithelial cells")
thyro2 <- SCTransform(thyro2, method = "glmGamPoi",assay ="Spatial",verbose = F)
thyro2 <- RunPCA(thyro2)
thyro2<- FindNeighbors(thyro2,dims = 1:30)

library(clustree)
thyro2 <- FindClusters(thyro2, resolution = seq(0.1,1,by=0.1),graph.name = "SCT_snn")
clustree(thyro2)
Idents(thyro2) <- thyro2@meta.data$SCT_snn_res.0.8
thyro2 <- RunTSNE(thyro2, dims = 1:30)


#c0, the well-differentiated cluster#
SpatialDimPlot(thyro2,
               pt.size.factor = 2.5,
               label = T, label.size = 4,label.color = "white",label.box = T,
               alpha = 0.8,
               cells.highlight = CellsByIdentities(thyro2,idents = c("c0")),
               cols.highlight = c("#F8766D","grey50"),
               images = NULL,
               facet.highlight = T,
)

#TDS score#
tds_gene <- list(c("Tg","Tpo","Dio1",
                    "Fhl1","Sorbs2"))
thyro2 <- AddModuleScore(
  object = thyro2,
  features = tds_gene,
  ctrl = 200,
  name = "TDS_Thyro2"
)

data<- FetchData(thyro2,vars = c("SCT_snn_res.0.8","TDS_Thyro21"))
data$SCT_snn_res.0.8 <- factor(data$SCT_snn_res.0.8,levels = c("c0","c1","c2","c3","c4","c5"))
ggplot(data, aes(x=SCT_snn_res.0.8,y=TDS_Thyro21)) +
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=10,hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL,title = "TDS score")+ geom_jitter(col="#00000033", pch=19,cex=2, position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0),aes(fill= SCT_snn_res.0.8),width=0.7,cex=0.8)+
  NoLegend()+theme(plot.title = element_text(hjust = 0.5))