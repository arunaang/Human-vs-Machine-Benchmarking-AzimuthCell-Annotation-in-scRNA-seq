### Single-Cell RNA seq Analysis + Azimuth + GSEA  #####

#Install Packages======
install.packages("dplyr")
install.packages("Seurat")
install.packages("readr")
install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("fgsea")
install.packages("ggtext")
install.packages("here")

#==================================================
# 1. Load Libraries & Data
#==================================================

# Load Libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(Azimuth)
library(readr)
library(ggplot2)
library(ggtext)
library(fgsea)
library(tidyverse)
# Load Data 
data2 <- Read10X(data.dir = "data/hg19")
# Print the first 6 rows by default
head(data2)
# no of rows vs no of colums 
dim(data2)
# Check class
class(data2)
dim(data2)
str(data2)
# Enable automatic calculation
options(Seurat.object.assay.calcn = TRUE)


# ==================================================
# 2. Create Seurat Object
# ==================================================

# Create a Seurat Object
data_seurat2 <- CreateSeuratObject(counts = data2, min.cells = 3 , min.features = 200)
# Check metadata again
head(data_seurat2@meta.data)
#changing row names and col names
gene_names<-rownames(data_seurat2)
cell<- colnames(data_seurat2)
# Calculate mitochondrial gene percentage
grep('^MT-',gene_names)
data_seurat2[["percent.mt"]] <- PercentageFeatureSet(data_seurat2, pattern = "^MT-")
data_seurat2@meta.data
#Visualizing all 3 cols 
VlnPlot(data_seurat2, c("nFeature_RNA","nCount_RNA","percent.mt"))


# ==================================================
# 3. Quality Control
# ==================================================

# QC 
data_sub2<- subset(data_seurat2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) 
VlnPlot(data_sub2, c("nFeature_RNA","nCount_RNA","percent.mt"))   
ggsave("figures/Vlnplot_subset.png", bg = "white")
dim(data_sub2)


# ==================================================
# 4. Normalization & PCA
# ==================================================

# Normalization
data_normalize<-NormalizeData(data_sub2)
# Finding variable genes!
data_var2<-FindVariableFeatures(data_normalize, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(data_var2)
ggsave("figures/Variable_feature_plot.png", bg = "white")

# ==================================================
# 5. Clustering
# ==================================================

#Scaling data
data_scale<-ScaleData(data_var2)
#Performing PCA
pca<- RunPCA(data_scale)
#visualizing the PC1 vs PC 50 ( max variance)
DimPlot(pca, dims = c(1,2), reduction = "pca")
ggsave("figures/PCA.png",plot = last_plot(), bg = "white")
#visualizing the PC which has most variance
ElbowPlot(pca)
ggsave("figures/Elbow_Plot.png",plot = last_plot(), bg = "white")
# Clustering the cells 
pbmc_data<- FindNeighbors(pca, dims = 1:8)
pbmc_data <- FindClusters(pbmc_data, resolution = 0.5)
umap_data2<-RunUMAP(pbmc_data, dims = 1:8)
DimPlot(umap_data2, reduction = "umap", group.by = "seurat_clusters")
ggsave("figures/Seurat_Clusters.png",plot = last_plot(), bg = "white")
# ==================================================
# 6. Marker Genes
# ==================================================

# Finding the markers
pbmc_marker<-FindAllMarkers(pbmc_data,only.pos = TRUE) # helps to find differencially expressed genes in one cluster vs the other clusters
#Note: #Findconservedmarkers()- if there are 2 groups and 2 condition and this will help separate 2 types of cell band based on the condition.
class(pbmc_marker) 
head(pbmc_marker)

pbmc_marker_top <-pbmc_marker%>%
  group_by(cluster)%>%
  filter(avg_log2FC >= 2, p_val_adj < 0.05)%>%
  slice_min(p_val_adj, n = 5)


# Print
print(pbmc_marker_top)
#Visualizing it using dotplot 
marker_genes2<- pbmc_marker_top$gene
marker_genes2_unique <- unique(marker_genes2)
gene_list_by_cluster <- split(pbmc_marker_top$gene, pbmc_marker_top$cluster)
DotPlot(pbmc_data, features = marker_genes2_unique) +                             
  scale_color_gradient(low = "black", high = "red") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 1, size = 8),
        axis.text.y = element_text(size = 10))+
  ylab("Clusters")

ggsave("figures/Differencially expressed genes in each Seurat Clusters.png",plot = last_plot(), bg = "white", width = 12, height = 7)
# Naming clusters based on gene expressed- through my findings 
DotPlot(pbmc_data, features = c("LDHB", "S100AB", "IL32", "CD79A", "CCL5","GZMB", "GP9"))+ RotatedAxis()

# Naming clusters based on most common PBMC cell types (genes predicted by google)  
DotPlot(pbmc_data, features = c("TNFRSF4", "CLEC4E", "LYZ", "GNLY", "FCGR3A","GP9")) + RotatedAxis()

# ==================================================
# 7. Manual Annotations
# ==================================================

# Assigning celltypes to clustrers:
new.cluster.ids <- c(
  "CD4_Memory", 
  "CD4_Naive", 
  "CD4_Effector", 
  "B_Cell", 
  "CD4_EffectorMemory", 
  "Monocyte", 
  "CD8_NK_Cytotoxic", 
  "Monocyte", 
  "Platelet")
names(new.cluster.ids) <- levels(umap_data2)
pbmc_rename <- RenameIdents(umap_data2, new.cluster.ids)
my_annots<- DimPlot(pbmc_rename, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position = "right")+
  ggtitle("Manual Annotations") +
  theme(
    legend.position = "none",
    
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
ggsave("figures/Assigning Celltypes to Clustrers.png",plot = last_plot(), bg = "white")





# ==================================================
# 8. Azimuth Comparison
# ==================================================

## AZIMUTH
data_azimu<- RunAzimuth(query =data_seurat2  ,reference = "pbmcref" )
# Plottig all 3 variations of Azimuth reference annotations
pi<-DimPlot(data_azimu, group.by = "predicted.celltype.l1", label = TRUE, label.size = 3 )+ NoLegend()
pii<-DimPlot(data_azimu, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 )+ NoLegend()
piii<-DimPlot(data_azimu, group.by = "predicted.celltype.l3", label = TRUE, label.size = 3 )+ NoLegend()
pi+pii+piii
ggsave("figures/azimuth_reference_fin.png",plot = last_plot(), height = 10, width = 22 )
# Prediction score based clustering
Idents(data_azimu)<- "predicted.celltype.l2"
FeaturePlot(data_azimu, features = "predicted.celltype.l2.score", label = TRUE, repel = TRUE)+ ggtitle("Umap based on predicted score and celltype" )
#Extract metadata
meta_data<-data_azimu@meta.data
# Create a data frame with just Cell Id and predicted.celltype.l2
export_data <- data.frame(Cell_ID = rownames(meta_data), Predicted_Celltype_l2 = meta_data$predicted.celltype.l2)
# Export to CSV
write.csv(export_data, "predicted_celltype_L2.csv", row.names = FALSE)
# Comparing azimuth annots to my UMAP
df_new<- pbmc_rename@meta.data 
df_new$umi<-rownames(df_new)
meta_data$umi<- rownames(meta_data)
df_new<-df_new%>%
  left_join(meta_data, by = "umi" )
dim(df_new)  
pbmc_rename@meta.data$azimuth_predicts<- df_new$predicted.celltype.l1
dim(pbmc_rename)
colnames(pbmc_rename@meta.data)
dim(data_azimu)
# Comparing my annots with Azimuth annots 
azimu_annots<- DimPlot(pbmc_rename, group.by = "azimuth_predicts", label = TRUE, pt.size = 0.5) +
  ggtitle("Azimuth_annotation") +
  theme( legend.position = "none",
         plot.title = element_text(hjust = 0.5, face = "bold")
  )
ggsave("figures/Assigning Celltypes to Clustrers.png",plot = last_plot(), bg = "white")
my_annots + azimu_annots
ggsave("figures/Comparison between Manual and Azimuth Annotations.png")


# ==================================================
# 9. GSEA Analysis
# ==================================================

#GSEA Analysis 
##prep data to load 

pbmc_marker_cluster0 <-FindMarkers(pbmc_data, ident.1 = 0, only.pos = TRUE) # helps to find differencially expressed genes in one cluster vs the other clusters
df_gsea_2 <- pbmc_marker_cluster0 %>% arrange(desc(avg_log2FC))
df_gsea_2$gene<- rownames(df_gsea_2)
df_gsea_subset2<- df_gsea_2[, c("gene", "avg_log2FC")]

#write as rnk file
write.table(
  df_gsea_subset2,
  file = "Results/cluster0_vs_others_gsea_fin.rnk",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
##Load tsv results from pre-ranked gsea
positive_na<-read_tsv("Results/gsea_report_for_na_pos_1763315471064.tsv")
positive_na_subset<- positive_na[1:10,]
negative_na<- read_tsv("Results/gsea_report_for_na_neg_1763315471064.tsv")
negative_na_subset<- negative_na[12:21,]
pos_neg_data<-rbind(positive_na_subset, negative_na_subset)
pos_neg_data$NES= as.numeric(pos_neg_data$NES)
pos_neg_data
##Coloring the text acc to the color of values 
# Remove NA FDR rows
pos_neg_data <- pos_neg_data[!is.na(pos_neg_data$`NES`), ]

# Make a safe color function
col_fun <- scales::col_numeric(
  palette = c("blue", "red"),
  domain = range(pos_neg_data$`NES`, na.rm = TRUE)
)

# Generate colors
pos_neg_data$color <- col_fun(pos_neg_data$`NES`)

# Fix any NA colors
pos_neg_data$color[is.na(pos_neg_data$color)] <- "black"

# Markdown labels
pos_neg_data$Pathways <- paste0(
  "<span style='color:", pos_neg_data$color, "'>",
  pos_neg_data$NAME,
  "</span>"
)

# Plot
ggplot(pos_neg_data, aes(x = NES, y = Pathways)) +
  geom_point(aes(size = SIZE, color = `FDR q-val`)) +
  scale_color_gradient(low = "hotpink", high = "cyan") +
  theme(axis.text.y = element_markdown())

ggsave("figures/GSEA analysis plot cluster_0 vs others.png",plot = last_plot(), bg = "white")