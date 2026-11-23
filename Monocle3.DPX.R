rm(list=ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
#library(presto)
library(ggplot2)
SEURATtotal<- readRDS(file = "./fibroblast.rds")


P2 <- DimPlot(DA, 
              reduction = "umap", 
              label = TRUE, 
              group.by = "seurat_clusters", 
              pt.size = 1.2, 
              label.size = 6) + 
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2, 
                   arrow = arrow(length = unit(0.2, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        legend.text = element_text(size = 14),  # 调整图例字体大小
        legend.title = element_text(size = 16)) # 调整图例标题字体大小


P2 <- DimPlot(group, 
              reduction = "umap", 
              label = TRUE, 
              group.by = "group", 
              pt.size = 1.2, 
              label.size = 6) + 
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2, 
                   arrow = arrow(length = unit(0.2, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        legend.text = element_text(size = 14),  # 调整图例字体大小
        legend.title = element_text(size = 16)) # 调整图例标题字体大小
# 提取Early stage的细胞
early_stage_cells <- subset(DA, stage == "Early stage")
ctrl <- subset(DA, group == "ctrl")
PD <- subset(DA, group == "PD/PDD")
# 提取“Early stage”和“Late stage”的细胞
stage <- subset(DA, stage %in% c("Early stage", "Late stage"))
group <- subset(DA, group %in% c("ctrl", "PD/PDD"))
Late_stage_cells <- subset(DA, stage == "Late stage")
library(SCP)
color <- c("gray", "coral2")#设置颜色  
plot8 <- FeaturePlot(ctrl, features = 'SLC44A1', cols = color, pt.size = 1) +  
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) +
  ggtitle("SLC44A1 in ctrl")  # 添加标题

# 显示图形
print(plot8)


plot9 <- FeaturePlot(PD, features = 'SLC44A1',cols = color, pt.size = 1)+  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) +
  ggtitle("SLC44A1 in PD")#加边框
combined_plot <- plot8 + plot9 + plot_layout(ncol = 2)

plot3 <- FeaturePlot(early_stage_cells, features = 'SLC44A1', cols = color, pt.size = 1) +  
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) +
  ggtitle("SLC44A1 in Early stage")  # 添加标题

plot4 <- FeaturePlot(Late_stage_cells, features = 'SLC44A1',cols = color, pt.size = 1)+  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) +
  ggtitle("SLC44A1 in Late stage")#加边框

combined_plot1 <- plot3 + plot4 + plot_layout(ncol = 2)
# 显示合并后的图形
print(combined_plot1)
# 保存为PDF文件
ggsave("G:/邓鹏熙/result/umap.stage.slc44a1.pdf", plot = combined_plot1, device = "pdf", width = 12, height = 6)


SEURATtotal <- SCTransform(SEURATtotal)
SEURATtotal <- RunPCA(SEURATtotal, features = VariableFeatures(object = SEURATtotal))
SEURATtotal <- RunUMAP(SEURATtotal, dims = 1:10)
SEURATtotal <- FindNeighbors(SEURATtotal, dims = 1:10)
SEURATtotal <- FindClusters(SEURATtotal, resolution = 0.2 )
SEURATtotal <- RunHarmony(SEURATtotal, group.by.vars="orig.ident",assay.use="SCT",max.iter.harmony=20)
SEURATtotal <- JoinLayers(SEURATtotal)
SEURATtotal.markers <- FindAllMarkers(SEURATtotal, only.pos = TRUE)
SEURATtotal.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(SEURATtotal.markers, file = "neuron.markers.csv", row.names = T)
SEURATtotal.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
plot3 <- DoHeatmap(SEURATtotal, features = top10$gene) + NoLegend()

ggsave("neuron_heatmap.pdf", plot3, width = 8, height = 6)
saveRDS(SEURATtotal,"neuron.RDS")

setwd("G:/邓鹏熙")
DA<- readRDS(file = "./DA.rds")

load("G:/邓鹏熙/cds.rData")
my_comparisons <- list(c("Early stage", "Late stage"))
my_comparisons <- list(c("ctrl", "PD/PDD"))
p1 <- VlnPlot(DA,features='SLC44A1',
              group.by="group"
              )
library(ggplot2)
library(ggpubr)
Idents(group)<- group$group
Idents(stage)<- stage$stage
# 绘制小提琴图
VlnPlot(group, features = 'SLC44A1', add.noise = FALSE) +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(color = 'black', face = "bold", size = 12), 
        axis.text.y = element_text(color = 'black', face = "bold"),       
        axis.title.y = element_text(color = 'black', face = "bold", size = 15),     
        panel.grid.major = element_blank(),     
        panel.grid.minor = element_blank(),      
        panel.border = element_rect(color = "black", size = 1.2, linetype = "solid"),     
        panel.spacing = unit(0.12, "cm"),      
        plot.title = element_text(hjust = 0.5, face = "bold.italic"),    
        legend.position = 'none') +  
  stat_compare_means(method = "t.test", hide.ns = FALSE,   
                     comparisons = my_comparisons,                   
                     label = "p.signif",                   
                     bracket.size = 0.8,                  
                     tip.length = 0, 
                     size = 6) +  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

library(monocle3)
#library(monocle)
library(dplyr)
library(Seurat)
library(patchwork)
data <- GetAssayData(DA, assay = 'RNA', layer = 'counts')

#其次获得meta信息
cell_metadata <- DA@meta.data
cell_metadata
#然后要获得gene_annotation 
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
table(DA@meta.data$stage)
#创建CDS对象
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 30)     #preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
plot_pc_variance_explained(cds) 
#umap降维
cds <- reduce_dimension(cds,preprocess_method = "PCA") #preprocess_method默认是PCA

#tSNE降维
cds <- reduce_dimension(cds, reduction_method="tSNE")
cds <- reduce_dimension(cds, method = "UMAP")
cds <- cluster_cells(cds, reduction_method = 'UMAP', cluster_method = "louvain")
plot_cells(cds, reduction_method="UMAP", color_cells_by="stage") 
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=F,
           group_label_size=6,
           cell_size=1.5)  


# 提取 "late stage" 和 "early stage" 的细胞
selected_cells <- subset(sc, subset = stage %in% c("late stage", "early stage"))
table(sc@meta.data$group)
#monocle2教程
#构建矩阵
expr_matrix <- as(as.matrix(GetAssayData(sc, assay = "RNA", slot = "counts")), 'sparseMatrix')
p_data <- sc@meta.data 
p_data$celltype <- sc@active.ident #表型信息
#基因信息
f_data <- data.frame(gene_short_name = row.names(sc),row.names = row.names(sc))
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
#构建cds对象，参数可参考默认设置
# 确保 expr_matrix 和 fd 的行名一致
fData <- data.frame(gene_short_name = rownames(expr_matrix), row.names = rownames(expr_matrix))
fd <- new("AnnotatedDataFrame", data = fData)

# 创建 CellDataSet
cds <- monocle::newCellDataSet(expr_matrix,
                               phenoData = pd,
                               featureData = fd,
                               lowerDetectionLimit = 0.5,
                               expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
#降维
cds <- monocle::reduceDimension(cds, max_components = 2,method = 'DDRTree')


AFD_genes <- c("SLC44A1","PINK1", "PARKIN", "FUNDC1", "BNIP3","BNIP3L")
AFD_genes <-"SLC44A1"
AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,
                       colData(cds)$stage %in% c("Early stage")]
AFD_lineage_cds <- order_cells(AFD_lineage_cds)
plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="pseudotime",
                         min_expr=0.5)


plot_cells(cds,
           genes=AFD_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
library(SCP)
ht<- readRDS(file = "./rename.rds")


celltype <-          c('0' = 'Neuron', 
                       '1' = 'Oligodendrocytes', 
                       '2' = 'Astrocytes', 
                       '3' = 'Neuron', 
                       '4' = 'Microglia', 
                       '5' = 'Neuron',
                       '6' = 'Neuron',
                       "7"="Endothelial cells",
                       "8"="Others","9"="Others","10"="Others","11"="Others")
ht[['celltype']] = base::unname(celltype[ht@meta.data$seurat_clusters])
P2 <- DimPlot(ht, 
              reduction = "umap", 
              label = TRUE, 
              group.by = "celltype", 
              pt.size = 1.2, 
              label.size = 6) + 
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2, 
                   arrow = arrow(length = unit(0.2, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        legend.text = element_text(size = 14),  # 调整图例字体大小
        legend.title = element_text(size = 16)) # 调整图例标题字体大小

p <- GroupHeatmap(
  srt = ht,
  features = c(
    "DCC", "RBFOX3", "GAD2","SLC17A6", "SLC6A3","TH","SLC18A2",#Neuron 
    "GFAP","AQP4", #Astrocytes
    "MOBP", "MAG","UGT8", # Oligodendrocytes
    "Fev", "Neurod1", # Pre-endocrine
    "Rbp4", "Pyy", # Endocrine
    "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
  ),
  group.by = "celltype",
  heatmap_palette = "YlOrRd",
  #cell_annotation = c("Phase", "G2M_score", "Cdh2"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE
)
print(ht$plot)
library(scop)
