# Aim: This script will generate a single Seurat object for all scRNA-seq datasets and also a 
# combined object for all scRNA-seq datasets 
# Input: count file from the cellranger pipeline. 
# Technical details: The scRNAs-seq datasets were generated using the 5' VDJ kit from 10x Genomics 
# and fastq files were generated using cellranger. 

# setting environment ------------------------------------------------------------------------------
# loading libraries and custom functions
pacman::p_load(tidyverse, edgeR, data.table, scales, gtools, Seurat, rio, DT, Matrix, Matrix.utils, monocle3, pheatmap)
source(file.path(my.path, "R", "utilities.R"))

# set data path
my.path = "/home/yanjun/Projects/Brain_organoids/Analysis/hBO/annotation_Fantom5/"


# Annotate clusters using cellbarcode idendity from previous Seurat object ------------------------
# Prepare old cell idendities for transfer to new Seurat objects, the new annotation should not change 
# much in aspect of protein coding genes, therefore we assume that the cell idendities are tranferrable. 

# load previous Seurat object
Seurat_all_old <- readRDS("/home/yanjun/Projects/Brain_organoids/Analysis/hBO/annotation_Fantom5/Data/hBO_Seurat_All.RDS")
Seurat_all_old_df <- as_tibble(Seurat_all_old@meta.data, rownames = "cellbarcode")

print("QC how many cells do you have per timepoint, check if that aligns with your current dataset")
print(Seurat_all_old_df %>% group_by(orig.ident) %>% count())
print("QC how many clusters do you have per timepoint (tally)")
print(Seurat_all_old_df %>% group_by(orig.ident, seurat_clusters) %>% count())
# In this case the old cluster number and cells/cluster are similar to the new objects

# pre-processing the annotation file --------------------------------------------------------------
annotation = read_tsv("/home/yanjun/Projects/Brain_organoids/Analysis/hBO/annotation_SCAFE_v3/data/pool_all/annotate/pool_all/bed/pool_all.CRE.annot.bed.gz", col_names = FALSE)
annotation = read_tsv("/home/yanjun/Projects/Brain_organoids/Analysis/hBO/annotation_SCAFE_v3/data/pool_all/annotate/pool_all/log/pool_all.CRE.info.tsv.gz")

# change underscore to hyphen (due to Seurat convention)
annotation$CREID <- sapply(annotation$CREID, gsub, 
                           pattern = "_", 
                           replacement = "-")




# generation of Seurat object for 40days timepoint ------------------------------------------------
# read raw data
dat.ref <- Read10X("/home/yanjun/Projects/Brain_organoids/Analysis/hBO/annotation_SCAFE_v3/data/pool_all/count/40days/matrix/")
colnames(dat.ref) <- paste0("40_Days_", colnames(dat.ref))

# create Seurat object
Seurat_40d <- CreateSeuratObject(counts = dat.ref)

# add gene names to data 
# Seurat_40d[["RNA"]]@meta.features$geneNameStr <- my_anno$geneNameStr

# getting rid of doublets, senescent cells, dead cells
Seurat_40d[["percent.mt"]] <- PercentageFeatureSet(Seurat_40d, pattern = "chrM")
VlnPlot(Seurat_40d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Remove dead cells
Seurat_40d <- subset(Seurat_40d, 
                     subset = nFeature_RNA > 800 & nFeature_RNA < 7500 & percent.mt < 20)

# Normalization
Seurat_40d <- NormalizeData(Seurat_40d)

# find differentially expressed genes 
Seurat_40d <- FindVariableFeatures(Seurat_40d)

# scale data
Seurat_40d <- ScaleData(Seurat_40d)

# find PCAs
Seurat_40d <- RunPCA(Seurat_40d,
                     features = VariableFeatures(object = Seurat_40d),
                     verbose = F)

# Heuristic method to identify how many PCAs to include
Seurat_40d <- JackStraw(Seurat_40d, dims = 50)
Seurat_40d <- ScoreJackStraw(Seurat_40d, dims = 1:50)
JackStrawPlot(Seurat_40d, dims = 1:50)
ElbowPlot(Seurat_40d)

# cluster cells
Seurat_40d <- FindNeighbors(Seurat_40d, dims = 1:50)
Seurat_40d <- FindClusters(Seurat_40d, resolution = 0.8)
Seurat_40d <- RunUMAP(Seurat_40d, dims = 1:50)


# DimPlot(Seurat_40d, reduction = "umap")

#extract cell metadata to combine with old cell type annotation
Seurat_40d_df <- as_tibble(Seurat_40d@meta.data, rownames = "cellbarcode")
Seurat_40d_df$cellbarcode <- strsplit2(Seurat_40d_df$cellbarcode, "-")[,1]
Seurat_40d_df = Seurat_40d_df %>% 
  left_join(Seurat_all_old_df %>% select(c("cellbarcode", "RNA_snn_res.0.8", "seurat_clusters")), 
            by = "cellbarcode",
            suffix = c("_tCRE", "_old"))

#QC how many cells are have no cell type annotation in the tCRE dataset
length(which(is.na(Seurat_40d_df$RNA_snn_res.0.8_old)))

#table showing how much the clustering has changes the cell type idendity
cell_type_tbl <- Seurat_40d_df %>% group_by(seurat_clusters_tCRE) %>% count(RNA_snn_res.0.8_old)
cell_type_tbl

#extracting keys with highest match
cell_type_tbl_short <- cell_type_tbl %>% group_by(seurat_clusters_tCRE) %>% top_n(1)
cell_type_tbl_short

# if cell annotation looks good, add to object
Seurat_40d@meta.data$RNA_snn_res.0.8_old <- Seurat_40d_df$RNA_snn_res.0.8_old
Seurat_40d@meta.data$seurat_clusters_old <- Seurat_40d_df$seurat_clusters_old

#update active ident in Seurat object
Seurat_40d <- SetIdent(Seurat_40d, value = Seurat_40d@meta.data$RNA_snn_res.0.8_old)
DimPlot(Seurat_40d, reduction = "umap")

# #generate vector according to old cell IDs
# current.cluster.ids <- c(0:7)
# new.cluster.ids <- c("40D_RGC_1", 
#                      "40D_PN_1", 
#                      "40D_NEC_1",
#                      "40D_NPC_1", 
#                      "40D_PN_2", 
#                      "40D_CP_1", 
#                      "40D_IP_1",
#                      "40D_NPC_2")

# obtain marker genes, takes time
Seurat_40d_markers = FindAllMarkers(Seurat_40d, only.pos = T)

# force to get all possible gene markers
Seurat_40d_markers_loose <- FindAllMarkers(Seurat_40d, 
                                           min.pct = 0,
                                           return.thresh = Inf,
                                           logfc.threshold = 0,
                                           min.cells.group = 0,
                                           only.pos = F)

# save data
write_csv(Seurat_40d_markers, "/home/yanjun/Projects/Brain_organoids/Analysis_v3_git/Brain_organoids/data/Seurat/Seurat_40d_markers.tsv")
write_csv(Seurat_40d_markers_loose, "/home/yanjun/Projects/Brain_organoids/Analysis_v3_git/Brain_organoids/data/Seurat/Seurat_40d_markers_loose.tsv")
saveRDS(Seurat_40d, "/home/yanjun/Projects/Brain_organoids/Analysis_v3_git/Brain_organoids/data/Seurat/40days.rds")









