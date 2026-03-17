set.seed(1234)

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(glmGamPoi)
library(presto)
library(spacedeconv)
library(SpatialExperiment)

plot_dir <- "spatial_plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(plot_obj, filename, width = 8, height = 6, dpi = 300) {
  ggplot2::ggsave(
    filename = file.path(plot_dir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi
  )
}

brca_sp <- readRDS("Data/brca_spatial_s4.rds")
print(brca_sp)

p <- SpatialPlot(brca_sp, alpha = 0)
save_plot(p, "01_spatial_plot.png")

head(brca_sp@meta.data)

p <- VlnPlot(
  brca_sp,
  features = c("nCount_Spatial", "nFeature_Spatial"),
  pt.size = 0.1
) + NoLegend()
save_plot(p, "02_qc_violin.png")

p <- SpatialFeaturePlot(
  brca_sp,
  features = "nCount_Spatial",
  image.alpha = 0
) +
  theme(legend.position = "right")
save_plot(p, "03_ncount_spatial.png")

p <- SpatialFeaturePlot(
  brca_sp,
  features = "nCount_Spatial",
  image.alpha = 0
) +
  theme(legend.position = "right")
save_plot(p, "04_nfeature_spatial_placeholder.png")

brca_sp <- SCTransform(
  brca_sp,
  assay = "Spatial",
  verbose = FALSE
)

p <- SpatialFeaturePlot(brca_sp, features = c("KRT8", "COL1A2", "CD8B"))
save_plot(p, "05_marker_triplet.png", width = 10, height = 4)

p <- SpatialFeaturePlot(brca_sp, features = "KRT8")
save_plot(p, "06_krt8_default.png")

p <- SpatialFeaturePlot(brca_sp, features = "KRT8", pt.size.factor = 1)
save_plot(p, "07_krt8_ptsize1.png")

p <- SpatialFeaturePlot(brca_sp, features = "KRT8", alpha = c(0.1, 1))
save_plot(p, "08_krt8_alpha.png")

brca_sp <- RunPCA(
  brca_sp,
  assay = "SCT",
  verbose = FALSE
)
brca_sp <- FindNeighbors(
  brca_sp,
  reduction = "pca",
  dims = 1:30
)
brca_sp <- FindClusters(
  brca_sp,
  resolution = 0.5,
  verbose = FALSE
)
brca_sp <- RunUMAP(
  brca_sp,
  reduction = "pca",
  dims = 1:30
)
p <- DimPlot(
  brca_sp,
  reduction = "umap",
  label = TRUE
)
save_plot(p, "09_umap_clusters.png")

p <- SpatialDimPlot(brca_sp, label = TRUE, label.size = 3)
save_plot(p, "10_spatial_clusters.png")

p <- SpatialDimPlot(
  brca_sp,
  cells.highlight = CellsByIdentities(object = brca_sp, idents = c(1, 4)),
  facet.highlight = TRUE,
  ncol = 2
)
save_plot(p, "11_spatial_highlight_clusters.png", width = 10, height = 5)

de_markers <- FindMarkers(brca_sp, ident.1 = 1, ident.2 = 4)
head(de_markers)
p <- SpatialFeaturePlot(
  object = brca_sp,
  features = rownames(de_markers)[1:4],
  alpha = c(0.1, 1),
  ncol = 2
)
save_plot(p, "12_de_markers_top4.png", width = 10, height = 8)

brca_sc <- readRDS("Data/brca_sc_s4.rds")

print(brca_sc)
head(brca_sc@meta.data)

cat(
  "Coarse resolution:",
  length(unique(brca_sc@meta.data$celltype_coarse)),
  "labels\n"
)
print(unique(brca_sc@meta.data$celltype_coarse))
cat(
  "Medium resolution",
  length(unique(brca_sc@meta.data$celltype)),
  "labels\n"
)
print(unique(brca_sc@meta.data$celltype))
cat(
  "Fine-grained resolution",
  length(unique(brca_sc@meta.data$celltype_fine)),
  "labels\n"
)
print(unique(brca_sc@meta.data$celltype_fine))

brca_sc <- NormalizeData(brca_sc)
brca_sc <- FindVariableFeatures(brca_sc, nfeatures = 2000)
brca_sc <- ScaleData(brca_sc, features = rownames(brca_sc))
brca_sc <- RunPCA(brca_sc)
brca_sc <- FindNeighbors(brca_sc, dims = 1:40)
brca_sc <- FindClusters(
  brca_sc,
  resolution = 0.5,
  random.seed = 1234
)
brca_sc <- RunUMAP(
  brca_sc,
  dims = 1:40,
  n.neighbors = 30,
  min.dist = 0.3,
  seed.use = 1234
)
brca_sp <- SCTransform(brca_sp, assay = "Spatial")
brca_sc <- SCTransform(brca_sc)
anchors <- FindTransferAnchors(
  reference = brca_sc,
  query = brca_sp,
  normalization.method = "SCT",
  reference.assay = "SCT",
  query.assay = "SCT"
)
brca_sp <- TransferData(
  anchorset = anchors,
  query = brca_sp,
  slot = "data",
  refdata = brca_sc@meta.data$celltype_coarse,
  prediction.assay = TRUE
)
DefaultAssay(brca_sp) <- "prediction.score.id"

p <- SpatialFeaturePlot(
  brca_sp,
  features = c(
    "Cancer Epithelial",
    "Endothelial",
    "Myeloid",
    "CAFs",
    "T-cells",
    "B-cells"
  ),
  ncol = 2
)
save_plot(p, "13_label_transfer_scores.png", width = 10, height = 8)

brca_sce <- as.SingleCellExperiment(brca_sc)

brca_spe <- readRDS("Data/brca_spe.rds")

brca_sce <- spacedeconv::preprocess(brca_sce)
brca_spe <- spacedeconv::preprocess(brca_spe)
brca_spe <- spacedeconv::normalize(brca_spe, method = "cpm")

brca_sdwls_sig <- spacedeconv::build_model(
  single_cell_obj = brca_sce,
  cell_type_col = "celltype_coarse",
  method = "spatialdwls",
  verbose = TRUE
)

brca_sdwls_res <- spacedeconv::deconvolute(
  spatial_obj = brca_spe,
  method = "spatialdwls",
  assay_sp = "cpm",
  signature = brca_sdwls_sig
)

p <- plot_spatial(
  spe = brca_sdwls_res,
  result = "spatialdwls_CAFs",
  sample_id = "brca",
  image_id = "lowres",
  title = "CAF",
  density = FALSE
)
save_plot(p, "14_spatialdwls_cafs.png")

p <- plot_spatial(
  spe = brca_sdwls_res,
  result = "spatialdwls_CAFs",
  sample_id = "brca",
  image_id = "lowres",
  title = "CAF",
  smooth = TRUE,
  density = FALSE
)
save_plot(p, "15_spatialdwls_cafs_smooth.png")

p <- plot_spatial(
  spe = brca_sdwls_res,
  result = "spatialdwls_Myeloid",
  sample_id = "brca",
  image_id = "lowres",
  title = "Myeloid",
  smooth = TRUE,
  density = FALSE
)
save_plot(p, "16_spatialdwls_myeloid_smooth.png")

p <- plot_spatial(
  spe = brca_sdwls_res,
  result = "spatialdwls_T.cells",
  sample_id = "brca",
  image_id = "lowres",
  title = "T cells",
  smooth = TRUE,
  density = FALSE
)
save_plot(p, "17_spatialdwls_tcells_smooth.png")

p <- spacedeconv::plot_comparison(
  brca_sdwls_res,
  cell_type_1 = "spatialdwls_Myeloid",
  cell_type_2 = "spatialdwls_CAFs",
  sample_id = "brca",
  image_id = "lowres",
  density = FALSE,
  smooth = TRUE,
  title = "Myeloid vs. CAFs",
  palette = "Purple-Green",
  reverse_palette = TRUE,
  shift_positive = FALSE,
  nDigits = 2
)
save_plot(p, "18_myeloid_vs_cafs.png")

p <- plot_most_abundant(
  brca_sdwls_res,
  method = "spatialdwls",
  sample_id = "brca",
  image_id = "lowres",
  title = "Most abundant",
  min_abundance = 0.05,
  palette = "Accent"
)
save_plot(p, "19_most_abundant.png")

tf_ref <- get_decoupleR_reference(
  method = "collectri",
  organism = "human"
)
brca_sdwls_res <- compute_activities(
  brca_sdwls_res,
  reference = tf_ref,
  method = "wmean",
  assay = "cpm"
)

pw_ref <- get_decoupleR_reference(
  method = "progeny",
  organism = "human"
)
brca_sdwls_res <- compute_activities(
  brca_sdwls_res,
  reference = pw_ref,
  method = "wmean",
  assay = "cpm"
)
p <- plot_spatial(
  spe = brca_sdwls_res,
  result = "progeny_JAK.STAT",
  sample_id = "brca",
  image_id = "lowres",
  title = "JAK-STAT",
  smooth = TRUE,
  density = FALSE
)
save_plot(p, "20_progeny_jak_stat.png")

p <- plot_spatial(
  spe = brca_sdwls_res,
  result = "progeny_TGFb",
  sample_id = "brca",
  image_id = "lowres",
  title = "TGFb",
  smooth = TRUE,
  density = FALSE
)
save_plot(p, "21_progeny_tgfb.png")
