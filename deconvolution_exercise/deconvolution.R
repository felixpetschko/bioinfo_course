## =========================================================
## Cell-type deconvolution analysis
## Guided analysis + homework
## =========================================================

set.seed(1234)

library(immunedeconv)
library(RColorBrewer)

data_dir <- "/gpfs/gpfs1/scratch/c9881023/shared/Data"

## =========================================================
## 1. Load Finotello et al. TPM data
## =========================================================
finotello_tpm <- readRDS(file.path(data_dir, "finotello_tpm.rds"))

dim(finotello_tpm)
head(finotello_tpm)

## =========================================================
## 2. MCP-counter analysis
## =========================================================
finotello_mcp <- deconvolute_mcp_counter(finotello_tpm)

dim(finotello_mcp)
head(finotello_mcp)
rownames(finotello_mcp)
apply(finotello_mcp, 2, sum)

## Neutrophils barplot
celltype <- "Neutrophils"
barplot(
  finotello_mcp[celltype, ],
  main = paste(celltype, "(MCP-counter)"),
  las = 2,
  ylab = "Cell-type score",
  col = "#cab2d6",
  border = FALSE
)

## CD8 T cells barplot
celltype <- "CD8 T cells"
barplot(
  finotello_mcp[celltype, ],
  main = paste(celltype, "(MCP-counter)"),
  las = 2,
  ylab = "Cell-type score",
  col = "#fdbf6f",
  border = FALSE
)

which.max(finotello_mcp[celltype, ])
colnames(finotello_mcp)[which.max(finotello_mcp[celltype, ])]

## =========================================================
## 3. quanTIseq analysis
## =========================================================
finotello_quanTIseq <- deconvolute_quantiseq(
  finotello_tpm,
  arrays = FALSE,
  tumor = FALSE,
  scale_mrna = TRUE
)

dim(finotello_quanTIseq)
head(finotello_quanTIseq)
rownames(finotello_quanTIseq)
apply(finotello_quanTIseq, 2, sum)

pal_12 <- brewer.pal(12, "Set3")

barplot(
  finotello_quanTIseq,
  col = pal_12,
  ylab = "Cell-type fractions",
  las = 2,
  main = "quanTIseq deconvolution results",
  border = FALSE
)

par(mar = c(5, 4, 4, 10))
barplot(
  finotello_quanTIseq,
  col = pal_12,
  ylab = "Cell-type fractions",
  las = 2,
  main = "quanTIseq deconvolution results",
  border = FALSE
)

legend(
  "right",
  legend = rownames(finotello_quanTIseq),
  fill = pal_12,
  cex = 0.7,
  bty = "n",
  border = FALSE,
  xpd = TRUE,
  inset = -0.45
)

## MCP-counter stacked plot (not meaningful)
par(mar = c(5, 4, 4, 10))
barplot(
  finotello_mcp,
  col = pal_12,
  ylab = "Cell-type scores",
  las = 2,
  main = "MCP-counter results - NOT MEANINGFUL",
  border = FALSE
)

legend(
  "right",
  legend = rownames(finotello_mcp),
  fill = pal_12,
  cex = 0.7,
  bty = "n",
  border = FALSE,
  xpd = TRUE,
  inset = -0.3
)

## =========================================================
## 4. Compare MCP-counter and quanTIseq
## =========================================================
cd8_mcp <- as.numeric(finotello_mcp["Cytotoxic lymphocytes", ])
cd8_qts <- as.numeric(finotello_quanTIseq["T.cells.CD8", ])

plot(
  cd8_mcp, cd8_qts,
  xlab = "MCP-counter (cell scores)",
  ylab = "quanTIseq (cell fractions)",
  pch = 19,
  cex = 1.5,
  main = "CD8+ T cells"
)

r_cd8 <- round(cor(cd8_mcp, cd8_qts), 2)
text(
  x = mean(range(cd8_mcp)),
  y = max(cd8_qts) * 0.9,
  labels = paste("r =", r_cd8)
)

## Neutrophils comparison
neut_mcp <- as.numeric(finotello_mcp["Neutrophils", ])
neut_qts <- as.numeric(finotello_quanTIseq["Neutrophils", ])

plot(
  neut_mcp, neut_qts,
  xlab = "MCP-counter (cell scores)",
  ylab = "quanTIseq (cell fractions)",
  pch = 19,
  cex = 1.5,
  main = "Neutrophils"
)

r_neut <- round(cor(neut_mcp, neut_qts), 2)
text(
  x = mean(range(neut_mcp)),
  y = max(neut_qts) * 0.9,
  labels = paste("r =", r_neut)
)

r_neut

## =========================================================
## 5. Load FACS ground truth
## =========================================================
finotello_truth <- readRDS(file.path(data_dir, "finotello_facs.rds"))

dim(finotello_truth)
head(finotello_truth)
rownames(finotello_truth)

ccells <- intersect(rownames(finotello_quanTIseq), rownames(finotello_truth))
csamples <- intersect(colnames(finotello_quanTIseq), colnames(finotello_truth))

finotello_quanTIseq_sel <- finotello_quanTIseq[ccells, csamples, drop = FALSE]
finotello_truth_sel <- finotello_truth[ccells, csamples, drop = FALSE]

colnames(finotello_quanTIseq_sel)[5]
colnames(finotello_truth_sel)[5]
rownames(finotello_quanTIseq_sel)[3]
rownames(finotello_truth_sel)[3]

## =========================================================
## 6. Per-cell-type comparison: quanTIseq vs ground truth
## =========================================================
old_par <- par(no.readonly = TRUE)
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))

for (ct in ccells) {
  x <- as.numeric(finotello_truth_sel[ct, ])
  y <- as.numeric(finotello_quanTIseq_sel[ct, ])
  
  r <- round(cor(x, y), 2)
  
  plot(
    x, y,
    xlab = "True cell fractions",
    ylab = "Estimated cell fractions",
    main = paste0(ct, "\n r = ", r),
    pch = 19,
    col = "#1f78b4",
    xlim = range(c(x, y)),
    ylim = range(c(x, y))
  )
  abline(0, 1, lty = 2)
}

par(old_par)

cor(
  as.numeric(finotello_truth_sel["T.cells.CD8", ]),
  as.numeric(finotello_quanTIseq_sel["T.cells.CD8", ])
)

cor(
  as.numeric(finotello_truth_sel["Tregs", ]),
  as.numeric(finotello_quanTIseq_sel["Tregs", ])
)

## =========================================================
## 7. Global comparison across all cell types
## =========================================================
finotello_truth_sel_v <- as.numeric(finotello_truth_sel)
finotello_quanTIseq_sel_v <- as.numeric(finotello_quanTIseq_sel)

labels <- rep(rownames(finotello_truth_sel), times = ncol(finotello_truth_sel))
colors <- pal_12[as.factor(labels)]

r_all <- round(cor(finotello_truth_sel_v, finotello_quanTIseq_sel_v), 2)

par(mar = c(5, 4, 4, 14))
plot(
  finotello_truth_sel_v,
  finotello_quanTIseq_sel_v,
  xlab = "True cell fractions",
  ylab = "Estimated cell fractions",
  main = "quanTIseq (all cell types)",
  xlim = c(0, 0.4),
  ylim = c(0, 0.4),
  pch = 19,
  col = colors
)

abline(0, 1)
text(0.05, 0.38, paste0("r = ", r_all))

legend(
  "right",
  legend = unique(labels),
  fill = unique(colors),
  cex = 0.7,
  y.intersp = 0.8,
  x.intersp = 0.7,
  bty = "n",
  border = FALSE,
  xpd = TRUE,
  inset = -0.45
)

## =========================================================
## 8. Homework dataset: Hoek et al.
## =========================================================
hoek_tpm <- readRDS(file.path(data_dir, "hoek_tpm.rds"))
hoek_truth <- readRDS(file.path(data_dir, "hoek_facs.rds"))

rownames(hoek_truth)

hoek_epic <- deconvolute_epic(
  hoek_tpm,
  tumor = FALSE,
  scale_mrna = TRUE
)

dim(hoek_epic)
head(hoek_epic)
rownames(hoek_epic)
apply(hoek_epic, 2, sum)

## Rename cell types to match ground truth
idx <- which(rownames(hoek_epic) == "Bcells")
if (length(idx) > 0) {
  rownames(hoek_epic)[idx] <- "B cells"
}

idx <- which(rownames(hoek_epic) == "NKcells")
if (length(idx) > 0) {
  rownames(hoek_epic)[idx] <- "NK cells"
}

## Add T cell row = CD4 + CD8
hoek_epic <- rbind(
  hoek_epic,
  "T cell" = hoek_epic["CD4_Tcells", ] + hoek_epic["CD8_Tcells", ]
)

rownames(hoek_epic)
rownames(hoek_truth)

pal_epic <- brewer.pal(8, "Set2")

par(mar = c(5, 4, 4, 10))
barplot(
  hoek_epic,
  col = pal_epic,
  ylab = "Cell-type fractions",
  las = 2,
  main = "EPIC deconvolution results",
  border = FALSE
)

legend(
  "right",
  legend = rownames(hoek_epic),
  fill = pal_epic,
  cex = 0.7,
  bty = "n",
  border = FALSE,
  xpd = TRUE,
  inset = -0.35
)

ccells_hoek <- intersect(rownames(hoek_epic), rownames(hoek_truth))
csamples_hoek <- intersect(colnames(hoek_epic), colnames(hoek_truth))

hoek_epic_sel <- hoek_epic[ccells_hoek, csamples_hoek, drop = FALSE]
hoek_truth_sel <- hoek_truth[ccells_hoek, csamples_hoek, drop = FALSE]

ccells_hoek
csamples_hoek

old_par <- par(no.readonly = TRUE)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

for (ct in ccells_hoek) {
  x <- as.numeric(hoek_truth_sel[ct, ])
  y <- as.numeric(hoek_epic_sel[ct, ])
  
  r <- round(cor(x, y), 2)
  
  plot(
    x, y,
    xlab = "True cell fractions",
    ylab = "Estimated cell fractions",
    main = paste0(ct, "\n r = ", r),
    pch = 19,
    col = "#33a02c",
    xlim = range(c(x, y)),
    ylim = range(c(x, y))
  )
  abline(0, 1, lty = 2)
}

par(old_par)

hoek_truth_sel_v <- as.numeric(as.matrix(hoek_truth_sel))
hoek_epic_sel_v <- as.numeric(as.matrix(hoek_epic_sel))

labels_hoek <- rep(rownames(hoek_truth_sel), times = ncol(hoek_truth_sel))
colors_hoek <- brewer.pal(max(3, length(unique(labels_hoek))), "Set2")[as.factor(labels_hoek)]

r_hoek_all <- round(cor(hoek_truth_sel_v, hoek_epic_sel_v), 2)

par(mar = c(5, 4, 4, 12))
plot(
  hoek_truth_sel_v,
  hoek_epic_sel_v,
  xlab = "True cell fractions",
  ylab = "Estimated cell fractions",
  main = "EPIC (all common cell types)",
  pch = 19,
  col = colors_hoek
)

abline(0, 1)
usr <- par("usr")
text(
  x = usr[1] + 0.1 * (usr[2] - usr[1]),
  y = usr[4] - 0.1 * (usr[4] - usr[3]),
  labels = paste0("r = ", r_hoek_all)
)

legend(
  "right",
  legend = unique(labels_hoek),
  fill = unique(colors_hoek),
  cex = 0.8,
  bty = "n",
  border = FALSE,
  xpd = TRUE,
  inset = -0.3
)

## =========================================================
## 9. Optional: numeric summary of per-cell-type correlations
## =========================================================
finotello_corrs <- sapply(ccells, function(ct) {
  cor(
    as.numeric(finotello_truth_sel[ct, ]),
    as.numeric(finotello_quanTIseq_sel[ct, ])
  )
})
finotello_corrs

hoek_corrs <- sapply(ccells_hoek, function(ct) {
  cor(
    as.numeric(hoek_truth_sel[ct, ]),
    as.numeric(hoek_epic_sel[ct, ])
  )
})
hoek_corrs
