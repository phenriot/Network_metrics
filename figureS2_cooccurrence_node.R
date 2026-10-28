# =============================================================================
# Figure S2 — Co-occurrence of node-level metrics across reviewed articles
# Cell [i, j] = number of references reporting both metric i and j
# Metrics shown individually (no aggregation), grouped by subgroup blocks.
# Output: FigureS2_cooccurrence_node.pdf
#
# Requires: ComplexHeatmap (Bioconductor), circlize
#   if (!require("BiocManager")) install.packages("BiocManager")
#   BiocManager::install("ComplexHeatmap")
#   install.packages("circlize")
# =============================================================================

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

source("load_metrics.R")

INPUT  <- "extraction_table.csv"
DICT   <- "metric_dictionary.csv"
LEVELS <- "Node-level"
CELL   <- 6.0   # mm per cell

# ── Load individual metrics (no aggregation) ─────────────────────────────────
d    <- load_metrics(INPUT, DICT, LEVELS, drop_unused = TRUE, drop_empty = TRUE)
mat  <- d$mat;  meta <- d$meta

# ── Split by supergroup + subgroup ───────────────────────────────────────────
sub_raw   <- meta$subgroup
split_key <- paste(meta$supergroup, sub_raw, sep = " | ")
split_f   <- factor(split_key, levels = unique(split_key))
slices    <- levels(split_f)

first       <- match(slices, split_key)
slice_super <- meta$supergroup[first]
super_lv    <- unique(meta$supergroup)

COL_MAX   <- "#2A5D78"
COL_MIN   <- "#CFE0EA"
COL_ZERO  <- "#F5F5F5"
COL_DIAG  <- "#D9D9D9"
COL_SUPER <- c("#C9DAE3", "#E7D8C3")
COL_SUB   <- "#F4F4F4"
fill_map  <- setNames(rep_len(COL_SUPER, length(super_lv)), super_lv)

n_col_slice <- as.integer(table(split_f)[slices])

# ── Co-occurrence matrix ──────────────────────────────────────────────────────
n_ref <- colSums(mat)
M     <- crossprod(mat);  storage.mode(M) <- "double"
diag(M) <- NA

pos     <- M[!is.na(M) & M > 0]
lo      <- if (length(pos)) min(pos) else 1
hi      <- if (length(pos)) max(pos) else 1
col_fun <- if (hi > lo) colorRamp2(c(0, lo, hi), c(COL_ZERO, COL_MIN, COL_MAX)) else
                        colorRamp2(c(0, lo),      c(COL_ZERO, COL_MIN))
thr     <- 0.55 * hi

# ── Measure text extents to size the annotation bands ────────────────────────
pdf(NULL)

slice_lab <- vapply(seq_along(slices), function(k)
  fit_rot_label(sub_raw[first[k]], n_col_slice[k] * CELL - 1), character(1))

sub_lines <- unlist(strsplit(slice_lab, "\n", fixed = TRUE))
sub_h <- max(vapply(sub_lines, function(s) convertWidth(
  grobWidth(textGrob(s, gp = gpar(fontsize = 7))), "mm", valueOnly = TRUE),
  numeric(1))) + 3

need_mm <- vapply(super_lv, function(g) convertWidth(
  grobWidth(textGrob(g, gp = gpar(fontsize = 10, fontface = "bold"))),
  "mm", valueOnly = TRUE), numeric(1))

lab_slice <- vapply(super_lv, function(g) {
  s     <- which(slice_super == g)
  avail <- n_col_slice[s] * CELL - 2
  cand  <- if (any(avail >= 0.8 * need_mm[[g]])) s[avail >= 0.8 * need_mm[[g]]] else
    s[which.max(n_col_slice[s])]
  cand[which.min(abs(cand - mean(range(s))))]
}, integer(1))

fits     <- lapply(seq_along(super_lv), function(i)
  fit_label(super_lv[i], n_col_slice[lab_slice[i]] * CELL - 2))
lab_text <- setNames(vapply(fits, `[[`, character(1), "text"), super_lv)
lab_size <- setNames(vapply(fits, `[[`, numeric(1), "size"), super_lv)
super_h  <- max(6, max(vapply(fits, `[[`, numeric(1), "lines")) * 3.1 + 2)

invisible(dev.off())

# ── Annotations ──────────────────────────────────────────────────────────────
top <- HeatmapAnnotation(
  super = anno_block(
    panel_fun = function(index, nm) {
      si <- match(nm, slices); g <- slice_super[si]
      grid.rect(gp = gpar(fill = fill_map[[g]], col = "grey35"))
      if (si == lab_slice[[g]])
        grid.text(lab_text[[g]], gp = gpar(fontsize = lab_size[[g]], fontface = "bold"))
    }, height = unit(super_h, "mm")),
  sub = anno_block(
    labels = slice_lab, labels_rot = 90,
    gp = gpar(fill = COL_SUB, col = "grey60"),
    labels_gp = gpar(fontsize = 7), height = unit(sub_h, "mm")),
  `n ref.` = anno_barplot(
    n_ref, gp = gpar(fill = COL_MAX, col = NA),
    height = unit(11, "mm"),
    axis_param = list(gp = gpar(fontsize = 6))),
  annotation_name_gp   = gpar(fontsize = 7),
  annotation_name_side = "left",
  show_annotation_name = c(super = FALSE, sub = FALSE, `n ref.` = TRUE),
  gap = unit(0.5, "mm"))

left_ann <- rowAnnotation(
  super = anno_block(
    panel_fun = function(index, nm) {
      g <- slice_super[match(nm, slices)]
      grid.rect(gp = gpar(fill = fill_map[[g]], col = "grey35"))
    }, width = unit(3, "mm")),
  show_annotation_name = FALSE)

right_ann <- rowAnnotation(
  `n ref.` = anno_barplot(
    n_ref, gp = gpar(fill = COL_MAX, col = NA),
    width = unit(16, "mm"),
    axis_param = list(gp = gpar(fontsize = 6))),
  annotation_name_gp  = gpar(fontsize = 7),
  annotation_name_rot = 0)

# ── Heatmap ───────────────────────────────────────────────────────────────────
ht <- Heatmap(
  M, name = "n references",
  col     = col_fun, na_col = COL_DIAG,
  rect_gp = gpar(col = "white", lwd = 0.6),
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_split    = split_f, column_split = split_f,
  row_title    = NULL,    column_title  = NULL,
  row_gap      = unit(1.2, "mm"), column_gap = unit(1.2, "mm"),
  border       = TRUE, border_gp = gpar(col = "grey40", lwd = 0.8),
  row_names_side    = "left",   row_names_gp    = gpar(fontsize = 7),
  column_names_side = "bottom", column_names_gp = gpar(fontsize = 7),
  column_names_rot  = 90, column_names_max_height = unit(140, "mm"),
  width  = unit(CELL * ncol(M), "mm"),
  height = unit(CELL * nrow(M), "mm"),
  top_annotation   = top,
  left_annotation  = left_ann,
  right_annotation = right_ann,
  cell_fun = function(j, i, x, y, w, h, fill) {
    v <- M[i, j]
    if (!is.na(v) && v >= 1)
      grid.text(sprintf("%d", as.integer(v)), x, y,
                gp = gpar(fontsize = 5,
                          col = if (v >= thr) "white" else "grey25"))
  },
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 7), title_gp = gpar(fontsize = 8),
    legend_width = unit(35, "mm"), direction = "horizontal"))

# ── Export ────────────────────────────────────────────────────────────────────
write.csv(
  data.frame(metric = rownames(M), M, check.names = FALSE),
  "FigureS2_cooccurrence_node_matrix.csv", row.names = FALSE)

pdf(NULL)
drawn <- draw(ht, heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 2), "mm"))
w <- convertWidth( ComplexHeatmap:::width(drawn),  "in", valueOnly = TRUE)
h <- convertHeight(ComplexHeatmap:::height(drawn), "in", valueOnly = TRUE)
invisible(dev.off())

pdf("FigureS2_cooccurrence_node.pdf", width = w, height = h)
draw(ht, heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 2), "mm"))
invisible(dev.off())

message("OK -> FigureS2_cooccurrence_node.pdf  (",
        round(w * 2.54, 1), " x ", round(h * 2.54, 1), " cm)  |  ",
        nrow(M), " x ", ncol(M), " metrics")
