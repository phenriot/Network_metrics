# =============================================================================
# Figure 3 — UpSet panels: frequency of use and combinations of metric groups
# Combinations of Network-level (A) and Node-level (B) metrics.
# mode = combination, temporal metrics collapsed to one row
# Output: Figure3_upset_panels.pdf
#         Figure3_upset_network_combinations.csv
#         Figure3_upset_node_combinations.csv
#
# Requires: ComplexHeatmap (Bioconductor)
#   if (!require("BiocManager")) install.packages("BiocManager")
#   BiocManager::install("ComplexHeatmap")
# =============================================================================

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(grid)
})

metrics_dir <- "."
out_dir     <- getwd()   # output files go to the current working directory

# upset_core.R sources load_metrics.R with a relative path, so we temporarily
# set the working directory to metrics_dir before sourcing it.
old_wd <- setwd(metrics_dir)
source("upset_core.R")   # defines make_upset(), draw_upset(), upset_defaults()
setwd(old_wd)

# ── Build the two panels ──────────────────────────────────────────────────────

base_opt          <- upset_defaults()
base_opt$mode     <- "combination"   # one column per exact combination profile
base_opt$collapse <- "temporal"      # collapse temporal metrics into one row
base_opt$input    <- file.path(metrics_dir, "extraction_table.csv")

old_wd2 <- setwd(metrics_dir)        # load_metrics() is called with relative DICT

us <- list(
  A = make_upset(
    modifyList(base_opt, list(levels = "Network-level")),
    out_csv = file.path(out_dir, "Figure3_upset_network_combinations.csv")),
  B = make_upset(
    modifyList(base_opt, list(levels = "Node-level")),
    out_csv = file.path(out_dir, "Figure3_upset_node_combinations.csv"))
)

setwd(old_wd2)

# ── Page layout: two panels stacked, matrices left-aligned ───────────────────

w_i <- vapply(us, function(u) u$size[1], numeric(1))
h_i <- vapply(us, function(u) u$size[2], numeric(1))

# measure the row-name width of each panel to align matrices
pdf(NULL)
rn_i <- vapply(us, function(u) max(vapply(u$meta$label, function(s)
  convertWidth(grobWidth(textGrob(s, gp = gpar(fontsize = u$fs))),
               "in", valueOnly = TRUE), numeric(1))), numeric(1))
invisible(dev.off())

off_i <- max(rn_i) - rn_i      # horizontal offset so matrices line up
W     <- max(w_i + off_i)
H     <- sum(h_i)
x_i   <- off_i
y_i   <- 1 - (cumsum(h_i) - h_i) / H   # top of each panel (npc)

# ── Export ────────────────────────────────────────────────────────────────────

pdf(file.path(out_dir, "Figure3_upset_panels.pdf"), width = W, height = H)
grid.newpage()
for (i in seq_along(us)) {
  pushViewport(viewport(
    x      = unit(x_i[i], "in"),
    y      = unit(y_i[i], "npc"),
    width  = unit(w_i[i], "in"),
    height = unit(h_i[i], "in"),
    just   = c("left", "top")))
  draw_upset(us[[i]], newpage = FALSE)
  popViewport()
}
# panel letters placed after all panels are drawn (decorate_* moves the viewport)
upViewport(0)
grid.text(
  LETTERS[seq_along(us)],
  x    = unit(0, "npc") + unit(2, "mm"),
  y    = unit(y_i, "npc") - unit(2.5, "mm"),
  just = c("left", "top"),
  gp   = gpar(fontface = "bold", fontsize = 13))
invisible(dev.off())

message("OK -> Figure3_upset_panels.pdf  (",
        round(W * 2.54, 1), " x ", round(H * 2.54, 1), " cm)  |  ",
        paste(sprintf("%s=%d col.", LETTERS[seq_along(us)],
                      vapply(us, `[[`, numeric(1), "n_col")),
              collapse = ", "))
