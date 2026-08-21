# ═══════════════════════════════════════════════════════════════════════════
# BALL MAPPER — R IMPLEMENTATION
# ═══════════════════════════════════════════════════════════════════════════

# ── 1. PACKAGES ───────────────────────────────────────────────────────────
library(BallMapper)
library(igraph)
library(fields)

# ── 2. PATHS ──────────────────────────────────────────────────────────────
DATA_PATH <- "C:/Users/namra/gluten-free-spatial-inequality/data/processed/analytical_dataset_final.csv"
FIG_PATH  <- "C:/Users/namra/gluten-free-spatial-inequality/outputs/figures"
TAB_PATH  <- "C:/Users/namra/gluten-free-spatial-inequality/outputs/tables"

# ── 3. DATA ───────────────────────────────────────────────────────────────
df <- read.csv(DATA_PATH, stringsAsFactors = FALSE)
cat("Dataset loaded:", nrow(df), "rows,", ncol(df), "columns\n")

# ── 4. POINT CLOUD ────────────────────────────────────────────────────────
vars <- c("prescribing_rate_recent",
          "retail_accessibility_score",
          "imd_score_weighted")

X <- scale(df[, vars])

cat("\nPoint cloud: 106 areas x 3 axes\n")
cat("  Axis 1: prescribing rate\n")
cat("  Axis 2: retail accessibility score\n")
cat("  Axis 3: population-weighted IMD score\n")

# ── 5. EPSILON SENSITIVITY ────────────────────────────────────────────────
epsilons <- c(0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0)

sens <- data.frame(
  epsilon    = numeric(),
  nodes      = integer(),
  edges      = integer(),
  components = integer(),
  avg_size   = numeric()
)

for (eps in epsilons) {
  bm_e  <- BallMapper(
    points  = X,
    values  = data.frame(prescribing = df$prescribing_rate_recent),
    epsilon = eps
  )
  n_e   <- nrow(bm_e$vertices)
  ed_e  <- nrow(bm_e$edges)
  mem_e <- bm_e$points_covered_by_landmarks
  
  g_e <- if (ed_e > 0) {
    graph_from_data_frame(
      as.data.frame(bm_e$edges), directed = FALSE,
      vertices = data.frame(id = 1:n_e))
  } else {
    make_empty_graph(n = n_e, directed = FALSE)
  }
  
  sens <- rbind(sens, data.frame(
    epsilon    = eps,
    nodes      = n_e,
    edges      = ed_e,
    components = components(g_e)$no,
    avg_size   = mean(sapply(mem_e, length))
  ))
}

cat("\nEpsilon sensitivity:\n")
cat(sprintf("%10s %8s %8s %12s %10s\n",
            "Epsilon", "Nodes", "Edges", "Components", "Avg size"))
cat(strrep("-", 52), "\n")
for (i in 1:nrow(sens)) {
  mk <- if (sens$epsilon[i] == 1.5) "  << SELECTED" else ""
  cat(sprintf("%10.2f %8d %8d %12d %10.1f%s\n",
              sens$epsilon[i], sens$nodes[i], sens$edges[i],
              sens$components[i], sens$avg_size[i], mk))
}

# ── 5a. EPSILON SENSITIVITY PLOT ──────────────────────────────────────────
png(file.path(FIG_PATH, "BM_R_01_epsilon_sensitivity.png"),
    width = 1600, height = 700, res = 150)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2), bg = "white")

# Panel 1 — nodes and edges
plot(sens$epsilon, sens$nodes,
     type = "b", col = "#2a78d6", pch = 16, lwd = 2,
     ylim = c(0, max(sens$nodes, sens$edges) + 5),
     xlab = expression(paste("Epsilon (", epsilon, ")")),
     ylab = "Count",
     main = "Graph size vs epsilon")
lines(sens$epsilon, sens$edges,
      type = "b", col = "#1baf7a", pch = 15, lwd = 2)
abline(v = 1.5, col = "#e34948", lwd = 2, lty = 2)
legend("topright",
       legend = c("Nodes", "Edges",
                  expression(paste("Selected ", epsilon, " = 1.5"))),
       col = c("#2a78d6", "#1baf7a", "#e34948"),
       lwd = 2, pch = c(16, 15, NA),
       lty = c(1, 1, 2), cex = 0.85)

# Panel 2 — components
plot(sens$epsilon, sens$components,
     type = "b", col = "#4a3aa7", pch = 16, lwd = 2,
     xlab = expression(paste("Epsilon (", epsilon, ")")),
     ylab = "Connected components",
     main = "Graph connectivity vs epsilon")
abline(h = 1, col = "#73726c", lwd = 1, lty = 3)
abline(v = 1.5, col = "#e34948", lwd = 2, lty = 2)
text(x = 1.6,
     y = max(sens$components) * 0.75,
     labels = expression(paste(epsilon, " = 1.5")),
     col = "#73726c", cex = 0.75, pos = 4)
legend("topright",
       legend = c("Components",
                  expression(paste("Selected ", epsilon, " = 1.5"))),
       col = c("#4a3aa7", "#e34948"),
       lwd = 2, lty = c(1, 2), pch = c(16, NA), cex = 0.85)

dev.off()
cat("Saved BM_R_01_epsilon_sensitivity.png\n")

# ── 5b. BALL MAPPER GRAPHS AT FIVE EPSILON VALUES ─────────────────────────
epsilons_plot <- c(1.0, 1.25, 1.5, 1.75, 2.0)
titles_plot   <- c(
  expression(epsilon ~ "= 1.00"),
  expression(epsilon ~ "= 1.25"),
  expression(epsilon ~ "= 1.50  (selected)"),
  expression(epsilon ~ "= 1.75"),
  expression(epsilon ~ "= 2.00")
)

bm_list   <- list()
all_means <- c()

for (eps in epsilons_plot) {
  bm_e <- BallMapper(
    points  = X,
    values  = data.frame(prescribing = df$prescribing_rate_recent),
    epsilon = eps
  )
  bm_list[[as.character(eps)]] <- bm_e
  all_means <- c(all_means,
                 sapply(bm_e$points_covered_by_landmarks, function(m)
                   mean(df$prescribing_rate_recent[m])))
}

vmin_sh <- min(all_means)
vmax_sh <- max(all_means)

rdylbu_pal <- colorRampPalette(c(
  "#d73027","#f46d43","#fdae61","#fee090",
  "#e0f3f8","#abd9e9","#74add1","#4575b4"))(100)

colour_from_shared <- function(means) {
  norm <- pmax(0, pmin(1, (means - vmin_sh) / (vmax_sh - vmin_sh)))
  rdylbu_pal[ceiling(norm * 99 + 1)]
}

png(file.path(FIG_PATH, "BM_R_02_epsilon_graphs.png"),
    width = 7000, height = 4200, res = 150)

layout(
  matrix(c(1, 2, 3,
           4, 5, 0,
           6, 6, 6),
         nrow = 3, byrow = TRUE),
  widths  = c(1, 1, 1),
  heights = c(1, 1, 0.18)
)

par(bg = "white")

for (k in seq_along(epsilons_plot)) {
  eps   <- epsilons_plot[k]
  bm_e  <- bm_list[[as.character(eps)]]
  mem_e <- bm_e$points_covered_by_landmarks
  sz_e  <- sapply(mem_e, length)
  n_e   <- nrow(bm_e$vertices)
  ed_e  <- nrow(bm_e$edges)
  
  g_e <- if (ed_e > 0) {
    graph_from_data_frame(
      as.data.frame(bm_e$edges), directed = FALSE,
      vertices = data.frame(id = 1:n_e))
  } else {
    make_empty_graph(n = n_e, directed = FALSE)
  }
  
  means_e <- sapply(mem_e, function(m) mean(df$prescribing_rate_recent[m]))
  cols_e  <- colour_from_shared(means_e)
  
  set.seed(42)
  lay_e <- layout_with_fr(g_e, niter = 500)
  
  par(mar = c(0, 0, 3, 0))
  
  plot(g_e,
       layout             = lay_e,
       rescale            = TRUE,
       asp                = 0,
       vertex.size        = sqrt(sz_e) * 5,
       vertex.color       = cols_e,
       vertex.frame.color = "white",
       vertex.frame.width = 1,
       vertex.label       = as.character(0:(n_e - 1)),
       vertex.label.cex   = 2.5,            
       vertex.label.color = "black",
       vertex.label.font  = 2,
       edge.color         = "#555555",
       edge.width         = 2.5,
       main               = "")       
  
  # Add title separately — this respects cex
  title(main     = titles_plot[[k]],
        cex.main = 4.0,
        font.main = 2,
        line      = 1)
}

# Full-width colourbar using smallplot to force exact coordinates
par(mar = c(0, 0, 0, 0))
plot.new()
image.plot(legend.only = TRUE,
           zlim        = c(vmin_sh, vmax_sh),
           col         = rdylbu_pal,
           horizontal  = TRUE,
           smallplot   = c(0.02, 0.98, 0.4, 0.7),
           axis.args   = list(
             cex.axis = 1.8,
             mgp      = c(3, 1.0, 0)),
           legend.args = list(
             text = "Mean NHS GF prescribing rate (£ per 1,000 registered patients)",
             side = 3,
             line = 0.8,
             cex  = 2.2,
             font = 2))

dev.off()
cat("Saved BM_R_02_epsilon_graphs.png\n")

# ── 6. BUILD BALL MAPPER AT ε = 1.5 ───────────────────────────────────────
EPSILON <- 1.5

bm <- BallMapper(
  points  = X,
  values  = data.frame(prescribing = df$prescribing_rate_recent),
  epsilon = EPSILON
)

members   <- bm$points_covered_by_landmarks
node_size <- sapply(members, length)
n_nodes   <- nrow(bm$vertices)
n_edges   <- nrow(bm$edges)
total_mem <- sum(node_size)

cat(sprintf("\nBall Mapper at ε = %.1f\n", EPSILON))
cat("Nodes:", n_nodes, "\n")
cat("Edges:", n_edges, "\n")
cat("Total membership:", total_mem, "\n")

# Build graph
g <- graph_from_data_frame(
  d        = as.data.frame(bm$edges),
  directed = FALSE,
  vertices = data.frame(id = 1:n_nodes)
)

# Components
comp      <- components(g)
main_comp <- as.integer(names(which.max(table(comp$membership))))
iso_nodes <- which(comp$membership != main_comp)

cat("\nConnected components:", comp$no, "\n")
for (ci in 1:comp$no) {
  comp_nodes <- which(comp$membership == ci)
  cat(sprintf("  Component %d: %d node(s), %d areas\n",
              ci, length(comp_nodes),
              sum(node_size[comp_nodes])))
  if (length(comp_nodes) <= 3) {
    for (nid in comp_nodes) {
      areas <- gsub("NHS ", "", df$org_name[members[[nid]]])
      cat(sprintf("    Node %d: %s\n", nid-1, paste(areas, collapse = ", ")))
    }
  }
}

# Fixed layout — seed ensures reproducibility
set.seed(42)
layout_bm <- layout_with_fr(g, niter = 1000)

# Frame colours — red border for isolated nodes
frame_cols <- rep("white", n_nodes)
frame_widths <- rep(1, n_nodes)

# ── 7. NODE SUMMARY ───────────────────────────────────────────────────────
cat("\nNode summary at ε = 1.5:\n")
cat(sprintf("%-5s %-5s %-14s %-10s %-8s %-8s %-8s\n",
            "Node", "n", "Prescribing", "Retail",
            "IMD", "Mths0", "Slope"))
cat(strrep("-", 65), "\n")

for (i in 1:n_nodes) {
  m <- members[[i]]
  cat(sprintf("  %-3d  %3d   £%6.1f       %.3f    %.1f    %.1f   %+.4f\n",
              i-1, length(m),
              mean(df$prescribing_rate_recent[m]),
              mean(df$retail_accessibility_score[m]),
              mean(df$imd_score_weighted[m]),
              mean(df$months_zero[m]),
              mean(df$trend_slope[m])))
}

# ── 8. HELPER FUNCTION — PURPOSE-SPECIFIC PALETTES ────────────────────────
# Each variable gets a semantically appropriate diverging palette

make_pal <- function(colours, n = 100) colorRampPalette(colours)(n)

palettes <- list(
  prescribing = make_pal(c(
    "#d73027","#f46d43","#fdae61","#fee090",
    "#e0f3f8","#abd9e9","#74add1","#4575b4")),  # red=low, blue=high
  retail      = make_pal(c(
    "#d73027","#f46d43","#fdae61","#fee08b",
    "#d9ef8b","#91cf60","#1a9850")),             # red=low, green=high
  imd         = make_pal(c(
    "#4575b4","#74add1","#abd9e9","#e0f3f8",
    "#fee090","#fdae61","#f46d43","#d73027"))   # blue=low/good, red=high/bad
)

node_cols_pal <- function(values, members, pal) {
  means <- sapply(members, function(m) mean(values[m]))
  norm  <- pmax(0, pmin(1, (means - min(means)) / (max(means) - min(means))))
  list(
    cols = pal[ceiling(norm * 99 + 1)],
    min  = min(means),
    max  = max(means),
    pal  = pal
  )
}

# ── 9. TWO-PANEL COLOURED GRAPH (retail and IMD) ──────────────────────────
panels <- list(
  list(
    values = df$retail_accessibility_score,
    pal    = palettes$retail,
    label  = "GF Retail Accessibility Score\n(red = low, green = high)"
  ),
  list(
    values = df$imd_score_weighted,
    pal    = palettes$imd,
    label  = "Area Deprivation (IMD 2025)\n(blue = low, red = high)"
  )
)

png(file.path(FIG_PATH, "BM_R_03_two_panels.png"),
    width = 3200, height = 1400, res = 150)

# Simple 1x2 grid — one graph per cell, colourbar inside right margin
par(mfrow = c(1, 2), bg = "white")

for (k in seq_along(panels)) {
  panel <- panels[[k]]
  res   <- node_cols_pal(panel$values, members, panel$pal)
  
  # Wide right margin to fit colourbar inside the panel
  par(mar = c(2, 2, 4, 6))
  
  plot(g,
       layout             = layout_bm,
       vertex.size        = sqrt(node_size) * 7,
       vertex.color       = res$cols,
       vertex.frame.color = frame_cols,
       vertex.frame.width = frame_widths,
       vertex.label       = as.character(0:(n_nodes - 1)),
       vertex.label.cex   = 1.1,
       vertex.label.color = "black",
       vertex.label.font  = 2,
       edge.color         = "#73726c",
       edge.width         = 1.5,
       main               = panel$label,
       cex.main           = 1.0)
  
  # Colourbar placed in the right margin of the current panel
  image.plot(legend.only = TRUE,
             zlim        = c(res$min, res$max),
             col         = res$pal,
             smallplot   = c(0.88, 0.92, 0.15, 0.85),
             axis.args   = list(cex.axis = 0.85, mgp = c(3, 0.5, 0)),
             legend.args = list(text = "", side = 3))
}

dev.off()
cat("Saved BM_R_03_two_panels.png\n")

# ── 10. ANNOTATED SINGLE GRAPH (PRESCRIBING) ──────────────────────────────
res_presc <- node_cols_pal(
  df$prescribing_rate_recent, members, palettes$prescribing)

png(file.path(FIG_PATH, "BM_R_04_annotated.png"),
    width = 1800, height = 1200, res = 150)

layout(matrix(c(1, 2), nrow = 1), widths = c(1, 0.15))

par(mar = c(2, 2, 4, 1), bg = "white")
plot(g,
     layout             = layout_bm,
     vertex.size        = sqrt(node_size) * 9,
     vertex.color       = res_presc$cols,
     vertex.frame.color = frame_cols,
     vertex.frame.width = frame_widths,
     vertex.label       = as.character(0:(n_nodes - 1)),
     vertex.label.cex   = 1.3,
     vertex.label.color = "black",
     vertex.label.font  = 2,
     edge.color         = "#73726c",
     edge.width         = 2,
     main = paste0(
       "Ball Mapper (\u03b5 = ", EPSILON,
       ") \u2014 106 English Sub-ICB Areas\n",
       "Coloured by NHS GF prescribing rate",
       " | Red = low provision, Blue = high provision"))

par(mar = c(6, 0, 6, 4))
image.plot(legend.only  = TRUE,
           zlim         = c(res_presc$min, res_presc$max),
           col          = res_presc$pal,
           legend.width = 2,
           axis.args    = list(cex.axis = 1.0),
           legend.args  = list(
             text = "Prescribing rate (£)",
             side = 3, line = 1, cex = 0.9, font = 2))

dev.off()
cat("Saved BM_R_04_annotated.png\n")

# ── 11. NODE MEMBERSHIP CSV ────────────────────────────────────────────────
rows <- list()
for (i in 1:n_nodes) {
  m <- members[[i]]
  for (idx in m) {
    rows[[length(rows) + 1]] <- data.frame(
      node_id                    = i - 1,
      component                  = comp$membership[i],
      org_id                     = df$org_id[idx],
      org_name                   = df$org_name[idx],
      prescribing_rate_recent    = df$prescribing_rate_recent[idx],
      retail_accessibility_score = df$retail_accessibility_score[idx],
      imd_score_weighted         = df$imd_score_weighted[idx],
      months_zero                = df$months_zero[idx],
      trend_slope                = df$trend_slope[idx]
    )
  }
}

membership_r <- do.call(rbind, rows)
write.csv(
  membership_r,
  file.path(TAB_PATH, "BM_R_node_membership.csv"),
  row.names = FALSE
)
cat("Saved BM_R_node_membership.csv\n")

# ── 12. NODE SUMMARY CSV ──────────────────────────────────────────────────
node_summary <- data.frame(
  node_id    = 0:(n_nodes - 1),
  component  = comp$membership,
  n_areas    = node_size,
  mean_prescribing       = sapply(members, function(m)
    round(mean(df$prescribing_rate_recent[m]), 3)),
  mean_retail            = sapply(members, function(m)
    round(mean(df$retail_accessibility_score[m]), 3)),
  mean_imd               = sapply(members, function(m)
    round(mean(df$imd_score_weighted[m]), 3)),
  mean_months_zero       = sapply(members, function(m)
    round(mean(df$months_zero[m]), 3)),
  mean_trend_slope       = sapply(members, function(m)
    round(mean(df$trend_slope[m]), 4))
)

write.csv(node_summary,
          file.path(TAB_PATH, "BM_R_node_summary.csv"),
          row.names = FALSE)
cat("Saved BM_R_node_summary.csv\n")

# ── 13. ARI BETWEEN EPSILON VALUES ────────────────────────────────────────
# Requires mclust package for ARI
if (!requireNamespace("mclust", quietly = TRUE)) install.packages("mclust")
library(mclust)

# Primary node assignment function
primary_partition <- function(bm_obj, n_points) {
  labels <- rep(-1, n_points)
  for (i in seq_along(bm_obj$points_covered_by_landmarks)) {
    for (idx in bm_obj$points_covered_by_landmarks[[i]]) {
      if (labels[idx] == -1) labels[idx] <- i
    }
  }
  labels
}

eps_ari <- c(1.25, 1.5, 1.75, 2.0)
bm_ari_list <- list()

for (eps in eps_ari) {
  bm_ari_list[[as.character(eps)]] <- BallMapper(
    points  = X,
    values  = data.frame(prescribing = df$prescribing_rate_recent),
    epsilon = eps
  )
}

part_125 <- primary_partition(bm_ari_list[["1.25"]], nrow(df))
part_150 <- primary_partition(bm_ari_list[["1.5"]],  nrow(df))
part_175 <- primary_partition(bm_ari_list[["1.75"]], nrow(df))
part_200 <- primary_partition(bm_ari_list[["2"]],    nrow(df))

ari_125_150 <- adjustedRandIndex(part_125, part_150)
ari_150_175 <- adjustedRandIndex(part_150, part_175)
ari_175_200 <- adjustedRandIndex(part_175, part_200)

cat("\nAdjusted Rand Index — epsilon robustness:\n")
cat(sprintf("  ARI(ε=1.25 vs ε=1.50): %.3f\n", ari_125_150))
cat(sprintf("  ARI(ε=1.50 vs ε=1.75): %.3f\n", ari_150_175))
cat(sprintf("  ARI(ε=1.75 vs ε=2.00): %.3f\n", ari_175_200))
cat("  (Python values: 0.533, 0.402, 0.624 — confirm match)\n")

# ── 14. BINARY RETAIL ROBUSTNESS ──────────────────────────────────────────
# Recompute retail score using binary scheme (0=budget, 1=non-budget)
df$retail_binary <- (df$regular_stores + df$quality_stores) / df$total_stores

cat(sprintf("\nBinary retail correlation with original: %.3f\n",
            cor(df$retail_accessibility_score, df$retail_binary)))

# Build binary point cloud
X_binary <- scale(cbind(
  df$prescribing_rate_recent,
  df$retail_binary,
  df$imd_score_weighted
))

bm_binary <- BallMapper(
  points  = X_binary,
  values  = data.frame(prescribing = df$prescribing_rate_recent),
  epsilon = EPSILON
)

members_bin   <- bm_binary$points_covered_by_landmarks
node_size_bin <- sapply(members_bin, length)
n_nodes_bin   <- nrow(bm_binary$vertices)
n_edges_bin   <- nrow(bm_binary$edges)

cat(sprintf("Binary retail BM at ε=1.5: %d nodes, %d edges\n",
            n_nodes_bin, n_edges_bin))

# ARI between original and binary partitions
part_orig   <- primary_partition(bm, nrow(df))
part_binary <- primary_partition(bm_binary, nrow(df))
ari_binary  <- adjustedRandIndex(part_orig, part_binary)
cat(sprintf("ARI (original vs binary retail): %.3f\n", ari_binary))
cat("  (Python value: 0.605 — confirm match)\n")

# Check if N3 remains isolated in binary scheme
g_bin <- graph_from_data_frame(
  as.data.frame(bm_binary$edges), directed = FALSE,
  vertices = data.frame(id = 1:n_nodes_bin)
)
comp_bin <- components(g_bin)
cat(sprintf("Binary scheme components: %d\n", comp_bin$no))
for (ci in 1:comp_bin$no) {
  comp_nodes_bin <- which(comp_bin$membership == ci)
  if (length(comp_nodes_bin) <= 3) {
    for (nid in comp_nodes_bin) {
      areas <- gsub("NHS ", "", df$org_name[members_bin[[nid]]])
      cat(sprintf("  Isolated node %d: %s\n",
                  nid - 1, paste(areas, collapse = ", ")))
    }
  }
}

# Side-by-side comparison figure
set.seed(42)
layout_bin <- layout_with_fr(g_bin, niter = 1000)

# Shared colour scale across both graphs
all_means_rob <- c(
  sapply(members,     function(m) mean(df$prescribing_rate_recent[m])),
  sapply(members_bin, function(m) mean(df$prescribing_rate_recent[m]))
)
vmin_rob <- min(all_means_rob)
vmax_rob <- max(all_means_rob)

colour_rob <- function(members_list) {
  means <- sapply(members_list, function(m) mean(df$prescribing_rate_recent[m]))
  norm  <- pmax(0, pmin(1, (means - vmin_rob) / (vmax_rob - vmin_rob)))
  rdylbu_pal[ceiling(norm * 99 + 1)]
}

cols_orig <- colour_rob(members)
cols_bin  <- colour_rob(members_bin)

png(file.path(FIG_PATH, "BM_R_05_binary_retail_robustness.png"),
    width = 3600, height = 1600, res = 150)

layout(matrix(c(1, 2, 3), nrow = 1),
       widths = c(1, 1, 0.2))
par(bg = "white")

# Original
par(mar = c(1, 1, 4, 1))
plot(g,
     layout             = layout_bm,
     vertex.size        = sqrt(node_size) * 7,
     vertex.color       = cols_orig,
     vertex.frame.color = "white",
     vertex.label       = as.character(0:(n_nodes - 1)),
     vertex.label.cex   = 1.1,
     vertex.label.color = "black",
     vertex.label.font  = 2,
     edge.color         = "#73726c",
     edge.width         = 1.5,
     main               = "Original retail weighting\n(budget=0, regular=1, quality=2)")

# Binary
par(mar = c(1, 1, 4, 1))
plot(g_bin,
     layout             = layout_bin,
     vertex.size        = sqrt(node_size_bin) * 7,
     vertex.color       = cols_bin,
     vertex.frame.color = "white",
     vertex.label       = as.character(0:(n_nodes_bin - 1)),
     vertex.label.cex   = 1.1,
     vertex.label.color = "black",
     vertex.label.font  = 2,
     edge.color         = "#73726c",
     edge.width         = 1.5,
     main               = paste0("Binary retail weighting\n",
                                 "(budget=0, non-budget=1) | ARI=",
                                 round(ari_binary, 3)))

# Colourbar
par(mar = c(4, 0, 4, 3))
image.plot(legend.only   = TRUE,
           zlim          = c(vmin_rob, vmax_rob),
           col           = rdylbu_pal,
           legend.width  = 2,
           legend.shrink = 0.8,
           axis.args     = list(cex.axis = 0.9),
           legend.args   = list(
             text = "Prescribing rate (£)",
             side = 3, line = 1, cex = 0.9, font = 2))

dev.off()
cat("Saved BM_R_05_binary_retail_robustness.png\n")

# Save robustness summary
robustness_summary <- data.frame(
  comparison              = c("ε=1.25 vs ε=1.50",
                              "ε=1.50 vs ε=1.75",
                              "ε=1.75 vs ε=2.00",
                              "Original vs binary retail"),
  ARI                     = round(c(ari_125_150,
                                    ari_150_175,
                                    ari_175_200,
                                    ari_binary), 3),
  interpretation          = c("Moderate agreement",
                              "Moderate agreement",
                              "Strong agreement",
                              "Moderate-strong agreement")
)
write.csv(robustness_summary,
          file.path(TAB_PATH, "BM_R_robustness_summary.csv"),
          row.names = FALSE)
cat("Saved BM_R_robustness_summary.csv\n")
