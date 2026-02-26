library(ggplot2)
library(ggtree)
library(treeio)
library(viridis)

# Plot settings
tip_label_padding  <- 0.01
phasing_panel_prop <- 0.5
root_edge_split    <- 0.5
allele_mode        <- "vertical"
plot_filename      <- "example/example_vertical.png"

# Define the color palette for posterior probabilities
# (Note: in this code, this palette is for both the clade posteriors
# and phasing posteriors)
posterior_palette <- rev(viridis::turbo(256, begin = 0.1, end = 0.9))

# Example inputs
map_tree_file       <- "example/22B02_r1/22B02_r1_trees_relabeled_MAP.tre"
marginal_file       <- "example/22B02_r1/22B02_r1_label_marginals.tsv"
locus_manifest_file <- "example/22B02_r1/22B02_r1_locus_manifest.tsv"
outgroup            <- c("bog_022G05", "bog_022F01", "bog_022G01")

allelic_groups <- list(
  alpha = c("A", "a"),
  beta  = c("B", "b"),
  gamma = c("G")
)

# Read tree with annotations, then reroot/rebalance its phylo slot.
tree <- read.beast(map_tree_file)
tree@phylo <- root_and_rebalance_tree(
  tree@phylo,
  outgroup = outgroup,
  root_edge_split = root_edge_split
)

# Plot tree
p <- ggtree(tree, size = 0.35) +
  geom_tiplab(size = 2.2, align = TRUE, linetype = 0, offset = 0.002) +
  geom_treescale(x = 0, y = 0, width = 0.01, fontsize = 2.8)

# Add node posterior points.
p <- p + geom_point2(
  aes(subset = !isTip & !is.na(posterior), color = posterior),
  size = 1.6
)

# Add phasing marginals panel.
p <- p + add_phase_marginals(
  marginal_file = marginal_file,
  locus_manifest_file = locus_manifest_file,
  allelic_groups = allelic_groups,
  allele_mode = allele_mode,
  tip_label_padding = tip_label_padding
)

# Add
p <- p +
  scale_color_gradientn(
    colours = posterior_palette,
    limits = c(0, 1),
    aesthetics = c("colour", "fill"),
    na.value = "transparent",
    name = "posterior probability"
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.03, 0.97),
    legend.justification = c(0, 1),
    plot.margin = margin(t = 25, r = 10, b = 10, l = 10),
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text.x.top = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 0),
    axis.ticks.x.top = element_blank(),
    axis.title.x.top = element_blank()
  ) +
  set_phase_panel_width(phasing_panel_prop = phasing_panel_prop)

p

# Save a high-resolution PNG (helps avoid seam/banding artifacts in split cells)
ggsave(
  filename = plot_filename,
  plot = p,
  device = ragg::agg_png,
  width = 14,
  height = 9,
  units = "in",
  dpi = 600,
  bg = "white"
)



