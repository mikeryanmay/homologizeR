# Explicit phase files (10 loci).
phase_files <- c(
  "L9"   = "example/22B02_r1/22B02_r1L9_phase.log",
  "L118" = "example/22B02_r1/22B02_r1L118_phase.log",
  "L148" = "example/22B02_r1/22B02_r1L148_phase.log",
  "L215" = "example/22B02_r1/22B02_r1L215_phase.log",
  "L226" = "example/22B02_r1/22B02_r1L226_phase.log",
  "L228" = "example/22B02_r1/22B02_r1L228_phase.log",
  "L240" = "example/22B02_r1/22B02_r1L240_phase.log",
  "L288" = "example/22B02_r1/22B02_r1L288_phase.log",
  "L329" = "example/22B02_r1/22B02_r1L329_edited_phase.log",
  "L350" = "example/22B02_r1/22B02_r1L350_phase.log"
)

# File with sampled trees
tree_file <- "example/22B02_r1/22B02_r1trees.trees"

# Define which labels can swap.
allelic_tips <- list(
  alpha = c("A", "a"),
  beta  = c("B", "b"),
  gamma = c("G")
)

# End-to-end file workflow:
# - writes:
#   <stem>_partition_relabeled.log
#   <stem>_trees_relabeled.trees
#   <stem>_relabel_map.tsv
#   <stem>_locus_manifest.tsv
#   <stem>_label_marginals.tsv
output_prefix <- "example/22B02_r1/22B02_r1"

# Do the relabeling
relabelPhaseAndTrees(
  phase_files = phase_files,
  allelic_tips = allelic_tips,
  tree_file = tree_file,
  output_prefix = output_prefix,
  burnin = 0.25
)
