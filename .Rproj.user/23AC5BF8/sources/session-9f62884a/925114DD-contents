library(gtools)

# input/output for trees
tree_file   <- "scratch/22B02_r1/22B02_r1trees.trees"
output_file <- "scratch/22B02_r1/22B02_r1trees_phase_labeled.trees"

# phasing files
# phase_files <- c("L10"  = "scratch/22B02_r1/22B02_r1L10_phase.log",
#                  "L100" = "scratch/22B02_r1/22B02_r1L100_phase.log",
#                  "L103" = "scratch/22B02_r1/22B02_r1L103_phase.log")

phase_files <- mixedsort(list.files("scratch/22B02_r1", full.names = TRUE, pattern = "phase.log"))
phase_names <- sapply(strsplit(phase_files, "_"), function(x) x[3])
phase_names <- gsub("r1", "", phase_names)
names(phase_files) <- phase_names

phaseLabelTrees(tree_file, phase_files, output_file)


