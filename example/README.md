# README

This is a working example of plotting homologizer phasing results on a MAP tree.

## Installation

Open up the homogolizeR.Rproj in Rstudio and click on the `Build > Install Package` in the menu.

## Step 1: relabeling output

The first step is to relabel the outputs so there is a consistent meaning of "A", "B", etc across MCMC samples.
An example of this is shown in `example/relabel_outputs.R`.
The main thing you need to do is provide filepaths to the phasing files and tree files that are going to be relabeled.
Notice that the phasing file paths are named (e.g., `"L9" = "example/22B02_r1/22B02_r1L9_phase.log"`), which tells the code what the loci are called in downstream steps.
There is a burnin which I have set to 0.25, but adjust as necessary.
Running this code will generate a few output files, most of which you can ignore.
The main output file you care about is the `*_relabeled.trees` file, which you will use to get a MAP tree.

## Step 2: estimate the MAP tree

Now that the tree samples are relabeled, you need to get the MAP.
I use revbayes for this, but you could do anything really.
I have a `example/summarize.Rev` script that makes the MAP tree.

## Step 3: plot the phasing probabilities

The script `example/plot_phase.R` plots the phasing probabilities next to the MAP tree.
The main arguments are at the top of this script, which provides controls for plot panel sizes (`phasing_panel_prop`), tip-label padding (basically how much space is between the tip labels and the phasing probabilities, `tip_label_padding`), and how to divide the root branch lengths to "balance" how far forward/back the outgroup is (`root_edge_split`).
There is also an `allele_mode` that controls how to merge probabilities between allelic tips (the options are `none`, `vertical`, and `horizontal`); I like `vertical` as a default.
`plot_filename` says where to save the plot to (adjust the plot dimensions and other properties at the bottom of the script).
Finally, there's a place where you define the posterior probability color ramp, `posterior_palette`.
After that, you define where to find the tree files, marginal files, and sample manifest (all generated in step 1, with default names that are reflected here).


