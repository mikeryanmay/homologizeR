#' @importFrom ggplot2 ggplot_add
NULL

#' Read locus manifest table
#'
#' Expects columns `locus_copy`, `locus`, and `copy`.
#' Returns unique rows by `locus_copy`.
#'
#' @param path File path to manifest TSV.
#' @return A data frame.
read_locus_manifest <- function(path) {
  m <- utils::read.delim(path, sep = "\t", header = TRUE, check.names = FALSE)
  required <- c("locus_copy", "locus", "copy")
  missing <- setdiff(required, names(m))
  if (length(missing) > 0L) {
    stop("locus manifest is missing required columns: ", paste(missing, collapse = ", "))
  }
  m$locus_copy <- as.character(m$locus_copy)
  m$locus <- as.character(m$locus)
  m$copy <- as.character(m$copy)
  m[!duplicated(m$locus_copy), , drop = FALSE]
}

#' Root and rebalance a phylogeny
#'
#' Roots the tree on an outgroup and rebalances the two root-descendant
#' edge lengths using `root_edge_split`.
#'
#' @param tree A phylo object.
#' @param outgroup Character vector of outgroup tip labels.
#' @param root_edge_split Proportion in (0,1) assigned to first root edge.
#' @return A rooted/rebalanced phylo object.
root_and_rebalance_tree <- function(tree, outgroup, root_edge_split = 0.4) {
  tree <- ape::root(tree, outgroup = outgroup, resolve.root = TRUE)

  if (!is.numeric(root_edge_split) || length(root_edge_split) != 1L ||
      !is.finite(root_edge_split) || root_edge_split <= 0 || root_edge_split >= 1) {
    stop("root_edge_split must be a single numeric value in (0, 1).")
  }

  root_node <- ape::Ntip(tree) + 1L
  root_edge_idx <- which(tree$edge[, 1] == root_node)
  total_root_len <- sum(tree$edge.length[root_edge_idx])
  tree$edge.length[root_edge_idx[1]] <- total_root_len * root_edge_split
  tree$edge.length[root_edge_idx[2]] <- total_root_len * (1 - root_edge_split)
  tree
}

best_tip_assignment <- function(prob_mat) {
  k <- nrow(prob_mat)
  m <- ncol(prob_mat)
  best_score <- -Inf
  best_assign <- rep(NA_integer_, k)

  recurse <- function(i, used_cols, assign_cols, score) {
    if (i > k) {
      if (score > best_score) {
        best_score <<- score
        best_assign <<- assign_cols
      }
      return(invisible(NULL))
    }
    avail <- setdiff(seq_len(m), used_cols)
    for (j in avail) recurse(i + 1L, c(used_cols, j), c(assign_cols, j), score + prob_mat[i, j])
    invisible(NULL)
  }

  recurse(1L, integer(0), integer(0), 0)
  best_assign
}

.apply_phase_marginals <- function(
  p,
  marginal_file,
  locus_manifest_file,
  allelic_groups,
  allele_mode = c("none", "vertical", "horizontal"),
  panel = "Marginals",
  tip_label_padding = 0.03
) {
  allele_mode <- match.arg(allele_mode)

  marg <- utils::read.delim(marginal_file, sep = "\t", header = TRUE, check.names = FALSE)

  manifest <- read_locus_manifest(locus_manifest_file)
  idx <- match(manifest$locus_copy, marg$locus_copy)
  keep <- !is.na(idx)
  idx <- idx[keep]
  manifest_use <- manifest[keep, , drop = FALSE]
  marg <- marg[idx, , drop = FALSE]
  marg$locus <- manifest_use$locus
  marg$copy <- manifest_use$copy
  loci <- unique(manifest$locus)

  p <- p + ggplot2::expand_limits(x = max(p$data$x, na.rm = TRUE) + tip_label_padding)
  tip_y <- setNames(p$data$y[!is.na(p$data$label)], p$data$label[!is.na(p$data$label)])

  label_rows <- unlist(allelic_groups, use.names = FALSE)
  if (allele_mode == "none") {
    row_members <- stats::setNames(as.list(label_rows), label_rows)
  } else {
    row_members <- allelic_groups
  }
  row_names <- names(row_members)
  row_y_mid <- vapply(row_members, function(x) mean(tip_y[x]), numeric(1))
  row_id <- vapply(row_names, function(g) {
    mem <- row_members[[g]]
    mem[which.min(abs(tip_y[mem] - row_y_mid[[g]]))]
  }, character(1))

  cells <- list()
  split_cells <- list()
  mid_text <- list()

  for (rnm in row_names) {
    for (loc in loci) {
      subm <- marg[marg$locus == loc, , drop = FALSE]
      mem <- row_members[[rnm]]
      mem <- mem[order(tip_y[mem], decreasing = TRUE)]

      if (allele_mode == "none") {
        vals <- subm[[mem[1]]]
        if (length(vals) == 0L || all(is.na(vals))) next
        mx <- max(vals, na.rm = TRUE)
        if (!is.finite(mx) || mx <= 0) next
        j <- which(vals == mx)
        win <- paste(sort(unique(as.character(subm$copy[j]))), collapse = "/")
        cells[[length(cells) + 1L]] <- data.frame(
          id = mem[1], row = rnm, locus = loc, x_plot = match(loc, loci),
          prob = mx, copy = win, stringsAsFactors = FALSE
        )
        next
      }

      vals <- rowSums(subm[, mem, drop = FALSE])
      pos_idx <- which(!is.na(vals) & vals > 0)
      if (length(vals) == 0L || all(is.na(vals))) next

      if (length(mem) == 1L) {
        mx <- max(vals, na.rm = TRUE)
        if (!is.finite(mx) || mx <= 0) next
        j <- which(vals == mx)
        win <- paste(sort(unique(as.character(subm$copy[j]))), collapse = "/")
        cells[[length(cells) + 1L]] <- data.frame(
          id = mem[1], row = rnm, locus = loc, x_plot = match(loc, loci),
          prob = mx, copy = win, stringsAsFactors = FALSE
        )
        next
      }

      if (length(pos_idx) <= 1L) {
        if (length(pos_idx) == 1L) {
          j <- pos_idx[1]
          span_prob <- vals[j]
          span_copy <- as.character(subm$copy[j])
          mid_text[[length(mid_text) + 1L]] <- data.frame(
            id = row_id[[rnm]], row = rnm, locus = loc, x_plot = match(loc, loci),
            nudge_y = row_y_mid[[rnm]] - tip_y[row_id[[rnm]]], copy = span_copy,
            stringsAsFactors = FALSE
          )
        } else {
          span_prob <- NA_real_
        }
        for (tip in mem) {
          cells[[length(cells) + 1L]] <- data.frame(
            id = tip, row = rnm, locus = loc, x_plot = match(loc, loci),
            prob = span_prob, copy = "", stringsAsFactors = FALSE
          )
        }
      } else if (allele_mode == "horizontal") {
        ord <- order(vals, decreasing = TRUE, na.last = NA)
        k <- min(length(ord), length(mem))
        sel <- ord[seq_len(k)]
        tip_prob_mat <- as.matrix(subm[sel, mem, drop = FALSE])
        assign_cols <- best_tip_assignment(tip_prob_mat)

        assigned_copy_idx <- rep(NA_integer_, length(mem))
        for (ii in seq_len(k)) assigned_copy_idx[assign_cols[ii]] <- sel[ii]
        for (i in seq_along(mem)) {
          has_copy <- !is.na(assigned_copy_idx[i])
          copy_row <- assigned_copy_idx[i]
          cells[[length(cells) + 1L]] <- data.frame(
            id = mem[i], row = rnm, locus = loc, x_plot = match(loc, loci),
            prob = if (has_copy) as.numeric(subm[copy_row, mem[i]]) else NA_real_,
            copy = if (has_copy) as.character(subm$copy[copy_row]) else "",
            stringsAsFactors = FALSE
          )
        }
      } else if (allele_mode == "vertical") {
        ord <- order(vals, decreasing = TRUE, na.last = NA)
        k <- min(length(ord), length(mem))
        sel <- ord[seq_len(k)]
        sel_copy_num <- as.integer(sub("^c", "", subm$copy[sel]))
        sel <- sel[order(sel_copy_num)]
        sub_w <- 0.92 / k
        x_off <- (seq_len(k) - (k + 1) / 2) * sub_w
        ny <- row_y_mid[[rnm]] - tip_y[row_id[[rnm]]]
        for (ii in seq_len(k)) {
          for (tip in mem) {
            split_cells[[length(split_cells) + 1L]] <- data.frame(
              id = tip, row = rnm, locus = loc,
              x_plot = match(loc, loci) + x_off[ii],
              prob = vals[sel[ii]], split_k = k, stringsAsFactors = FALSE
            )
          }
          mid_text[[length(mid_text) + 1L]] <- data.frame(
            id = row_id[[rnm]], row = rnm, locus = loc,
            x_plot = match(loc, loci) + x_off[ii], nudge_y = ny,
            copy = as.character(subm$copy[sel[ii]]), stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (length(cells) > 0L) {
    heat_df <- do.call(rbind, cells)
    tile_df <- heat_df[!is.na(heat_df$prob) & heat_df$prob > 0, , drop = FALSE]
    text_df <- heat_df[heat_df$copy != "" & !is.na(heat_df$prob) & heat_df$prob > 0, , drop = FALSE]

    if (nrow(tile_df) > 0L) {
      p <- ggtree::facet_plot(
        p, panel = panel, data = tile_df, geom = ggplot2::geom_tile,
        mapping = ggplot2::aes(x = x_plot, fill = prob),
        height = 1.02, width = 0.96, colour = NA, linewidth = 0
      )
    }
    if (nrow(text_df) > 0L) {
      p <- ggtree::facet_plot(
        p, panel = panel, data = text_df, geom = ggplot2::geom_text,
        mapping = ggplot2::aes(x = x_plot, label = copy), size = 2.3
      )
    }
  }

  if (length(split_cells) > 0L) {
    split_df <- do.call(rbind, split_cells)
    for (kk in sort(unique(split_df$split_k))) {
      subd <- split_df[split_df$split_k == kk, , drop = FALSE]
      p <- ggtree::facet_plot(
        p, panel = panel, data = subd, geom = ggplot2::geom_tile,
        mapping = ggplot2::aes(x = x_plot, fill = prob),
        height = 1.02, width = (0.92 / kk) + 0.01, colour = NA, linewidth = 0
      )
    }
  }

  if (length(mid_text) > 0L) {
    mid_df <- do.call(rbind, mid_text)
    for (ny in sort(unique(mid_df$nudge_y))) {
      subd <- mid_df[mid_df$nudge_y == ny, , drop = FALSE]
      p <- ggtree::facet_plot(
        p, panel = panel, data = subd, geom = ggplot2::geom_text,
        mapping = ggplot2::aes(x = x_plot, label = copy), nudge_y = ny, size = 2.3
      )
    }
  }

  p <- p + ggplot2::scale_x_continuous(
    breaks = seq_along(loci),
    labels = loci,
    position = "top",
    expand = ggplot2::expansion(mult = c(0.01, 0.01))
  )

  p
}

#' Add phase-marginal heatmap panel to a ggtree plot
#'
#' Use with `+` on an existing `ggtree` object.
#'
#' @param marginal_file File path to marginal-probability table.
#' @param locus_manifest_file File path to locus manifest table.
#' @param allelic_groups Named list of allelic tip groups.
#' @param allele_mode One of `"none"`, `"vertical"`, `"horizontal"`.
#' @param panel Facet panel name.
#' @param tip_label_padding Extra x-space for tip labels.
#' @return A ggplot-addable layer object.
#' @export
add_phase_marginals <- function(
  marginal_file,
  locus_manifest_file,
  allelic_groups,
  allele_mode = c("none", "vertical", "horizontal"),
  panel = "Marginals",
  tip_label_padding = 0.03
) {
  allele_mode <- match.arg(allele_mode)
  structure(
    list(
      marginal_file = marginal_file,
      locus_manifest_file = locus_manifest_file,
      allelic_groups = allelic_groups,
      allele_mode = allele_mode,
      panel = panel,
      tip_label_padding = tip_label_padding
    ),
    class = "phase_marginals_layer"
  )
}

#' @rdname add_phase_marginals
#' @method ggplot_add phase_marginals_layer
#' @export
ggplot_add.phase_marginals_layer <- function(object, plot, object_name) {
  .apply_phase_marginals(
    p = plot,
    marginal_file = object$marginal_file,
    locus_manifest_file = object$locus_manifest_file,
    allelic_groups = object$allelic_groups,
    allele_mode = object$allele_mode,
    panel = object$panel,
    tip_label_padding = object$tip_label_padding
  )
}

#' Set relative width of phasing panel
#'
#' Use with `+` as a final step after styling.
#'
#' @param phasing_panel_prop Proportion in (0,1) for Marginals panel width.
#' @return A ggplot-addable layer object.
#' @export
set_phase_panel_width <- function(phasing_panel_prop) {
  structure(list(phasing_panel_prop = phasing_panel_prop), class = "phase_panel_width_layer")
}

apply_phase_panel_width <- function(p, phasing_panel_prop) {
  if (!is.numeric(phasing_panel_prop) || length(phasing_panel_prop) != 1L ||
      !is.finite(phasing_panel_prop) || phasing_panel_prop <= 0 || phasing_panel_prop >= 1) {
    stop("phasing_panel_prop must be a single numeric value in (0, 1).")
  }
  ggtree::facet_widths(
    p,
    c(Tree = 1 - phasing_panel_prop, Marginals = phasing_panel_prop)
  )
}

#' @rdname set_phase_panel_width
#' @method ggplot_add phase_panel_width_layer
#' @export
ggplot_add.phase_panel_width_layer <- function(object, plot, object_name) {
  apply_phase_panel_width(plot, object$phasing_panel_prop)
}
