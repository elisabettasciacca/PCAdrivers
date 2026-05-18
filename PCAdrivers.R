# PCAdrivers
# [code adapted from https://github.com/KatrionaGoldmann/BioOutputs/blob/master/R/bio_drivers.R]
# Computes associations between PCs and variables, then visualises them as a
# tile heatmap.
# Standalone (non-federated) equivalent of dsPCAdrivers + dsPCAdriversClient.
#
# Dependencies: stats (base), ggplot2
# 
# PCAdrivers --------------------------------------------------------------------

#' Plot drivers of variation
#'
#' Computes associations between principal components and a set of variables,
#' then visualises the results as a \eqn{-\log_{10}(p)} tile heatmap.
#'
#' The function expects a pre-computed matrix of PC scores as input, so PCA
#' can be performed with any tool before calling \code{PCAdrivers()}. Common
#' workflows:
#' \itemize{
#'   \item \strong{stats} (\code{prcomp}): pass \code{pca$x}, optionally
#'     subsetting columns first, e.g. \code{pca$x[, 1:10]}.
#'   \item \strong{PCAtools} (\code{pca()}): pass \code{p$rotated}.
#'   \item Any other tool: pass any numeric matrix with samples as rows and
#'     PCs as columns. Column names are used as PC labels; if absent they are
#'     set to \code{PC1}, \code{PC2}, \ldots
#' }
#'
#' For numeric variables a correlation test is used (Pearson or Spearman).
#' For categorical variables a one-way test is used (ANOVA or Kruskal-Wallis).
#' Variables with zero variance, fewer non-NA values than \code{na_drop_threshold},
#' or that look like ID columns (all values unique, non-numeric) are silently
#' dropped before testing.
#'
#' @param scores A numeric matrix or data frame of PC scores with
#'   \strong{samples as rows and PCs as columns}. Column names are used as
#'   labels on the plot. If absent they are set to \code{PC1}, \code{PC2},
#'   \ldots
#'   \itemize{
#'     \item From \code{stats::prcomp}: pass \code{pca$x} or a column subset
#'       such as \code{pca$x[, 1:10]}.
#'     \item From \code{PCAtools::pca}: pass \code{p$rotated}.
#'   }
#' @param vars A data frame containing the variables to associate with the PCs
#'   (e.g. clinical covariates). Must have the same number of rows as
#'   \code{scores}.
#' @param parametric Logical. Use parametric tests (Pearson, ANOVA)? If
#'   \code{FALSE}, uses Spearman and Kruskal-Wallis. Default \code{TRUE}.
#' @param na_drop_threshold Integer. Minimum number of non-NA values required
#'   for a variable to be included. Default \code{4}.
#' @param p_adj Character. P value adjustment method passed to
#'   \code{\link[stats]{p.adjust}}. One of \code{"holm"}, \code{"hochberg"},
#'   \code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"},
#'   \code{"fdr"}, or \code{"none"}. Default \code{NULL} (no adjustment).
#' @param sig_cutoff Numeric. Significance threshold for outlining tiles.
#'   Default \code{0.05}.
#' @param max_col Numeric. Maximum value for the colour scale. If \code{NULL},
#'   uses the observed maximum. Default \code{NULL}.
#' @param title Character. Plot title. Default \code{"Drivers of Variation"}.
#' @param legend Character. Legend position passed to \code{ggplot2}. Default
#'   \code{"right"}.
#' @param label Logical. Print \eqn{-\log_{10}(p)} values on tiles? Default
#'   \code{FALSE}.
#' @param transpose_plot Logical. Swap axes (PCs on y, variables on x)? Default
#'   \code{FALSE}.
#' @param drop_insignificant_x Logical. Drop PCs with no significant
#'   associations from the plot? Default \code{FALSE}.
#' @param drop_insignificant_y Logical. Drop variables with no significant
#'   associations from the plot? Default \code{FALSE}.
#' @param return_data Logical. Return the results data frame instead of a plot?
#'   Default \code{FALSE}.
#'
#' @return When \code{return_data = FALSE} (default), a \code{ggplot2} object.
#'   When \code{return_data = TRUE}, a list with:
#'   \describe{
#'     \item{\code{results}}{Data frame with columns \code{Feature}, \code{PC},
#'       \code{pvalue}, \code{Association} (\eqn{-\log_{10}(p)}), and
#'       \code{Significant}.}
#'     \item{\code{metadata}}{List with \code{n_observations}, \code{pc_names},
#'       \code{var_names}, \code{parametric}, \code{p_adj}.}
#'   }
#'
#' @examples
#' set.seed(42)
#' expr     <- matrix(rnorm(200 * 50), nrow = 200, ncol = 50)
#' colnames(expr) <- paste0("gene", seq_len(50))
#'
#' clinical <- data.frame(
#'   age   = rnorm(200, 50, 10),
#'   sex   = sample(c("M", "F"), 200, replace = TRUE),
#'   batch = factor(sample(1:3, 200, replace = TRUE))
#' )
#'
#' # Using stats::prcomp
#' pca <- prcomp(expr, center = TRUE, scale. = FALSE)
#' PCAdrivers(scores = pca$x[, 1:5], vars = clinical)
#'
#' # Using PCAtools (not run)
#' # p <- PCAtools::pca(t(expr), metadata = clinical)
#' # PCAdrivers(scores = p$rotated, vars = clinical)
#'
#' PCAdrivers(scores = pca$x[, 1:5], vars = clinical,
#'            parametric = FALSE, p_adj = "BH", label = TRUE)
#'
#' @importFrom stats cor.test lm anova kruskal.test p.adjust
#' @importFrom ggplot2 ggplot aes geom_tile geom_text coord_equal
#'   scale_fill_gradientn scale_colour_manual guides guide_legend labs
#'   theme_bw theme element_text
#' @export
PCAdrivers <- function(scores,
                       vars,
                       parametric           = TRUE,
                       na_drop_threshold    = 4L,
                       p_adj                = NULL,
                       sig_cutoff           = 0.05,
                       max_col              = NULL,
                       title                = "Drivers of Variation",
                       legend               = "right",
                       label                = FALSE,
                       transpose_plot       = FALSE,
                       drop_insignificant_x = FALSE,
                       drop_insignificant_y = FALSE,
                       return_data          = FALSE) {
  
  # Input validation -------------------------------------------------------------
  
  if (!is.matrix(scores) && !is.data.frame(scores)) {
    stop("'scores' must be a numeric matrix or data frame of PC scores ",
         "(samples × PCs).", call. = FALSE)
  }
  scores <- as.matrix(scores)
  if (!is.numeric(scores)) {
    stop("'scores' must be numeric.", call. = FALSE)
  }
  if (!is.data.frame(vars) && !is.matrix(vars)) {
    stop("'vars' must be a data frame or matrix.", call. = FALSE)
  }
  if (nrow(scores) != nrow(vars)) {
    stop("'scores' and 'vars' must have the same number of rows (samples).",
         call. = FALSE)
  }
  if (!is.numeric(sig_cutoff) || length(sig_cutoff) != 1 ||
      sig_cutoff <= 0 || sig_cutoff >= 1) {
    stop("'sig_cutoff' must be a single numeric value in (0, 1).", call. = FALSE)
  }
  
  # PC names: use column names if present, otherwise PC1, PC2, ...
  if (is.null(colnames(scores))) {
    colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
  }
  pc_names <- colnames(scores)
  
  # Clean variables
  vars_clean <- .clean_vars(vars, na_threshold = na_drop_threshold)
  
  if (ncol(vars_clean) == 0) {
    stop("No valid variables remaining after filtering.", call. = FALSE)
  }
  
  var_names <- colnames(vars_clean)
  
  # Compute associations
  combinations <- expand.grid(
    feature = factor(var_names, levels = var_names),
    pc      = factor(pc_names,  levels = pc_names),
    stringsAsFactors = FALSE
  )
  combinations$feature <- as.character(combinations$feature)
  combinations$pc      <- as.character(combinations$pc)
  
  pvals <- mapply(
    FUN       = .single_association,
    feat_name = combinations$feature,
    pc_name   = combinations$pc,
    MoreArgs  = list(scores = scores, vars = vars_clean,
                     parametric = parametric),
    SIMPLIFY  = TRUE
  )
  
  # Build results ----------------------------------------------------------------
  
  results <- data.frame(
    Feature = factor(combinations$feature, levels = var_names),
    PC      = factor(combinations$pc,      levels = pc_names),
    pvalue  = pvals,
    stringsAsFactors = FALSE
  )
  
  if (!is.null(p_adj)) {
    results$pvalue <- p.adjust(results$pvalue, method = p_adj)
  }
  
  results$Association <- -log10(results$pvalue)
  results$Significant <- !is.na(results$pvalue) & results$pvalue <= sig_cutoff
  results$Feature     <- as.character(results$Feature)
  results$PC          <- as.character(results$PC)
  
  # Return data or plot ----------------------------------------------------------
  
  if (return_data) {
    return(list(
      results  = results,
      metadata = list(
        n_observations = nrow(scores),
        pc_names       = pc_names,
        var_names      = var_names,
        parametric     = parametric,
        p_adj          = p_adj
      )
    ))
  }
  
  if (drop_insignificant_x) {
    keep_pc <- unique(results$PC[results$Significant])
    results  <- results[results$PC %in% keep_pc, ]
  }
  if (drop_insignificant_y) {
    keep_feat <- unique(results$Feature[results$Significant])
    results   <- results[results$Feature %in% keep_feat, ]
  }
  
  .build_plot(results, sig_cutoff = sig_cutoff, p_adj = p_adj,
              max_col = max_col, title = title, legend = legend,
              transpose_plot = transpose_plot, label = label)
}


# Internal helpers -------------------------------------------------------------

# Remove zero-variance, data-sparse, and ID-like columns from a variables data
# frame before association testing.
.clean_vars <- function(vars, na_threshold) {
  
  # Replace Inf / -Inf with NA
  vars[] <- lapply(vars, function(x) {
    if (is.numeric(x)) replace(x, is.infinite(x), NA) else x
  })
  
  # Drop zero-variance columns
  has_var <- vapply(vars, function(x) length(unique(x[!is.na(x)])) > 1, logical(1))
  if (any(!has_var)) {
    message("Dropping zero-variance column(s): ",
            paste(colnames(vars)[!has_var], collapse = ", "))
  }
  vars <- vars[, has_var, drop = FALSE]
  if (ncol(vars) == 0) return(vars)
  
  # Drop columns with too few non-NA values
  enough_data <- colSums(!is.na(vars)) >= na_threshold
  if (any(!enough_data)) {
    message("Dropping column(s) with fewer than ", na_threshold, " non-NA values: ",
            paste(colnames(vars)[!enough_data], collapse = ", "))
  }
  vars <- vars[, enough_data, drop = FALSE]
  if (ncol(vars) == 0) return(vars)
  
  # Drop ID-like columns: non-numeric where every non-missing value is unique
  is_numeric <- vapply(vars, is.numeric, logical(1))
  all_unique  <- vapply(vars, function(x) {
    nn <- sum(!is.na(x))
    nn > 0 && length(unique(x[!is.na(x)])) == nn
  }, logical(1))
  is_id <- all_unique & !is_numeric
  
  if (any(is_id)) {
    message("Dropping ID-like column(s): ",
            paste(colnames(vars)[is_id], collapse = ", "))
  }
  vars[, !is_id, drop = FALSE]
}


# Compute the association p value between one PC and one variable.
.single_association <- function(feat_name, pc_name, scores, vars,
                                parametric) {
  
  pc_vals   <- scores[, pc_name]
  feat_vals <- vars[[feat_name]]
  
  complete  <- !is.na(pc_vals) & !is.na(feat_vals)
  
  if (sum(complete) < 3L) {
    message("Skipping ", feat_name, " ~ ", pc_name,
                         ": fewer than 3 complete cases.")
    return(NA_real_)
  }
  
  pc_vals   <- pc_vals[complete]
  feat_vals <- feat_vals[complete]
  
  if (is.numeric(feat_vals)) {
    method     <- if (parametric) "pearson" else "spearman"
    test       <- cor.test(pc_vals, feat_vals, method = method)
    return(test$p.value)
  }
  
  # Categorical variable
  feat_vals <- as.factor(feat_vals)
  
  if (nlevels(feat_vals) < 2L) {
    message("Skipping ", feat_name, " ~ ", pc_name,
                         ": fewer than 2 levels after subsetting.")
    return(NA_real_)
  }
  
  if (parametric) {
    fit <- lm(pc_vals ~ feat_vals)
    return(anova(fit)[1L, "Pr(>F)"])
  } else {
    return(kruskal.test(pc_vals ~ feat_vals)$p.value)
  }
}


# Build the ggplot2 tile heatmap.
.build_plot <- function(results, sig_cutoff, p_adj, max_col, title, legend,
                        transpose_plot, label) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to produce the plot. ",
         "Install it with install.packages('ggplot2').", call. = FALSE)
  }
  
  if (is.null(max_col)) {
    max_col <- max(ceiling(results$Association[is.finite(results$Association)]),
                   na.rm = TRUE)
    if (!is.finite(max_col) || max_col == 0) max_col <- 1
  }
  
  results$x <- if (transpose_plot) results$Feature else results$PC
  results$y <- if (transpose_plot) results$PC       else results$Feature
  
  leg_lab <- if (!is.null(p_adj)) {
    expression(-log[10](italic(p)[adj]))
  } else {
    expression(-log[10](italic(p)))
  }
  
  p <- ggplot2::ggplot(
    results,
    ggplot2::aes(x = x, y = y,
                 fill   = Association,
                 colour = Significant)
  ) +
    ggplot2::geom_tile(linewidth = 0.8, width = 0.9, height = 0.9) +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_gradientn(
      colours = c("white", "dodgerblue1", "dodgerblue3", "dodgerblue4"),
      name    = leg_lab,
      limits  = c(0, max_col),
      na.value = "grey80"
    ) +
    ggplot2::scale_colour_manual(
      values = c("TRUE" = "black", "FALSE" = "grey88"),
      labels = c(
        "TRUE"  = paste(ifelse(is.null(p_adj), "p", "p adj"), "\u2264", sig_cutoff),
        "FALSE" = paste(ifelse(is.null(p_adj), "p", "p adj"), ">",      sig_cutoff)
      ),
      name = ""
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(override.aes = list(fill = "white"))
    ) +
    ggplot2::labs(title = title, x = "", y = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title     = ggplot2::element_text(hjust = 0.5),
      axis.text.x    = ggplot2::element_text(
        angle = if (transpose_plot) 315 else 0,
        hjust = if (transpose_plot) 0   else 0.5
      ),
      legend.position = legend
    )
  
  if (label) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = ifelse(is.finite(Association),
                                  round(Association, 2), "NA")),
      colour = "black", size = 2.8
    )
  }
  
  p
}

