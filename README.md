# PCAdrivers

R function to visualise which variables drive variation in your
data. It computes associations between pre-computed principal components and a
set of variables (e.g. clinical covariates), then displays the results as a
−log₁₀(p) tile heatmap.

Numeric variables are tested with Pearson correlation (or Spearman if
`parametric = FALSE`). Categorical variables are tested with one-way ANOVA
(or Kruskal-Wallis). Tiles are outlined in black when the association reaches
the significance threshold.

## Installation

No package installation needed. Just source the script:

```r
source("https://raw.githubusercontent.com/elisabettasciacca/PCAdrivers/main/PCAdrivers.R")
```

The only non-base dependency is [ggplot2](https://ggplot2.tidyverse.org/),
which is needed only to produce the plot (`return_data = FALSE`).

```r
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
```

## Usage

`PCAdrivers()` expects a matrix of PC scores as input — PCA is computed
**outside** the function, so you are free to use any tool and any
preprocessing you like. Only the scores matrix is passed in.

```r
PCAdrivers(
  scores,                    # numeric matrix: samples × PCs
  vars,                      # data frame of variables to test
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
  return_data          = FALSE
)
```

### Scores input

| Source | What to pass |
|---|---|
| `stats::prcomp` | `pca$x` or a column subset, e.g. `pca$x[, 1:10]` |
| `PCAtools::pca` | `p$rotated` |
| Any other tool | Any numeric matrix with samples as rows and PCs as columns |

Column names of `scores` are used as labels on the plot. If absent they are
set to `PC1`, `PC2`, …

## Examples

### Example 1 — basic usage with `prcomp`

```r
source("PCAdrivers.R")

set.seed(42)
expr <- matrix(rnorm(200 * 500), nrow = 200, ncol = 500)
colnames(expr) <- paste0("gene", seq_len(500))

clinical <- data.frame(
  age      = rnorm(200, 50, 10),
  sex      = sample(c("M", "F"), 200, replace = TRUE),
  batch    = factor(sample(1:3, 200, replace = TRUE)),
  lab_site = factor(sample(c("London", "Rome", "Berlin"), 200, replace = TRUE))
)

pca <- prcomp(expr, center = TRUE, scale. = FALSE)

PCAdrivers(scores = pca$x[, 1:10], vars = clinical)
```

### Example 2 — non-parametric tests, p value adjustment, and returning data

Useful when variables are non-normally distributed or the data contain
outliers. `return_data = TRUE` gives back the full results table alongside
plot metadata, so you can do further downstream work.

```r
res <- PCAdrivers(
  scores      = pca$x[, 1:10],
  vars        = clinical,
  parametric  = FALSE,   # Spearman + Kruskal-Wallis
  p_adj       = "BH",    # Benjamini-Hochberg correction
  sig_cutoff  = 0.05,
  return_data = TRUE
)

# results table
head(res$results)
#>   Feature  PC      pvalue Association Significant
#> 1     age PC1 0.623941760   0.2046851       FALSE
#> 2     sex PC1 0.450981189   0.3459871       FALSE
#> 3   batch PC1 0.007293785   2.1369540        TRUE
#> ...

# metadata
res$metadata$n_observations
```

> **Note:** when `return_data = TRUE` no plot is produced. Use the default
> `return_data = FALSE` to get the heatmap, or build your own from
> `res$results`.

### Example 3 — with PCAtools

```r
library(PCAtools)

# PCAtools expects features × samples
p <- pca(t(expr), metadata = clinical, removeVar = 0.1)

PCAdrivers(
  scores     = p$rotated,   # samples × PCs data frame
  vars       = clinical,
  sig_cutoff = 0.01,
  p_adj      = "bonferroni",
  label      = TRUE         # print -log10(p) values on tiles
)
```

## Arguments

| Argument | Type | Default | Description |
|---|---|---|---|
| `scores` | matrix / data frame | — | PC scores, samples × PCs |
| `vars` | data frame | — | Variables to associate with PCs |
| `parametric` | logical | `TRUE` | Pearson + ANOVA (`TRUE`) or Spearman + Kruskal-Wallis (`FALSE`) |
| `na_drop_threshold` | integer | `4` | Minimum non-NA values required per variable |
| `p_adj` | character | `NULL` | P value adjustment method (see `?p.adjust`); `NULL` for none |
| `sig_cutoff` | numeric | `0.05` | Significance threshold for tile outlines |
| `max_col` | numeric | `NULL` | Maximum value for colour scale; `NULL` uses observed maximum |
| `title` | character | `"Drivers of Variation"` | Plot title |
| `legend` | character | `"right"` | Legend position (`"right"`, `"bottom"`, `"left"`, `"top"`) |
| `label` | logical | `FALSE` | Print −log₁₀(p) values on tiles |
| `transpose_plot` | logical | `FALSE` | Swap axes (PCs on y, variables on x) |
| `drop_insignificant_x` | logical | `FALSE` | Drop PCs with no significant associations |
| `drop_insignificant_y` | logical | `FALSE` | Drop variables with no significant associations |
| `return_data` | logical | `FALSE` | Return results list instead of printing the plot |

## Output

When `return_data = FALSE` (default), a `ggplot2` object is returned and
printed.

When `return_data = TRUE`, a list with two elements:

- `$results` — data frame with columns `Feature`, `PC`, `pvalue`,
  `Association` (−log₁₀(p)), `Significant`
- `$metadata` — list with `n_observations`, `pc_names`, `var_names`,
  `parametric`, `p_adj`

## Variable filtering

Before testing, the following columns are automatically dropped from `vars`:

- zero-variance columns (only one unique value)
- columns with fewer non-NA values than `na_drop_threshold`
- ID-like columns (non-numeric, all values unique)

## Related packages

This script is the standalone equivalent of the
[dsPCAdrivers](https://github.com/YOUR_USERNAME/dsPCAdrivers) +
[dsPCAdriversClient](https://github.com/YOUR_USERNAME/dsPCAdriversClient)
package pair, which implements the same analysis in a federated
[DataSHIELD](https://www.datashield.org/) setting.

## Licence

GPL-3
