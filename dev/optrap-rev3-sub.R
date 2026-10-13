# Reviewer of lane optrap, final check: is predictor_resid_sd()'s
# subsample of evenly spaced rows representative? On small designs the
# cap (internal option frmtmb.unit_ls_cells) is lowered so that the
# subsample is taken (stride about 10), and the unit is compared with
# the one from every row (cap Inf), frame only. Designs: rows sorted by
# group with large group effects; rows sorted by the response; a factor
# level present only in the last rows, with a large effect; a factor
# whose levels cycle with the stride; heteroscedastic rows sorted by x.
.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
pu <- frmtmb:::predictor_units
unit_at <- function(fo, d, cells) {
  old <- options(frmtmb.unit_ls_cells = cells)
  on.exit(options(old))
  u <- pu(frm(fo, data = d, dry_run = "frame"))
  if (length(u)) u[[1]] else 1
}
X <- list()
for (s in 1:20) {
  set.seed(s)
  n <- 1000
  designs <- list()
  d <- data.frame(g = factor(rep(1:100, each = 10)), x = rnorm(n))
  d$y <- 1 + 0.5 * d$x + rnorm(100, 0, 5)[d$g] + rnorm(n)
  designs$sorted_group <- list(fo = y ~ x + (1 | g), d = d)
  d2 <- d[order(d$y), ]
  designs$sorted_y <- list(fo = y ~ x + (1 | g), d = d2)
  d3 <- data.frame(x = rnorm(n), f = factor(c(rep("a", n - 5), rep("z", 5))))
  d3$y <- 1 + 0.5 * d3$x + 1e4 * (d3$f == "z") + rnorm(n)
  designs$last_level <- list(fo = y ~ x + f, d = d3)
  d4 <- data.frame(x = rnorm(n), f = factor(rep(1:10, length.out = n)))
  d4$y <- 1 + 0.5 * d4$x + c(0, 50, -50, 20, 0, 10, -30, 5, 40, -10)[d4$f] +
    rnorm(n, 0, c(1, 5, 1, 5, 1, 5, 1, 5, 1, 5)[d4$f])
  designs$cycle_level <- list(fo = y ~ x + f, d = d4)
  d5 <- data.frame(x = sort(rnorm(n)))
  d5$y <- 1 + 0.5 * d5$x + rnorm(n, 0, exp(1.5 * d5$x))
  designs$sorted_hetero <- list(fo = y ~ x, d = d5)
  for (nm in names(designs)) {
    ds <- designs[[nm]]
    p <- ncol(model.matrix(lme4::nobars(ds$fo), ds$d))
    X[[length(X) + 1L]] <- data.frame(
      design = nm, seed = s, full = unit_at(ds$fo, ds$d, Inf),
      sub = unit_at(ds$fo, ds$d, 100 * p))
  }
}
X <- do.call(rbind, X)
X$log2ratio <- log2(X$sub / X$full)
write.table(X, "dev/optrap-rev3-log/sub.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE)
for (nm in unique(X$design)) {
  x <- X[X$design == nm, ]
  cat(sprintf("%-14s | full units %s | subsample units %s | |log2 ratio| > 1: %d of 20, max %g\n",
              nm, paste(sort(unique(x$full)), collapse = ","),
              paste(sort(unique(x$sub)), collapse = ","),
              sum(abs(x$log2ratio) > 1), max(abs(x$log2ratio))))
}
