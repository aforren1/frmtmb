# Lane wt-reunc, punch round 1 (M1): the CATEGORICAL half of the
# finite-difference route, which the cost script does not reach.
#
# A categorical fit has one linear predictor per non-reference
# category, so re_row_support() sums several Z matrices and
# re_b_batches() has to attribute a difference across all of them. A
# cumulative fit has one. The bound and the batch must change no answer
# here either, and a distributional fit (sigma ~ x with its own group
# term) is the other multi-predictor shape.
#
#   Rscript dev/reunc-fdcat.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(lib)) lib <- "C:/Users/adf44/source/r/reunc2-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")

cmp <- function(tag, a, b) {
  a <- unname(as.vector(a))
  b <- unname(as.vector(b))
  cat(sprintf("  %-34s identical %-5s max ulp %.1f\n", tag,
              identical(a, b),
              max(abs(a - b) / (.Machine$double.eps * abs(a)))))
}

set.seed(12)
ng <- 25
d <- data.frame(g = factor(rep(seq_len(ng), each = 6)),
                x = rnorm(ng * 6))
u <- rnorm(ng, 0, 0.8)
e <- function() rlogis(ng * 6)
lat2 <- 0.5 + 0.7 * d$x + u[d$g] + e()
lat3 <- -0.3 - 0.5 * d$x + u[d$g] + e()
d$y <- factor(max.col(cbind(0, lat2, lat3)), labels = c("a", "b", "c"))
fit <- suppressWarnings(frm(bf(y ~ x + (1 | g)) + categorical(), data = d))
nd <- d[c(1, 61), c("x", "g")]

gov <- frmtmb:::re_governed_b(fit)
used <- frmtmb:::re_used_b(fit, nd, "y", FALSE)
bts <- frmtmb:::re_b_batches(fit, NULL, "y", FALSE, gov)
cat("\ncategorical, ", ng, " levels: kept ", length(gov),
    ", loaded by the 2 rows ", length(used), ", batches ",
    if (is.null(bts)) "refused" else length(bts), "\n", sep = "")
cmp("2 newdata rows, bounded+batched",
    frmtmb:::fit_fd_se(fit, function(x) frmtmb:::fitted_point(x, nd),
                       b_idx = gov),
    fitted(fit, newdata = nd)[, "Est.Error", ])
cmp("in sample, batched",
    frmtmb:::fit_fd_se(fit, function(x) frmtmb:::fitted_point(x),
                       b_idx = gov),
    fitted(fit)[, "Est.Error", ])

# a group term on a second linear predictor of the SAME response, which
# is the shape re_row_support() sums over
set.seed(13)
d$h <- factor(rep(1:5, length.out = nrow(d)))
f2 <- suppressWarnings(frm(bf(y ~ x + (1 | g), muc ~ x + (1 | h)) +
                             categorical(), data = d))
g2 <- frmtmb:::re_governed_b(f2)
b2 <- frmtmb:::re_b_batches(f2, NULL, "y", FALSE, g2)
cat("\ntwo predictors, different grouping factors: kept ", length(g2),
    ", batches ", if (is.null(b2)) "refused" else length(b2), "\n",
    sep = "")
cmp("in sample, batched",
    frmtmb:::fit_fd_se(f2, function(x) frmtmb:::fitted_point(x),
                       b_idx = g2),
    fitted(f2)[, "Est.Error", ])
nd2 <- d[c(1, 61), c("x", "g", "h")]
cmp("2 newdata rows, bounded+batched",
    frmtmb:::fit_fd_se(f2, function(x) frmtmb:::fitted_point(x, nd2),
                       b_idx = g2),
    fitted(f2, newdata = nd2)[, "Est.Error", ])

# the SAME grouping factor on two predictors: one row then loads two
# columns that live in two different blocks, which is the case a
# per-block batch has to keep separate
f3 <- suppressWarnings(frm(bf(y ~ x + (1 | g), muc ~ x + (1 | g)) +
                             categorical(), data = d))
g3 <- frmtmb:::re_governed_b(f3)
b3 <- frmtmb:::re_b_batches(f3, NULL, "y", FALSE, g3)
cat("\nsame factor on two predictors: kept ", length(g3), ", batches ",
    if (is.null(b3)) "refused" else length(b3), "\n", sep = "")
cmp("in sample, batched",
    frmtmb:::fit_fd_se(f3, function(x) frmtmb:::fitted_point(x),
                       b_idx = g3),
    fitted(f3)[, "Est.Error", ])
