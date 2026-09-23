# Lane wt-reunc, punch round 1 (M1): a FACTOR-SMOOTH block on the
# finite-difference route.
#
# re_governed_b() keeps a smooth whose basis is indexed by a grouping
# factor, so s(x, g, bs = "fs") joins the differenced vector on an
# ordinal or categorical fit. re_row_support(), which the batch reads
# to attribute a row's difference to a level, builds its indicator from
# each predictor's Z in sample and from re_design_matrix() on newdata.
#
# The newdata branch USED NOT TO carry the smooth parts, and a batch
# whose rows all say "no owner" writes a derivative of zero for the
# whole block instead of refusing, so the smooth's columns left the
# variance. That is FIXED: re_row_support() adds the smooth parts on
# newdata, and dev/reunc-log/fdsmooth-prefix.txt is the run that had
# the defect. This script is the guard that it stays fixed, and it
# checks the answer against the one-at-a-time reference, which is the
# only thing that settles it.
#
#   Rscript dev/reunc-fdsmooth.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(lib)) lib <- "C:/Users/adf44/source/r/reunc2-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")

cmp <- function(tag, a, b) {
  a <- unname(as.vector(a))
  b <- unname(as.vector(b))
  cat(sprintf("  %-28s identical %-5s max rel %.3g  max ulp %.1f\n", tag,
              identical(a, b), max(abs(a - b) / abs(a)),
              max(abs(a - b) / (.Machine$double.eps * abs(a)))))
}

set.seed(21)
ng <- 6
d <- data.frame(g = factor(rep(seq_len(ng), each = 25)),
                x = runif(ng * 25, -2, 2))
amp <- rnorm(ng, 0, 1)
lat <- amp[d$g] * sin(d$x) + 0.5 * d$x + rlogis(nrow(d))
d$y <- factor(cut(lat, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
              ordered = TRUE)
fit <- suppressWarnings(frm(bf(y ~ s(x, g, bs = "fs", k = 5)) +
                              cumulative(), data = d))
nd <- d[c(3, 40, 90), c("x", "g")]

gov <- frmtmb:::re_governed_b(fit)
used_in <- frmtmb:::re_used_b(fit, NULL, "y", FALSE)
used_nd <- frmtmb:::re_used_b(fit, nd, "y", FALSE)
bt_in <- frmtmb:::re_b_batches(fit, NULL, "y", FALSE, gov)
bt_nd <- frmtmb:::re_b_batches(fit, nd, "y", FALSE, used_nd)
cat("\nfactor smooth, ", ng, " levels: kept ", length(gov),
    ", loaded in sample ", length(used_in), ", by 3 newdata rows ",
    length(used_nd), "\n", sep = "")
cat("  batches in sample: ",
    if (is.null(bt_in)) "refused" else length(bt_in),
    "; on newdata: ", if (is.null(bt_nd)) "refused" else length(bt_nd),
    "\n", sep = "")

cmp("in sample",
    frmtmb:::fit_fd_se(fit, function(x) frmtmb:::fitted_point(x),
                       b_idx = gov),
    fitted(fit)[, "Est.Error", ])
cmp("3 newdata rows",
    frmtmb:::fit_fd_se(fit, function(x) frmtmb:::fitted_point(x, nd),
                       b_idx = gov),
    fitted(fit, newdata = nd)[, "Est.Error", ])

# a factor smooth BESIDE an ordinary group term, which is the shape
# where one block can batch and the other cannot
set.seed(22)
d$h <- factor(rep(1:5, length.out = nrow(d)))
f2 <- suppressWarnings(frm(bf(y ~ s(x, g, bs = "fs", k = 5) + (1 | h)) +
                             cumulative(), data = d))
g2 <- frmtmb:::re_governed_b(f2)
nd2 <- d[c(3, 40, 90), c("x", "g", "h")]
u2 <- frmtmb:::re_used_b(f2, nd2, "y", FALSE)
cat("\nsmooth plus (1 | h): kept ", length(g2), ", loaded by 3 rows ",
    length(u2), "\n", sep = "")
cmp("in sample",
    frmtmb:::fit_fd_se(f2, function(x) frmtmb:::fitted_point(x),
                       b_idx = g2),
    fitted(f2)[, "Est.Error", ])
cmp("3 newdata rows",
    frmtmb:::fit_fd_se(f2, function(x) frmtmb:::fitted_point(x, nd2),
                       b_idx = g2),
    fitted(f2, newdata = nd2)[, "Est.Error", ])
