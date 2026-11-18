# rel084 release review, slopetrig fix round 1: autoscale's ratio is read
# over the rows a block covers (any of its Z columns in that predictor,
# its intercept included, nonzero). New designs: gr(by =) with three
# unequal levels, (0 + dose | g) with dose 0 on 60 percent of rows,
# (1 + dose | g) on the same column, and one fit of the first design on
# both builds against the centered fit.
#   REVLIB=r23|r22 Rscript dev/reviews/rel084-review-slopetrig.R
source("dev/reviews/rel084-review-common.R")
ns <- asNamespace("frmtmb")
ratio_of <- function(x) sqrt(mean(x^2)) / stats::sd(x)
plan_ratios <- function(fo, d) {
  fr <- frm(fo, data = d, dry_run = "frame")
  pl <- tryCatch(ns$autoscale_plan(fr), error = function(e) NULL)
  zz <- pl$y.mu$z %||% pl[[1]]$z
  got <- vapply(zz, function(e) e$ratio %||% NA_real_, 0)
  dec <- !is.null(ns$autoscale_decide(fr, frmtmb_control(), FALSE))
  list(got = got, dec = dec)
}
cat("threshold autoscale_slope_ratio:", ns$autoscale_slope_ratio, "\n")
set.seed(8441)
G <- 24
g <- factor(rep(seq_len(G), times = rep(c(8, 12, 20), each = 8)))
n <- length(g)
h <- factor(c("a", "b", "c")[ (as.integer(g) - 1) %/% 8 + 1 ])
z <- rnorm(n)
d <- data.frame(g = g, h = h, x = 1990 + 2 * z + as.integer(h))
d$y <- 1 + rnorm(G, 0, 0.6)[g] + (0.4 + rnorm(G, 0, 0.4)[g]) * z + rnorm(n)
want <- vapply(levels(h), function(l) ratio_of(d$x[d$h == l]), 0)
allr <- ratio_of(d$x)
r <- plan_ratios(y ~ 1 + (1 + x | gr(g, by = h)), d)
cat(sprintf("A gr(by = h) 3 levels: plan %s | covered %s | every row %.3f | engage %s\n",
            paste(sprintf("%.2f", sort(r$got)), collapse = " "),
            paste(sprintf("%.2f", sort(want)), collapse = " "), allr, r$dec))
cat(sprintf("  max rel diff plan vs covered: %.2e\n",
            max(abs(sort(r$got) / sort(want) - 1))))
dose <- ifelse(runif(n) < 0.6, 0, 5 + rnorm(n))
d$dose <- dose
r <- plan_ratios(y ~ 1 + (0 + dose | g), d)
cat(sprintf("B (0 + dose | g): plan %s | nonzero rows %.3f | every row %.3f | engage %s\n",
            paste(sprintf("%.3f", r$got), collapse = " "),
            ratio_of(dose[dose != 0]), ratio_of(dose), r$dec))
r <- plan_ratios(y ~ 1 + (1 + dose | g), d)
cat(sprintf("C (1 + dose | g): plan %s | every row %.3f | engage %s\n",
            paste(sprintf("%.3f", r$got), collapse = " "), ratio_of(dose), r$dec))
# one fit, design A: default against the centered fit
dc <- d
dc$xc <- (d$x - mean(d$x)) / sd(d$x)
fd <- conds_of(frm(y ~ 1 + (1 + x | gr(g, by = h)), data = d))
fc <- conds_of(frm(y ~ 1 + (1 + xc | gr(g, by = h)), data = dc))
cat(sprintf("A fit default ll %.6f [%s] | centered ll %.6f [%s] | short %.4f\n",
            as.numeric(logLik(fd$v)), cls_line(fd), as.numeric(logLik(fc$v)),
            cls_line(fc), as.numeric(logLik(fc$v)) - as.numeric(logLik(fd$v))))
cat("DONE\n")
