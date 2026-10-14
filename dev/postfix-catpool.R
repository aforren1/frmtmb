# Lane postfix, fix round 2: pooled category probabilities of a
# frm_multiple() result. Four properties on a cumulative and a
# categorical case, for fitted() and conditional_effects():
# (1) the category estimates sum to 1, (2) each estimate lies inside its
# own interval, (3) fitted() and conditional_effects() agree on the same
# grid, (4) every bound lies in [0, 1]. Also how often the rule had to
# move an interval end to the estimate, and the gap between the
# response-scale estimate and the inverse logit of the mean logit.
.libPaths(c("C:/Users/adf44/source/r/wt-postfix-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb)})
cat("lib:", find.package("frmtmb"), "\n")
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))
props <- function(lbl, fm, nd) {
  fo <- fitted(fm, newdata = nd)
  E <- fo[, "Estimate", , drop = FALSE]
  L <- fo[, "Q2.5", , drop = FALSE]
  U <- fo[, "Q97.5", , drop = FALSE]
  s <- apply(fo[, "Estimate", , drop = FALSE], 1, sum)
  cat(lbl, "fitted: rows", nrow(fo), " categories", dim(fo)[3], "\n")
  cat("  (1) sum of estimates - 1, largest |.|:", max(abs(s - 1)), "\n")
  cat("  (2) estimates outside their interval:", sum(E < L | E > U),
      " ends equal to the estimate:", sum(E == L | E == U), "\n")
  cat("  (4) bounds range:", range(fo[, 3:4, ]), "\n")
  # the logit-scale centre, for the size of the difference
  lg <- lapply(fm$fits, function(f) qlogis(fitted(f, newdata = nd)[, 1, ]))
  ctr <- plogis(Reduce(`+`, lg) / length(lg))
  cat("  response-scale mean vs inverse logit of the mean logit:",
      "largest |diff|", max(abs(fo[, 1, ] - ctr)), "\n")
  invisible(fo)
}
ce_props <- function(lbl, fm) {
  ce <- conditional_effects(fm, "x", categorical = TRUE)[[1]]
  s <- tapply(ce$estimate__, ce$effect1__, sum)
  cat(lbl, "conditional_effects: rows", nrow(ce), "\n")
  cat("  (1) sum of estimates - 1, largest |.|:", max(abs(s - 1)), "\n")
  cat("  (2) estimates outside their band:",
      sum(ce$estimate__ < ce$lower__ | ce$estimate__ > ce$upper__), "\n")
  cat("  (4) bounds range:", range(ce$lower__, ce$upper__), "\n")
  xs <- sort(unique(ce$effect1__))
  fo <- fitted(fm, newdata = data.frame(x = xs))
  i <- match(ce$effect1__, xs)
  lab <- dimnames(fo)[[3]]
  k <- match(as.character(ce$cats__), substr(lab, 7L, nchar(lab) - 1L))
  stopifnot(!anyNA(k))
  cat("  (3) fitted vs conditional_effects on the curve's grid: est",
      rel(fo[cbind(i, 1, k)], ce$estimate__), " lower",
      rel(fo[cbind(i, 3, k)], ce$lower__), " upper",
      rel(fo[cbind(i, 4, k)], ce$upper__), "\n")
}
## cumulative, the review's case (dev/postfix-rev2-rubin.R)
set.seed(1)
dd <- data.frame(x = rnorm(60))
dd$y <- 1 + dd$x + rnorm(60)
imo <- lapply(1:3, function(i) {
  d <- dd; d$o <- cut(d$y + rnorm(60, 0, .3), c(-Inf, 0, 1, 2, Inf),
                      ordered_result = TRUE); d })
fmo <- frm_multiple(bf(o ~ x), family = cumulative(), data = imo)
props("cumulative, x = -3, 0, 3", fmo, data.frame(x = c(-3, 0, 3)))
props("cumulative, in sample", fmo, NULL)
ce_props("cumulative", fmo)
## categorical, 3 categories, 5 imputations of a covariate
set.seed(2)
n <- 150
x <- rnorm(n)
eta2 <- 0.3 + 1.0 * x; eta3 <- -0.4 - 0.8 * x
pr <- cbind(1, exp(eta2), exp(eta3)); pr <- pr / rowSums(pr)
yc <- factor(apply(pr, 1, function(p) sample(c("a", "b", "c"), 1, prob = p)))
imc <- lapply(1:5, function(i) {
  d <- data.frame(y = yc, x = x)
  miss <- seq(1, n, by = 5)
  d$x[miss] <- x[miss] + rnorm(length(miss), 0, 0.5)
  d })
fmc <- frm_multiple(bf(y ~ x), family = categorical(), data = imc)
props("categorical, x = -3, 0, 3", fmc, data.frame(x = c(-3, 0, 3)))
props("categorical, in sample", fmc, NULL)
# conditional_effects() of one categorical fit has no analytic band, and
# the pooled one refuses band = "boot", so property (3) has no curve here
cat("categorical conditional_effects:", tryCatch({
  ce_props("categorical", fmc); "ran"
}, error = conditionMessage), "\n")
