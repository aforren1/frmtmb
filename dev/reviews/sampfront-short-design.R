# Short check of lane sampfront's fix round: attack integrated_design_zero()
# and the watch of conditional_effects(method = "posterior_predict").
#   Rscript dev/reviews/sampfront-short-design.R [base]
# Each case: y ~ (0 + x | g) + T at newdata x = 0, so the slope term reads
# nothing and the question is whether T is seen. "SILENT WRONG" is a laplace
# answer that differs from the full draws it comes from (lap_pair, seed 1).
args <- commandArgs(TRUE)
base <- length(args) && identical(args[1], "base")
.libPaths(c(if (!base) "C:/Users/adf44/source/r/wt-sampfront-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb); library(frmtmb.sample); library(frmtmb.spatial)
})
cat("arm", if (base) "base" else "lane", find.package("frmtmb"),
    find.package("frmtmb.sample"), "\n")
lap_pair <- function(fit, n = 6L, seed = 1L) {
  tpl <- fit$frame[["par_template"]]
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  lab <- frmtmb::brms_par_labels(fit)
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, 0.05),
              n, dimnames = list(NULL, lab))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  inner <- setdiff(lab, frmtmb::brms_par_labels(fit, include_random = FALSE))
  full <- structure(list(stanfit = NULL, draws = M, fit = fit),
                    class = "frmtmb_draws")
  lap <- full
  lap$draws <- M[, setdiff(colnames(M), inner), drop = FALSE]
  list(full = full, lap = lap)
}
q <- function(...) suppressWarnings(suppressMessages(frm(...)))
cmp <- function(lab, p, f) {
  cat("\n-- ", lab, "\n", sep = "")
  if (inherits(p, "error")) {
    cat("  fit failed:", substr(conditionMessage(p), 1, 150), "\n")
    return(invisible())
  }
  a <- tryCatch(suppressWarnings(f(p$lap)), error = function(e) e)
  b <- tryCatch(suppressWarnings(f(p$full)), error = function(e) e)
  if (inherits(b, "error")) {
    cat("  full draws failed:", substr(conditionMessage(b), 1, 150), "\n")
    return(invisible())
  }
  if (inherits(a, "error")) {
    cat("  laplace: REFUSED:", substr(gsub("\\s+", " ",
                                           conditionMessage(a)), 1, 90), "\n")
    return(invisible())
  }
  va <- unlist(a); vb <- unlist(b)
  d <- if (length(va) == length(vb)) max(abs(va - vb), na.rm = TRUE) else NA
  cat("  laplace answered: NA cells", sum(is.na(va)), "of", length(va),
      "(full", sum(is.na(vb)), "); max |laplace - full| =", format(d),
      if (is.na(d) || d > 0 || sum(is.na(va)) > sum(is.na(vb)))
        " SILENT WRONG" else " (equal)", "\n")
}
fitp <- function(...) tryCatch(lap_pair(q(...)), error = function(e) e)

set.seed(5)
G <- 6; n <- 60
d <- data.frame(g = factor(rep(1:G, length.out = n)), x = rnorm(n),
                z = runif(n, -2, 2), w = rnorm(n))
d$h <- factor(rep(1:5, each = 12))
d$f <- factor(rep(c("a", "b"), 30))
d$g1 <- factor(sample(1:G, n, TRUE)); d$g2 <- factor(sample(1:G, n, TRUE))
d$loc <- factor(rep(1:10, 6))
u <- rnorm(G); v <- rnorm(10)
d$y <- 0.3 * d$x + u[d$g] + sin(2 * d$z) + v[d$loc] + rnorm(n, 0, 0.4)
d$cnt <- rpois(n, exp(0.2 + 0.3 * u[d$g]))
nd <- data.frame(x = 0, z = c(-1, 0.5, 1.5), w = c(-1, 0, 1),
                 g = factor(c(1, 2, 3), levels = 1:G),
                 h = factor(c(1, 2, 3), levels = 1:5),
                 f = factor(c("a", "b", "a"), levels = c("a", "b")),
                 g1 = factor(c(1, 2, 3), levels = 1:G),
                 g2 = factor(c(4, 5, 6), levels = 1:G),
                 loc = factor(c(1, 2, 3), levels = 1:10))
ep <- function(o) posterior_epred(o, newdata = nd)

cmp("control: (0 + x | g) alone at x = 0",
    fitp(bf(y ~ x + (0 + x | g)), data = d), ep)
cmp("+ s(z)", fitp(bf(y ~ x + (0 + x | g) + s(z)), data = d), ep)
cmp("+ gp(z) exact", fitp(bf(y ~ x + (0 + x | g) + gp(z)), data = d), ep)
cmp("+ gp(z, k = 6)", fitp(bf(y ~ x + (0 + x | g) + gp(z, k = 6)), data = d),
    ep)
cmp("+ s(z, by = w), w varies", fitp(bf(y ~ x + (0 + x | g) + s(z, by = w)),
                                     data = d), ep)
cmp("+ (1 | mm(g1, g2))", fitp(bf(y ~ x + (0 + x | g) + (1 | mm(g1, g2))),
                               data = d), ep)
cmp("+ (1 | gr(h, by = f))",
    fitp(bf(y ~ x + (0 + x | g) + (1 | gr(h, by = f))), data = d), ep)
A <- matrix(0, 10, 10, dimnames = list(1:10, 1:10))
for (i in 1:9) A[i, i + 1] <- A[i + 1, i] <- 1
cmp("+ car(A, gr = loc)",
    fitp(bf(y ~ x + (0 + x | g) + car(A, gr = loc)), data = d,
         data2 = list(A = A)), ep)
cmp("sigma ~ (1 | h): epred of a gaussian does not read sigma",
    fitp(bf(y ~ x + (0 + x | g), sigma ~ (1 | h)), data = d), ep)
cmp("zi ~ (1 | h), zero_inflated_poisson",
    fitp(bf(cnt ~ x + (0 + x | g), zi ~ (1 | h)),
         family = zero_inflated_poisson(), data = d), ep)
cmp("shared id: mu (0 + x | i | g), sigma (1 | i | g)",
    fitp(bf(y ~ x + (0 + x | i | g), sigma ~ (1 | i | g)), data = d), ep)
d$y2 <- d$y + rnorm(n)
cmp("mvbf: y (0 + x | p | g), y2 (1 | p | g), resp = 'y'",
    fitp(mvbf(bf(y ~ x + (0 + x | p | g)), bf(y2 ~ (1 | p | g)),
              rescor = FALSE), data = d),
    function(o) posterior_epred(o, newdata = nd, resp = "y"))
cmp("mvbf: same, both responses",
    fitp(mvbf(bf(y ~ x + (0 + x | p | g)), bf(y2 ~ (1 | p | g)),
              rescor = FALSE), data = d),
    function(o) posterior_epred(o, newdata = nd))
cmp("new levels: (0 + x | g) at x = 0, level 7, allow_new_levels",
    fitp(bf(y ~ x + (0 + x | g)), data = d),
    function(o) posterior_epred(o, newdata = transform(nd, g = factor(7)),
                                allow_new_levels = TRUE))
cmp("new levels: + (1 | h) at a new h, allow_new_levels",
    fitp(bf(y ~ x + (0 + x | g) + (1 | h)), data = d),
    function(o) { set.seed(9); posterior_epred(o, allow_new_levels = TRUE,
      newdata = transform(nd, h = factor(9))) })

# point 2: conditional_effects(method = "posterior_predict") on a model
# whose expected value does not read the integrated values but whose
# predictive draw does (sigma, a shape)
ce <- function(o) {
  set.seed(3)
  r <- conditional_effects(o, effects = "x", re_formula = NULL,
                           method = "posterior_predict",
                           conditions = data.frame(h = factor(2, levels = 1:5)))
  r[[1]][, c("estimate__", "lower__", "upper__")]
}
cmp("CE predict: bf(y ~ x, sigma ~ (1 | h))",
    fitp(bf(y ~ x, sigma ~ (1 | h)), data = d), ce)
cmp("CE predict: negbinomial, shape ~ (1 | h)",
    fitp(bf(cnt ~ x, shape ~ (1 | h)), family = negbinomial(), data = d), ce)
cmp("posterior_predict: bf(y ~ x, sigma ~ (1 | h))",
    fitp(bf(y ~ x, sigma ~ (1 | h)), data = d),
    function(o) { set.seed(3); posterior_predict(o, newdata = nd) })
cat("\nDONE\n")
