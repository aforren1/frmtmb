# Reviewer, lane sampfix, script 01: attack draws_laplace_probe().
#
# For each model, full draws (perturbed ML estimates, lap_pair() as in
# test-laplace-draws.R) and the same draws with every integrated column
# removed, which is the laplace layout. A call either refuses on the
# laplace object, or computes; a computed answer must be identical to
# the full draws' answer, and a call that reads b must refuse.
#
#   Rscript dev/sampfix-rev-01-probe.R      (data seed 77, draws seed 1)
.libPaths(c("C:/Users/adf44/source/r/wt-sampfix-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})

lap_pair <- function(fit, n = 6L, seed = 1L, tweak = NULL) {
  tpl <- fit$frame[["par_template"]]
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  lab <- frmtmb::brms_par_labels(fit)
  stopifnot(length(est) == length(lab))
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, 0.05),
              n, dimnames = list(NULL, lab))
  if (!is.null(tweak)) M <- tweak(M)
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  inner <- setdiff(lab, frmtmb::brms_par_labels(fit, include_random = FALSE))
  outer <- setdiff(colnames(M), inner)
  full <- structure(list(stanfit = NULL, draws = M, fit = fit),
                    class = "frmtmb_draws")
  lap <- full
  lap$draws <- M[, outer, drop = FALSE]
  stopifnot(frmtmb.sample:::draws_is_laplace(lap))
  list(full = full, lap = lap)
}

# expect = "compute" (must equal full) or "refuse"
chk <- function(lab, p, f, expect, seed = 5L) {
  set.seed(seed)
  a <- tryCatch(suppressWarnings(f(p$lap)), error = function(e) e)
  set.seed(seed)
  b <- tryCatch(suppressWarnings(f(p$full)), error = function(e) e)
  got <- if (inherits(a, "error")) {
    if (grepl("integrates them out", conditionMessage(a), fixed = TRUE))
      "refuse" else paste0("ERROR[", substr(conditionMessage(a), 1, 80), "]")
  } else "compute"
  detail <- ""
  if (identical(got, "compute")) {
    nf <- sum(!is.finite(unlist(a)))
    detail <- sprintf("identical-to-full=%s nonfinite=%d", identical(a, b), nf)
    if (!identical(a, b) && !inherits(b, "error")) {
      detail <- paste0(detail, sprintf(" maxabsdiff=%.3g",
                                       max(abs(unlist(a) - unlist(b)),
                                           na.rm = TRUE)))
    }
  }
  full_nf <- if (inherits(b, "error")) "full:ERROR" else
    sprintf("full-nonfinite=%d", sum(!is.finite(unlist(b))))
  verdict <- if (identical(got, expect) &&
                 (expect == "refuse" || grepl("full=TRUE", detail)))
    "PASS" else "FAIL"
  cat(sprintf("%-4s %-44s expect=%-7s got=%-7s %s %s\n", verdict, lab,
              expect, got, detail, full_nf))
}

set.seed(77)
G <- 8; n <- 10
dd <- data.frame(g = factor(rep(seq_len(G), each = n)))
dd$x <- rnorm(nrow(dd))
u <- rnorm(G, 0, 0.6)[dd$g]
dd$y <- 0.5 + 0.4 * dd$x + u + rnorm(nrow(dd), 0, 0.7)
dd$y2 <- -0.3 + 0.2 * dd$x + rnorm(nrow(dd), 0, 0.5)
dd$ys <- 0.3 + 0.4 * dd$x + rnorm(nrow(dd), 0, exp(-0.2 + u))
dd$yb <- rbinom(nrow(dd), 1, plogis(0.2 + 0.8 * dd$x + u))
dd$o <- cut(0.8 * dd$x + u + rlogis(nrow(dd)), c(-Inf, -0.7, 0.6, Inf),
            labels = FALSE)
dd$o <- factor(dd$o, ordered = TRUE)
dd$pos <- exp(0.2 + 0.1 * dd$x + 0.3 * u + rnorm(nrow(dd), 0, 0.2))
nd <- data.frame(x = c(-1, 0, 1), g = factor(c(1, 2, 3), levels = 1:G))
nd0 <- data.frame(x = c(0, 0), g = factor(c(1, 2), levels = 1:G))

q <- function(...) suppressWarnings(suppressMessages(frm(...)))

cat("== (1 | g) gaussian\n")
p <- lap_pair(q(bf(y ~ x + (1 | g)), family = gaussian(), data = dd))
chk("epred default", p, function(d) posterior_epred(d), "refuse")
chk("epred re NA", p, function(d) posterior_epred(d, re_formula = NA), "compute")
chk("linpred re NA", p, function(d) posterior_linpred(d, re_formula = NA), "compute")
chk("predict re NA", p, function(d) posterior_predict(d, re_formula = NA), "compute")
chk("epred newdata re NA", p, function(d) posterior_epred(d, newdata = nd, re_formula = NA), "compute")
chk("epred dpar sigma", p, function(d) posterior_epred(d, dpar = "sigma"), "compute")
chk("fitted re NA", p, function(d) fitted(d, re_formula = NA), "compute")
chk("predict() re NA", p, function(d) predict(d, re_formula = NA), "compute")
chk("residuals re NA", p, function(d) residuals(d, re_formula = NA), "compute")
chk("predictive_error re NA", p, function(d) predictive_error(d, re_formula = NA), "compute")
chk("bayes_R2 re NA", p, function(d) bayes_R2(d, re_formula = NA), "compute")
chk("pp_check re NA (data)", p, function(d) {
  g <- pp_check(d, ndraws = 3, re_formula = NA); g$data }, "compute")
chk("conditional_effects default", p, function(d)
  conditional_effects(d, effects = "x", resolution = 5), "compute")
chk("hypothesis sd", p, function(d)
  hypothesis(d, "sd_g__Intercept > 0.1", class = NULL)$hypothesis, "compute")
chk("epred newdata new level, allow_new_levels", p, function(d)
  posterior_epred(d, newdata = data.frame(x = 1, g = "new"),
                  allow_new_levels = TRUE), "refuse")
chk("epred newdata new level, gaussian", p, function(d)
  posterior_epred(d, newdata = data.frame(x = 1, g = "new"),
                  allow_new_levels = TRUE,
                  sample_new_levels = "gaussian"), "compute")

cat("== (0 + x | g) at x = 0 in newdata\n")
p <- lap_pair(q(bf(y ~ x + (0 + x | g)), family = gaussian(), data = dd))
chk("epred newdata x=0 default re", p, function(d)
  posterior_epred(d, newdata = nd0), "compute")

cat("== random effect in sigma only\n")
p <- lap_pair(q(bf(ys ~ x, sigma ~ 1 + (1 | g)), family = gaussian(), data = dd))
chk("epred default (mu reads no b)", p, function(d) posterior_epred(d), "compute")
chk("linpred default mu", p, function(d) posterior_linpred(d), "compute")
chk("epred dpar sigma", p, function(d) posterior_epred(d, dpar = "sigma"), "refuse")
chk("predict default", p, function(d) posterior_predict(d), "refuse")
chk("predict re NA", p, function(d) posterior_predict(d, re_formula = NA), "compute")
chk("epred newdata default", p, function(d) posterior_epred(d, newdata = nd), "compute")

cat("== multivariate, (1 | g) in y only\n")
p <- lap_pair(q(bf(y ~ x + (1 | g)) + bf(y2 ~ x) + set_rescor(FALSE),
                family = gaussian(), data = dd))
chk("epred resp y default", p, function(d) posterior_epred(d, resp = "y"), "refuse")
chk("epred resp y2 default", p, function(d) posterior_epred(d, resp = "y2"), "compute")
chk("predict resp y2 default", p, function(d) posterior_predict(d, resp = "y2"), "compute")
chk("epred resp y re NA", p, function(d) posterior_epred(d, resp = "y", re_formula = NA), "compute")

cat("== cumulative with (1 | g)\n")
p <- lap_pair(q(bf(o ~ x + (1 | g)), family = cumulative(), data = dd))
chk("epred default", p, function(d) posterior_epred(d), "refuse")
chk("epred re NA", p, function(d) posterior_epred(d, re_formula = NA), "compute")
chk("predict re NA", p, function(d) posterior_predict(d, re_formula = NA), "compute")
chk("hypothesis Intercept[1] < Intercept[2]", p, function(d)
  hypothesis(d, "Intercept[1] < Intercept[2]")$hypothesis, "compute")

cat("== bernoulli with (1 | g)\n")
p <- lap_pair(q(bf(yb ~ x + (1 | g)), family = bernoulli(), data = dd))
chk("epred default", p, function(d) posterior_epred(d), "refuse")
chk("epred re NA", p, function(d) posterior_epred(d, re_formula = NA), "compute")
chk("linpred re NA", p, function(d) posterior_linpred(d, re_formula = NA), "compute")

cat("== lognormal (1 | g) in mu and sigma\n")
p <- lap_pair(q(bf(pos ~ x + (1 | g), sigma ~ (1 | g)), family = lognormal(), data = dd))
chk("epred re NA", p, function(d) posterior_epred(d, re_formula = NA), "compute")
chk("epred re ~(1|g)", p, function(d) posterior_epred(d, re_formula = ~ (1 | g)), "refuse")

cat("== smooth only\n")
p <- lap_pair(q(bf(y ~ s(x, k = 5)), family = gaussian(), data = dd))
chk("epred re NA (smooth read)", p, function(d) posterior_epred(d, re_formula = NA), "refuse")
chk("epred dpar sigma", p, function(d) posterior_epred(d, dpar = "sigma"), "compute")
chk("conditional_effects", p, function(d) conditional_effects(d, effects = "x", resolution = 5), "refuse")

cat("== smooth + (1 | g)\n")
p <- lap_pair(q(bf(y ~ s(x, k = 5) + (1 | g)), family = gaussian(), data = dd))
chk("epred re NA (smooth kept)", p, function(d) posterior_epred(d, re_formula = NA), "refuse")

cat("== probe false negative: NA^0 in a nonlinear body, k = 0 on draw 1\n")
dd$yn <- 1 + 0.3 * dd$x + exp(0.3 + u)^0.8 + rnorm(nrow(dd), 0, 0.3)
fnl <- q(bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1,
            nl = TRUE), family = gaussian(), data = dd)
cat("nl fit converged:", isTRUE(fnl$opt$convergence == 0), "\n")
kcol <- grep("^b_k_Intercept$", brms_par_labels(fnl), value = TRUE)
cat("k column:", kcol, "\n")
p <- lap_pair(fnl)
chk("nl epred default, k free", p, function(d) posterior_epred(d), "refuse")
p0 <- lap_pair(fnl, tweak = function(M) { M[1, kcol] <- 0; M })
chk("nl epred default, k == 0 on draw 1", p0, function(d) posterior_epred(d), "refuse")
set.seed(5)
r0 <- tryCatch(posterior_epred(p0$lap), error = function(e) e)
if (!inherits(r0, "error")) {
  cat(sprintf("  laplace result: draw 1 non-finite %d of %d; draws 2..6 non-finite %d of %d\n",
              sum(!is.finite(r0[1, ])), ncol(r0), sum(!is.finite(r0[-1, ])),
              length(r0[-1, ])))
}
chk("nl epred default, k == 0, draw_ids = 2:6", p0, function(d)
  posterior_epred(d, draw_ids = 2:6), "refuse")
cat("DONE\n")
