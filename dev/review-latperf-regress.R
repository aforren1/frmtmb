# Review of lane latperf, attack 6: the shared readers on fits WITHOUT a
# latent residual, lane build against rellib-r18, digit for digit. Run both
# arms pinned to one core (dev/review-latperf-regress-run.ps1,
# dev/rtmb-pitfalls.md item 21).
# Usage: REV_ARM=lane|base Rscript dev/review-latperf-regress.R <outfile>
arm <- Sys.getenv("REV_ARM")
base <- c("C:/Users/adf44/source/r/rellib-r18",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lane") c("C:/Users/adf44/source/r/wt-latperf-lib", base)
          else base)
suppressPackageStartupMessages(library(frmtmb))
outf <- commandArgs(TRUE)[[1]]
con <- file(outf, "w")
say <- function(...) writeLines(paste0(...), con)
say("lib: ", dirname(find.package("frmtmb")), " frmtmb ",
    format(packageVersion("frmtmb")))
set.seed(20261106)
d <- expand.grid(week = 1:6, subj = factor(1:30))
d <- d[order(d$subj, d$week), ]
n <- nrow(d)
d$x <- stats::rnorm(n)
d$z <- stats::runif(n)
d$u <- stats::runif(n, 0, 10)
d$f <- factor(sample(c("a", "b", "c"), n, TRUE))
re <- rep(stats::rnorm(30, 0, 0.4), each = 6)
d$y <- 1 + 0.5 * d$x + re + stats::rnorm(n, 0, 0.7)
d$y0 <- 1 + 0.5 * d$x + stats::rnorm(n, 0, 0.7)
d$yp <- 1 + 0.3 * d$u - 0.02 * d$u^2 + re + stats::rnorm(n, 0, 0.5)
d$cnt <- stats::rpois(n, exp(0.7 + 0.3 * d$x + re))
d$zc <- ifelse(stats::runif(n) < 0.25, 0L, d$cnt)
d$nb <- stats::rnbinom(n, mu = exp(0.7 + 0.3 * d$x + re), size = exp(1 + d$z))
d$bin <- stats::rbinom(n, 1, stats::plogis(0.3 * d$x + re))
d$ord <- cut(0.8 * d$x + re + stats::rlogis(n), c(-Inf, -1, 0, 1, Inf),
             ordered_result = TRUE)
d$pos <- stats::rgamma(n, 3, 3 / exp(0.3 * d$x + sin(4 * d$z)))
d$lin <- 1 + 0.4 * d$z + stats::rnorm(n, 0, 0.5)
specs <- list(
  lm = quote(frm(y ~ x * f + z, data = d)),
  lmer = quote(frm(y ~ x + (x | subj), data = d)),
  glmer = quote(frm(cnt ~ x + (1 | subj), family = poisson(), data = d)),
  zip = quote(frm(bf(zc ~ x + (1 | subj), zi ~ z), data = d,
                  family = zero_inflated_poisson())),
  nbshape = quote(frm(bf(nb ~ x + (1 | subj), shape ~ z), data = d,
                      family = negbinomial())),
  bern = quote(frm(bin ~ x + (1 | subj), family = bernoulli(), data = d)),
  ordinal = quote(frm(ord ~ x + (1 | subj), family = cumulative(), data = d)),
  csord = quote(frm(ord ~ cs(x) + (1 | subj), family = acat(), data = d)),
  gamsmooth = quote(frm(pos ~ x + s(z) + (1 | subj), data = d,
                        family = Gamma(link = "log"))),
  # the linear truth runs the smooth's sd to 0: a lost standard error
  lostsd = quote(frm(lin ~ s(z), data = d)),
  # no subject effect in y0: (1 | subj) at its boundary
  boundary = quote(frm(y0 ~ x + (1 | subj), data = d)),
  # an ill-conditioned raw polynomial: the covariance held as a factor
  rawpoly = quote(frm(yp ~ poly(u, 6, raw = TRUE) + (1 | subj), data = d)),
  rawsmooth = quote(frm(yp ~ poly(u, 5, raw = TRUE) + s(z), data = d)),
  mv = quote(frm(bf(y ~ x + (1 | subj)) + bf(cnt ~ x, family = poisson()),
                 data = d)),
  reml = quote(frm(y ~ x + (1 | subj), data = d, REML = TRUE)),
  gp = quote(frm(y ~ x + gp(z, k = 8), data = d)),
  ar1re = quote(frm(cnt ~ x + ar1(factor(week) + 0 | subj),
                    family = poisson(), data = d)),
  gaussar = quote(frm(y ~ x + ar(week, subj, cov = TRUE), data = d))
)
fmt <- function(v) {
  v <- suppressWarnings(as.numeric(unlist(v)))
  paste(sprintf("%.17g", v), collapse = " ")
}
run <- function(label, expr) {
  out <- tryCatch(withCallingHandlers(expr,
    warning = function(w) {
      say("  W ", label, ": ", substr(conditionMessage(w), 1, 100))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage")),
    error = function(e) {
      say("  E ", label, ": ", substr(conditionMessage(e), 1, 140))
      NULL
    })
  if (!is.null(out)) say("  ", label, ": ", out)
}
nd <- d[c(1, 7, 50, 100, 180), ]
nd$subj <- factor(c("1", "2", "new", "9", "new2"))
for (nm in names(specs)) {
  say("### ", nm)
  fit <- tryCatch(suppressWarnings(suppressMessages(eval(specs[[nm]]))),
                  error = function(e) {
                    say("  FIT ERROR ", conditionMessage(e))
                    NULL
                  })
  if (is.null(fit)) next
  sdr <- frmtmb:::sdr_of(fit)
  say("  flags: root ", !is.null(sdr$cov_fixed_root), " prop ",
      !is.null(sdr$cov_fixed_prop), " lost ", length(sdr$se_lost))
  run("fitted", fmt(fitted(fit)))
  run("fitted_nd", fmt(fitted(fit, newdata = nd, allow_new_levels = TRUE)))
  run("linpred", fmt(frm_linpred(fit, se.fit = TRUE)))
  run("residuals", fmt(residuals(fit)))
  run("predict", {
    set.seed(5)
    fmt(predict(fit, newdata = nd[1:3, ], allow_new_levels = TRUE))
  })
  run("confint", fmt(confint(fit)))
  run("ranef", fmt(lapply(ranef(fit), function(x) unlist(x))))
  run("summary", paste(utils::capture.output(print(summary(fit))),
                       collapse = " | "))
  run("jointcov", fmt(frm_joint_cov(fit)$V))
  run("lpbasis", {
    lb <- frm_lp_basis(fit, newdata = nd, allow_new_levels = TRUE)
    fmt(list(as.matrix(lb$A), lb$V, lb$coef_pos))
  })
  if (nm %in% c("gamsmooth", "lostsd", "rawsmooth")) {
    run("csmooths", fmt(lapply(conditional_smooths(fit), function(x) {
      x[, c("estimate__", "se__")]
    })))
  }
  if (requireNamespace("emmeans", quietly = TRUE) &&
        nm %in% c("lm", "glmer", "zip", "rawpoly", "lostsd")) {
    run("emmeans", fmt(as.data.frame(emmeans::emmeans(
      fit, if (nm == "lm") ~ f else ~ 1))[, c("emmean", "SE")]))
  }
}
close(con)
