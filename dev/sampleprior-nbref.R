# Lane sampleprior, fix round 1, item 8: the negbinomial agreement case
# against a deterministic reference, and brms 2.23.0's own chains on it
# over several seeds.
#
# The reference is brms's posterior of log(shape) by quadrature: brms's
# model (centered Intercept with its default student_t, flat slope,
# shape ~ inv_gamma(0.4, 0.3)) on a grid over (Intercept, slope,
# log(shape)), with the Jacobian of log(shape). No sampler is involved.
#   Rscript dev/sampleprior-nbref.R <lane|r9> [brms]
a <- commandArgs(TRUE)
arm <- a[1]
lib <- if (identical(arm, "r9")) character(0) else
  "C:/Users/adf44/source/r/wt-sampleprior-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sampleprior"
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
source(file.path(wt, "dev/sampleprior-agree-common.R"))
cat("arm", arm, "lib:", find.package("frmtmb.sample"), "\n")
cs <- agree_cases()$negbinomial
d <- cs$data

nb_reference <- function(d, loc, scl) {
  xc <- d$x - mean(d$x)
  fit <- suppressWarnings(stats::glm(y ~ xc, data = d,
                                     family = stats::poisson()))
  se <- sqrt(diag(stats::vcov(fit)))
  mu <- mean(d$y)
  ug <- seq(-8, 40, by = 0.05)
  # per log(shape), a grid over (Intercept, slope) as wide as that
  # shape's spread: the negative binomial's variance is mu + mu^2 / shape,
  # so its standard errors are the Poisson ones times sqrt(1 + mu / shape)
  lpu <- vapply(ug, function(u) {
    k <- sqrt(1 + mu / exp(u))
    ag <- seq(stats::coef(fit)[1] - 12 * k * se[1],
              stats::coef(fit)[1] + 12 * k * se[1], length.out = 49)
    bg <- seq(stats::coef(fit)[2] - 12 * k * se[2],
              stats::coef(fit)[2] + 12 * k * se[2], length.out = 49)
    ab <- expand.grid(a = ag, b = bg)
    eta <- outer(ab$a, rep(1, nrow(d))) + outer(ab$b, xc)
    ll <- rowSums(matrix(stats::dnbinom(rep(d$y, each = nrow(ab)),
                                        size = exp(u),
                                        mu = exp(as.vector(eta)),
                                        log = TRUE), nrow(ab)))
    lp <- ll + stats::dt((ab$a - loc) / scl, 3, log = TRUE)
    m <- max(lp)
    # the cell area changes with u, and belongs in the integral
    m + log(sum(exp(lp - m))) + log(diff(ag)[1]) + log(diff(bg)[1])
  }, 1)
  # inv_gamma(0.4, 0.3) on shape and the Jacobian of u = log(shape)
  lpu <- lpu + 0.4 * log(0.3) - lgamma(0.4) - 1.4 * ug - 0.3 / exp(ug) + ug
  w <- exp(lpu - max(lpu))
  w <- w / sum(w)
  cdf <- cumsum(w)
  q <- function(p) stats::approx(cdf, ug, xout = p, ties = "ordered")$y
  list(mean = sum(w * ug), sd = sqrt(sum(w * ug^2) - sum(w * ug)^2),
       q = c(q(0.05), q(0.5), q(0.95)),
       edge = w[c(1L, length(w))])
}
# brms's own default on the centered intercept, read off brms
pri <- with(as.data.frame(brms::default_prior(
  y ~ x, data = d, family = brms::negbinomial())), prior[class == "Intercept"])
nums <- as.numeric(strsplit(sub("^student_t[(](.*)[)]$", "\\1", pri),
                            ",")[[1]])
t0 <- proc.time()[["elapsed"]]
ref <- nb_reference(d, nums[2], nums[3])
cat(sprintf(paste0("reference (quadrature; Intercept %s): mean %.4f ",
                   "sd %.4f q05/50/95 %s, grid-edge weights %.2g %.2g, ",
                   "%.1f s\n"), pri, ref$mean, ref$sd,
            paste(sprintf("%.3f", ref$q), collapse = " "), ref$edge[1],
            ref$edge[2], proc.time()[["elapsed"]] - t0))
for (seed in c(3L, 101L, 202L)) {
  lf <- log(agree_frm(cs, seed, 2L, 2000L))
  mf <- agree_mcse(lf, 2)
  w90 <- function(v) diff(stats::quantile(v, c(0.05, 0.95)))
  cat(sprintf(paste0("frmtmb seed %3d: mean %.4f mcse %.4f rhat %.3f ",
                     "ess %.0f q05/50/95 %s | z vs reference %.2f, 90%% ",
                     "width ratio %.3f\n"), seed, mean(lf), mf,
              posterior::rhat(matrix(lf, ncol = 2)),
              posterior::ess_bulk(matrix(lf, ncol = 2)),
              paste(sprintf("%.3f", stats::quantile(lf, c(.05, .5, .95))),
                    collapse = " "), (mean(lf) - ref$mean) / mf,
              w90(lf) / diff(ref$q[c(1, 3)])))
}
if (identical(a[2], "brms")) {
  for (sp in list(c(3L, 2L, 2000L), c(101L, 2L, 2000L), c(202L, 2L, 2000L),
                  c(101L, 4L, 6000L))) {
    lb <- log(suppressWarnings(agree_brms(cs, sp[1], sp[2], sp[3],
                                          file.path(wt, "dev/stan-cache"))))
    m <- matrix(lb, ncol = sp[2])
    cat(sprintf(paste0("brms seed %3d, %d x %d: mean %.4f rhat %.3f ess %.0f ",
                       "q05/50/95 %s\n"), sp[1], sp[2], sp[3], mean(lb),
                posterior::rhat(m), posterior::ess_bulk(m),
                paste(sprintf("%.3f", stats::quantile(lb, c(.05, .5, .95))),
                      collapse = " ")))
  }
}
