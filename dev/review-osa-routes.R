# Review lane osa: the route residuals(type = "osa") takes on cases the
# lane's audit did not fit: mixtures, se(), cens(), a multivariate fit,
# user families written with frmtmb_family(), and the time cost on a
# larger design.
#
#   Rscript dev/review-osa-routes.R before|after
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "after") "C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
cat("BEGIN review-osa-routes", arm, "\n")
seen <- new.env()
suppressMessages(trace(RTMB::oneStepPredict, quote({
  seen$method <- list(...)$method
  seen$discrete <- list(...)$discrete
  seen$support <- !is.null(list(...)$discreteSupport)
}), print = FALSE, where = asNamespace("RTMB")))
run <- function(label, fit, ref = NULL, ...) {
  seen$method <- NA
  seen$discrete <- NA
  seen$support <- NA
  if (inherits(fit, "error")) {
    cat(sprintf("%-38s FIT ERROR %s\n", label, substr(conditionMessage(fit), 1, 120)))
    return(invisible())
  }
  wn <- character(0)
  r <- tryCatch(withCallingHandlers(
    residuals(fit, type = "osa", ...)[, "Estimate"],
    warning = function(w) {
      wn <<- c(wn, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  info <- sprintf("method %-15s discrete %-5s support %-5s", seen$method,
                  seen$discrete, seen$support)
  if (inherits(r, "error")) {
    cat(sprintf("%-38s %s ERROR %s\n", label, info,
                substr(gsub("\n", " ", conditionMessage(r)), 1, 160)))
    return(invisible())
  }
  dif <- if (!is.null(ref)) sprintf(" max|r - ref| %.3g", max(abs(r - ref(fit)), na.rm = TRUE)) else ""
  cat(sprintf("%-38s %s nonfinite %d mean %.3f sd %.3f%s%s\n", label, info,
              sum(!is.finite(r)), mean(r, na.rm = TRUE), sd(r, na.rm = TRUE),
              dif, if (length(wn)) paste(" WARN:", substr(wn[1], 1, 100)) else ""))
}
fitq <- function(...) tryCatch(suppressWarnings(frm(...)), error = function(e) e)
set.seed(501)
n <- 160
d <- data.frame(x = rnorm(n), g = factor(rep(1:16, each = 10)))
d$cnt <- ifelse(runif(n) < 0.4, rpois(n, 1), rpois(n, 8))
d$yc <- ifelse(runif(n) < 0.5, rnorm(n, -2), rnorm(n, 2))
d$se <- runif(n, 0.2, 0.8)
d$yg <- 1 + 0.5 * d$x + rnorm(n, 0, sqrt(1 + d$se^2))
d$cen <- ifelse(d$yg > 2.5, 1, 0)
d$ygc <- pmin(d$yg, 2.5)
d$cnt2 <- rpois(n, exp(0.5 + 0.3 * d$x))

run("mixture(poisson, poisson)",
    fitq(bf(cnt ~ 1), family = mixture(poisson(), poisson()), data = d))
run("mixture(gaussian, gaussian)",
    fitq(bf(yc ~ 1), family = mixture(gaussian(), gaussian()), data = d))
run("gaussian se(se, sigma = TRUE)",
    fitq(yg | se(se, sigma = TRUE) ~ x, family = gaussian(), data = d))
run("student se(se, sigma = TRUE)",
    fitq(yg | se(se, sigma = TRUE) ~ x, family = student(), data = d),
    ref = function(f) {
      dp <- lapply(frmtmb:::eval_dpars(f)[[1]], as.numeric)
      qnorm(pt((d$yg - dp$mu) / sqrt(dp$sigma^2 + d$se^2), dp$nu))
    })
run("gaussian cens(cen) right at 2.5",
    fitq(ygc | cens(cen) ~ x, family = gaussian(), data = d))
run("poisson with offset",
    fitq(cnt2 ~ x + offset(log(1 + 0 * x)), family = poisson(), data = d))
run("mvbf(poisson, gaussian)",
    fitq(mvbf(bf(cnt2 ~ x, family = poisson()),
              bf(yg ~ x, family = gaussian())), data = d))
# a user family with a density only, on a count: no CDF, no support
upois <- tryCatch(frmtmb_family("upois", dpars = "mu",
                                links = list(mu = "log"), type = "discrete",
                                lpdf = function(y, dpars, aterms) {
                                  RTMB::dpois(y, dpars[["mu"]], log = TRUE)
                                },
                                init_dpars = list(mu = function(y, aterms) mean(y) + 0.1)),
                  error = function(e) e)
if (inherits(upois, "error")) {
  cat("user family upois:", conditionMessage(upois), "\n")
} else {
  fu <- fitq(cnt2 ~ x, family = upois, data = d)
  run("user discrete family, lpdf only", fu)
  run("  ... with discreteSupport = 0:40", fu, discreteSupport = 0:40,
      ref = function(f) {
        mu <- exp(f$estimates$beta[[1]] + f$estimates$beta[[2]] * d$x)
        U <- local({ set.seed(123); runif(n) })
        qnorm(ppois(d$cnt2, mu) - U * dpois(d$cnt2, mu))
      })
}
# a user family that supplies a probability-scale lcdf
upois2 <- tryCatch(frmtmb_family("upois2", dpars = "mu",
                                 links = list(mu = "log"), type = "discrete",
                                 lpdf = function(y, dpars, aterms) {
                                   RTMB::dpois(y, dpars[["mu"]], log = TRUE)
                                 },
                                 lcdf = function(q, dpars, aterms) {
                                   RTMB::ppois(q, dpars[["mu"]])
                                 },
                                 init_dpars = list(mu = function(y, aterms) mean(y) + 0.1)),
                   error = function(e) e)
if (inherits(upois2, "error")) {
  cat("user family upois2:", conditionMessage(upois2), "\n")
} else {
  fu2 <- fitq(cnt2 ~ x, family = upois2, data = d)
  run("user discrete family with lcdf", fu2, ref = function(f) {
    mu <- exp(f$estimates$beta[[1]] + f$estimates$beta[[2]] * d$x)
    U <- local({ set.seed(123); runif(n) })
    qnorm(ppois(d$cnt2, mu) - U * dpois(d$cnt2, mu))
  })
}
cat("END review-osa-routes\n")
