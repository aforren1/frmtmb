# Short check of lane priors2's fix round, item 3 (F3), part 2: bounds
# below the floor R's link objects clamp their inverse to (eps through
# logit, probit, cloglog and log), on the user's path, and the rest of
# the unit scan that a refusal stopped in priors2-short-f3.R.
# Rscript priors2-short-f3b.R <lane|base>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-priors2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
ib <- frmtmb:::internal_bound
cat("make.link clamps: logit linkinv(-40) =", make.link("logit")$linkinv(-40),
    " log linkinv(-40) =", make.link("log")$linkinv(-40),
    " plogis(-40) =", plogis(-40), "\n")
set.seed(1)
dz <- data.frame(z = rnorm(120))
dz$cnt <- rpois(120, 2) * rbinom(120, 1, 0.5)
dz$y <- 1 + 0.3 * dz$z + rnorm(120)
dz$cnt0 <- rpois(120, 2)
fit_line <- function(lab, expr) {
  r <- tryCatch(suppressWarnings(suppressMessages(expr)),
                error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-30s ERROR %s\n", lab, substr(conditionMessage(r), 1, 150)))
  } else {
    cat(sprintf("%-30s code %d objective %.15g par %s\n", lab,
                r$opt$convergence, r$opt$objective,
                paste(format(r$opt$par, digits = 10), collapse = " ")))
  }
}
for (u in c(1e-12, 1e-15, 1e-17, 1e-20)) {
  fit_line(paste("zi ub =", u), frm(cnt0 ~ z, data = dz,
                                    family = zero_inflated_poisson(),
                                    prior = set_prior("", class = "zi",
                                                      ub = u)))
  fit_line(paste("zi uniform(0,", u, ")"),
           frm(cnt0 ~ z, data = dz, family = zero_inflated_poisson(),
               prior = set_prior(sprintf("uniform(0, %g)", u), class = "zi")))
}
fit_line("zi lb = 1 - 1e-15", frm(cnt ~ z, data = dz,
                                  family = zero_inflated_poisson(),
                                  prior = set_prior("", class = "zi",
                                                    lb = 1 - 1e-15)))
fit_line("sigma ub = 1e-20", frm(y ~ z, data = dz,
                                 prior = set_prior("", class = "sigma",
                                                   ub = 1e-20)))
fit_line("shape ub = 1e-20", frm(cnt ~ z, data = dz, family = negbinomial(),
                                 prior = set_prior("", class = "shape",
                                                   ub = 1e-20)))
if (arm == "lane") {
  run <- function(name, vs, lab) {
    pm <- list(scale = "link", link = make.link(name))
    out <- 0L
    ref <- 0L
    for (side in c("lb", "ub")) for (v in vs) {
      b <- tryCatch(ib(v, pm, side), error = function(e) NA_real_)
      if (is.na(b)) {
        ref <- ref + 1L
        next
      }
      y <- make.link(name)$linkinv(b)
      if (!isTRUE(if (side == "lb") y >= v else y <= v)) out <- out + 1L
    }
    cat(sprintf("%-8s %-26s ends %5d  outside %d  refused %d\n", name, lab,
                2L * length(vs), out, ref))
  }
  set.seed(3)
  run("logit", runif(5000), "runif(5000)")
  run("probit", runif(5000), "runif(5000)")
  run("probit", 1 - 10^-seq(2, 15, by = 0.01), "1 - 10^-seq(2, 15, .01)")
  run("cloglog", runif(5000), "runif(5000)")
  run("log", exp(runif(5000, -700, 700)), "exp(U(-700, 700))")
  run("log", exp(runif(5000, -30, 30)), "exp(U(-30, 30))")
  run("log", 1 + (-2000:2000) * .Machine$double.eps, "1 + k eps")
  pm_sd <- list(scale = "sd")
  v <- exp(runif(5000, -50, 50))
  bad <- sum(exp(vapply(v, ib, 0, pm = pm_sd, which = "lb")) < v) +
    sum(exp(vapply(v, ib, 0, pm = pm_sd, which = "ub")) > v)
  cat("sd scale, 10000 ends, outside", bad, "\n")
  # termination: the worst case starts at b = 0 with a step of 2^-1074
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(frmtmb:::bound_inside(0, 2, stats::plogis, "lb"),
                error = conditionMessage)
  cat("bound_inside(0, 2, plogis, lb):", r, "\n  elapsed s:",
      proc.time()[["elapsed"]] - t0, "\n")
}
