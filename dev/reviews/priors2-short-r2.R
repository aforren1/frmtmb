# Re-check of lane priors2 after fix round 2: open_end_probe(), the held
# NA/NaN warning, and internal_bound() above the range.
# Rscript priors2-short-r2.R <lane|base>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-priors2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
show <- function(lab, expr) {
  ws <- character(0)
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, substr(conditionMessage(w), 1, 100))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) e)
  if (inherits(r, "error")) {
    cat(lab, ": ERROR ", substr(conditionMessage(r), 1, 170), "\n", sep = "")
  } else {
    cat(lab, ": code ", r$opt$convergence, " objective ",
        sprintf("%.15g", r$opt$objective), " par ",
        paste(sprintf("%.10g", r$opt$par), collapse = " "), "  hex ",
        paste(sprintf("%a", c(r$opt$objective, r$opt$par)), collapse = ","),
        "\n", sep = "")
  }
  for (w in ws) cat("    [warning]", w, "\n")
  invisible(r)
}
sp <- set_prior

cat("\n## P1: exponential('identity') on a tiny scale. Objective",
    "(n + 0.5) log b + S / b + b falls for decades as b drops, then turns",
    "up at b* = S / (n + 0.5) (to first order)\n")
set.seed(9)
n <- 5
base_y <- rexp(n)
for (sc in c(1e-4, 1e-9, 1e-12, 1e-30, 1e-100, 1e-300)) {
  e <- data.frame(y = base_y * sc, one = 1)
  S <- sum(e$y)
  bstar <- (-(n + 0.5) + sqrt((n + 0.5)^2 + 4 * S)) / 2
  bstar <- if (is.finite(bstar) && bstar > 0) bstar else S / (n + 0.5)
  cat(sprintf("scale %g: b* = %.10g\n", sc, S / (n + 0.5)))
  show(sprintf("  P1 gamma(0.5,1) lb=0, scale %g", sc),
       frm(y ~ 0 + one, data = e, family = exponential(link = "identity"),
           prior = sp("gamma(0.5, 1)", class = "b", coef = "one", lb = 0)))
  show(sprintf("  P1 gamma(0.5,1) lb=b*/10, scale %g", sc),
       frm(y ~ 0 + one, data = e, family = exponential(link = "identity"),
           prior = sp("gamma(0.5, 1)", class = "b", coef = "one",
                      lb = S / (n + 0.5) / 10)))
}

cat("\n## P2: a pole that the likelihood does not stop: gaussian slope the",
    "data put below 0, gamma(a, 1) with lb = 0, a near 1. The objective",
    "falls by (1 - a) log(b0 / 2^-1074), about (1 - a) 730\n")
set.seed(1)
m <- 120
g <- data.frame(z = rnorm(m), gr = factor(rep(1:10, each = 12)))
g$yneg <- 1 + g$z * 0.3 + rnorm(10)[g$gr] * 0.3 + rnorm(m, sd = 0.3) -
  0.8 * g$z
for (a in c("0.5", "0.9", "0.99", "0.999", "0.9999")) {
  show(sprintf("  P2 gamma(%s, 1) lb=0", a),
       frm(yneg ~ z, data = g,
           prior = sp(sprintf("gamma(%s, 1)", a), class = "b", coef = "z",
                      lb = 0)))
}
show("  P2 beta(0.999, 2) lb=0 ub=1",
     frm(yneg ~ z, data = g, prior = sp("beta(0.999, 2)", class = "b",
                                        coef = "z", lb = 0, ub = 1)))
show("  P2 beta(2, 0.5) lb=0 ub=1 (pole at ub, data below 0)",
     frm(yneg ~ z, data = g, prior = sp("beta(2, 0.5)", class = "b",
                                        coef = "z", lb = 0, ub = 1)))

cat("\n## P3: interior fits beside an open end, warnings and hex\n")
set.seed(1)
nn <- 150
d <- data.frame(x = rnorm(nn), z = rnorm(nn))
d$y <- 1 + 0.4 * d$x - 0.3 * d$z + rnorm(nn, sd = 0.7)
for (pr in c("lognormal(0, 1)", "gamma(2, 2)", "inv_gamma(3, 1)",
             "gamma(0.5, 1)")) {
  show(sprintf("  P3 %s lb=0 on x (inside)", pr),
       frm(y ~ x + z, data = d, prior = sp(pr, class = "b", coef = "x",
                                           lb = 0)))
}
show("  P3 beta(2, 2) lb=0 ub=1 on x (inside)",
     frm(y ~ x + z, data = d, prior = sp("beta(2, 2)", class = "b",
                                         coef = "x", lb = 0, ub = 1)))
set.seed(1)
dc <- data.frame(z = rnorm(120))
dc$y <- 1 + 0.3 * dc$z + rnorm(120, sd = 0.3)
show("  P3 C gamma(0.5,1) lb=0 gaussian (review's)",
     frm(y ~ z, data = dc, prior = sp("gamma(0.5, 1)", class = "b",
                                      coef = "z", lb = 0)))
# poisson("identity") with the slope's start pushed onto 0: the box run
# evaluates log(0) on its way
set.seed(11)
pd <- data.frame(x = runif(80, 0.5, 2))
pd$y <- rpois(80, 1.5 * pd$x)
show("  P3 poisson(identity) lognormal(0,1) lb=0",
     frm(y ~ 0 + x, data = pd, family = poisson(link = "identity"),
         prior = sp("lognormal(0, 1)", class = "b", coef = "x", lb = 0)))
show("  P3 poisson(identity) no prior",
     frm(y ~ 0 + x, data = pd, family = poisson(link = "identity")))
show("  P3 poisson(identity) uniform(0, 5)",
     frm(y ~ 0 + x, data = pd, family = poisson(link = "identity"),
         prior = sp("uniform(0, 5)", class = "b", coef = "x")))

cat("\n## P4: bounds above the inverse link's range\n")
set.seed(1)
dz <- data.frame(z = rnorm(120))
dz$cnt <- rpois(120, 2) * rbinom(120, 1, 0.5)
show("  P4 zi lb = 2", frm(cnt ~ z, data = dz,
                           family = zero_inflated_poisson(),
                           prior = sp("", class = "zi", lb = 2)))
show("  P4 zi lb = -1", frm(cnt ~ z, data = dz,
                            family = zero_inflated_poisson(),
                            prior = sp("", class = "zi", lb = -1)))
show("  P4 zi no prior", frm(cnt ~ z, data = dz,
                             family = zero_inflated_poisson()))
show("  P4 zi lb = 1", frm(cnt ~ z, data = dz,
                           family = zero_inflated_poisson(),
                           prior = sp("", class = "zi", lb = 1)))
ib <- frmtmb:::internal_bound
for (lk in c("logit", "probit", "cloglog", "cauchit", "log")) {
  if (arm != "lane" && lk != "log") next
  pm <- list(scale = "link", link = frmtmb:::get_link(lk))
  for (v in c(2, 1.5, -1, 0)) {
    r <- tryCatch(format(ib(v, pm, "lb"), digits = 17),
                  error = function(e) paste("ERROR", substr(conditionMessage(e),
                                                            1, 60)))
    cat(sprintf("  P4 %-8s lb = %4s -> %s\n", lk, v, r))
  }
}
