# Review of lane famend, item 3560: prior_start_inside() on base and lane.
#   Rscript dev/reviews/famend-review-prior.R <base|lane>
a <- commandArgs(trailingOnly = TRUE)
lanelib <- "C:/Users/adf44/source/r/wt-famend-lib"
base <- "C:/Users/adf44/source/r/rellib-r21"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (a[1] == "lane") lanelib, base, ulib))
suppressPackageStartupMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "\n")

cap <- function(expr) {
  w <- character(); m <- character()
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    },
    message = function(x) {
      m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage")
    }),
    error = function(e) structure(conditionMessage(e), class = "err"))
  list(v = v, w = w, m = m)
}
show <- function(label, r) {
  cat("\n==", label, "\n")
  for (x in r$w) cat("   W:", substr(gsub("\\s+", " ", x), 1, 200), "\n")
  for (x in r$m) cat("   M:", substr(gsub("\\s+", " ", x), 1, 200), "\n")
  f <- r$v
  if (inherits(f, "err")) {
    cat("   ERR:", substr(f, 1, 200), "\n")
    return(invisible(NULL))
  }
  st <- f$obj$env$par[names(f$obj$env$par) != "u"]
  cat("   start", paste(signif(f$start$beta %||% st, 7), collapse = " "), "\n")
  cat("   code", f$opt$convergence, " objective",
      sprintf("%.10f", f$opt$objective), " par",
      paste(sprintf("%.7g", f$opt$par), collapse = " "), "\n")
}
`%||%` <- function(a, b) if (is.null(a)) b else a

set.seed(77)
n <- 120
d <- data.frame(x = rnorm(n), fac = factor(rep(c("a", "b"), n / 2)))
d$ypos <- 1 + 0.3 * d$x + rnorm(n)
d$yneg <- 1 - 0.5 * d$x + rnorm(n)
d$ymid <- 1 + 0.75 * d$x + rnorm(n, 0, 0.5)
P <- function(...) set_prior(..., class = "b", coef = "x")

cases <- list(
  "gamma(2,1) slope>0 no bound (guard absent)" =
    list(f = ypos ~ x, p = P("gamma(2, 1)")),
  "gamma(2,1) slope<0 no bound" =
    list(f = yneg ~ x, p = P("gamma(2, 1)")),
  "gamma(2,1) slope<0 lb=0.1 ub=1 (center 2 outside box)" =
    list(f = yneg ~ x, p = P("gamma(2, 1)", lb = 0.1, ub = 1)),
  "gamma(2,1) slope>0 lb=0.1 ub=1 (guard absent)" =
    list(f = ypos ~ x, p = P("gamma(2, 1)", lb = 0.1, ub = 1)),
  "beta(2,2) slope 0.75 lb=0.6 ub=0.9 (guard absent)" =
    list(f = ymid ~ x, p = P("beta(2, 2)", lb = 0.6, ub = 0.9)),
  "beta(2,2) slope<0 lb=0.6 ub=0.9 (center 0.5 outside box)" =
    list(f = yneg ~ x, p = P("beta(2, 2)", lb = 0.6, ub = 0.9)),
  "uniform(2,3) slope 0.3" =
    list(f = ypos ~ x, p = P("uniform(2, 3)")),
  "exponential(1) slope<0" =
    list(f = yneg ~ x, p = P("exponential(1)")),
  "lognormal(0,1) slope<0 lb=1e-8" =
    list(f = yneg ~ x, p = P("lognormal(0, 1)", lb = 1e-8)),
  "student_t(3,0,2.5) slope<0 (guard absent)" =
    list(f = yneg ~ x, p = P("student_t(3, 0, 2.5)")),
  "normal(0,1) b lb=0 slope<0 (guard absent: half-normal)" =
    list(f = yneg ~ x, p = P("normal(0, 1)", lb = 0)),
  "inv_gamma(2,1) on sigma (guard absent)" =
    list(f = ypos ~ x, p = set_prior("inv_gamma(2, 1)", class = "sigma")),
  "gamma(2,1) on Intercept, y mean -3" =
    list(f = yi ~ x, p = set_prior("gamma(2, 1)", class = "Intercept")),
  "gamma(2,1) slope<0, caller start beta = c(1, -0.5)" =
    list(f = yneg ~ x, p = P("gamma(2, 1)"), st = list(beta = c(1, -0.5))),
  "gamma(2,1) sigma ~ x slope b of sigma, true slope -0.5" =
    list(f = bf(ys ~ 1, sigma ~ x), p = set_prior("gamma(2, 1)", class = "b",
                                               dpar = "sigma"))
)
d$yi <- -3 + 0.2 * d$x + rnorm(n)
d$ys <- rnorm(n, 0, exp(-0.5 * d$x))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  show(nm, cap(frm(cs$f, data = d, prior = cs$p, start = cs$st)))
}

# a nonlinear model with a lognormal prior on a rate that the default
# start puts at 0
set.seed(5)
dn <- data.frame(x = runif(100, 0, 5))
dn$y <- 3 * exp(-0.7 * dn$x) + rnorm(100, 0, 0.1)
pr <- c(set_prior("normal(2, 2)", nlpar = "a"),
        set_prior("lognormal(0, 1)", nlpar = "b"))
show("nl a * exp(-b x), lognormal on b",
     cap(frm(bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE), data = dn,
             prior = pr)))
pr2 <- c(set_prior("normal(2, 2)", nlpar = "a"),
         set_prior("gamma(2, 2)", nlpar = "b"))
show("nl a * exp(-b x), gamma(2,2) on b",
     cap(frm(bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE), data = dn,
             prior = pr2)))
show("nl, gamma(2,2) on b, caller start b = -1",
     cap(frm(bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE), data = dn,
             prior = pr2, start = list(b = -1))))
cat("\ndone\n")
