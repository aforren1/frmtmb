# Item 1 and 2: a group whose likelihood curvature at the mode is tiny but
# positive. rate = <dH, Hlik> / <Hlik, Hlik> then has a tiny denominator, and
# hscale = exp(rate * dtheta) is a cliff the optimizer meets within a small
# step. Student t: two rows of a group pull the mode apart, and the curvature
# crosses zero as the gap grows. Tune the gap so that it sits just above zero.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs()
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(11)
G <- 15L
g <- rep(seq_len(G), each = 3)
x <- rnorm(length(g))
u <- rnorm(G, 0, 0.7)
base <- 0.3 + 0.5 * x + u[g] + 0.3 * rt(length(g), 4)
mk <- function(gap) {
  y <- base
  y[1:3] <- c(0.3 + 0.5 * x[1] - gap, 0.3 + 0.5 * x[2] + gap, 0.3 + 0.5 * x[3])
  data.frame(y = y, x = x, g = factor(g))
}
ff <- bf(y ~ x + (1 | g), nu = 3)
state <- function(gap) {
  d <- mk(gap)
  fr <- frm(ff, family = student(), data = d, dry_run = "frame")
  nll <- ns$build_objective(fr)
  tpl <- fr[["par_template"]]
  lap <- RTMB::MakeADFun(nll, tpl, random = "b", map = fr[["map"]], silent = TRUE)
  opt <- suppressWarnings(nlminb(lap$par, lap$fn, lap$gr))
  lay <- ns$imp_layout(fr[["re_blocks"]])
  plan <- ns$imp_plan(lap, fr, lay, setNames(opt$par, names(lap$par)), 200L, 1L)
  list(hl = plan$hlik[1, 1, ], rate = plan$hlik_rate, opt = opt, lap = lap)
}
gaps <- seq(0.2, 3, by = 0.2)
tb <- t(vapply(gaps, function(gp) {
  s <- state(gp)
  c(gap = gp, hl1 = s$hl[1], maxrate = max(abs(s$rate[1, ])))
}, numeric(3)))
print(round(tb, 5))
# bisect on the first sign change, in the direction of a tiny positive value
lo <- max(tb[tb[, "hl1"] > 0, "gap"])
hi <- lo + 0.2
for (i in 1:40) {
  mid <- (lo + hi) / 2
  v <- state(mid)$hl[1]
  if (v > 0) lo <- mid else hi <- mid
  if (hi - lo < 1e-9) break
}
cat("gap", lo, "hl1", state(lo)$hl[1], "\n")

exact_for <- function(d) {
  dl <- split(d, d$g)
  gh <- gh_rule(81L)
  function(par) {
    sd <- exp(par[4])
    tot <- 0
    for (dg in dl) {
      fn <- function(u) {
        vapply(u, function(ui) {
          s <- exp(par[3])
          sum(dt((dg$y - par[1] - par[2] * dg$x - ui) / s, 3, log = TRUE) -
                log(s))
        }, 0) + dnorm(u, 0, sd, log = TRUE)
      }
      f <- function(x) exp(fn(x))
      br <- seq(-10 * sd - 8, 10 * sd + 8, length.out = 61L)
      v <- 0
      for (i in seq_len(length(br) - 1L)) {
        v <- v + integrate(f, br[i], br[i + 1L], rel.tol = 1e-8, abs.tol = 0,
                           stop.on.error = FALSE)$value
      }
      tot <- tot - log(v)
    }
    tot
  }
}
for (eps in c(0, 1e-2, 1e-1)) {
  gp <- lo - eps
  s <- state(gp)
  cat("\n==== gap", gp, " hl1", format(s$hl[1], digits = 4),
      " rate[1,]", format(s$rate[1, ], digits = 4), "\n")
  d <- mk(gp)
  w <- character(0)
  fit <- withCallingHandlers(
    tryCatch(frm(ff, family = student(), data = d, importance = 200L,
                 verbose = TRUE),
             error = function(e) e),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  if (inherits(fit, "error")) {
    cat("ERROR:", conditionMessage(fit), "\n")
    next
  }
  cat("warnings:", length(w), substr(paste(w, collapse = " | "), 1, 300), "\n")
  print(fit$importance[c("rounds", "moves", "grads", "capped", "ess_min",
                         "mcse")])
  ex <- exact_for(d)
  o <- nlminb(unname(s$opt$par), ex, control = list(iter.max = 60))
  H <- numDeriv::hessian(ex, o$par)
  se <- sqrt(diag(solve(H)))
  cat("exact optimum:", format(o$par, digits = 5), "nll",
      format(o$objective, digits = 8), "\n")
  cat("fit:", format(unname(fit$opt$par), digits = 5), "nll",
      format(fit$opt$objective, digits = 8), "\n")
  cat("laplace:", format(s$opt$par, digits = 5), "\n")
  cat("z fit:", format((unname(fit$opt$par) - o$par) / se, digits = 3), "\n")
  cat("z lap:", format((s$opt$par - o$par) / se, digits = 3), "\n")
  cat("exact nll at fit - optimum:",
      format(ex(unname(fit$opt$par)) - o$objective, digits = 4), "\n")
}
