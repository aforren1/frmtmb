# Review of lane quadesc: where the rule's gradient turns NaN on a run-off
# restart of the rule's stationary escape (backlog LOW filed by the lane,
# "not isolated"). Base build, ALPHAX at k (default 9).
#  1. The fit, with the rule's escape restarts recorded (fn and gr points).
#  2. At the restart's first point with a non-finite gradient (x*) and at
#     its lowest finite-gradient point (xb): the rule's value and gradient.
#  3. A replica of the model in plain RTMB (b ~ N(0, exp(theta)), as
#     covstruct_core$us), taped with the same GK spec at the rule's own
#     anchor template; checked against the rule at x* and xb. Then one
#     replica per group: each group's rule value, gradient finiteness, and
#     exact value (stats::integrate, split at the cliffs).
#   Rscript dev/reviews/quadesc-review-iso.R [k]
a <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(RTMB)})
k <- if (length(a)) as.numeric(a[1]) else 9
cat("lib:", dirname(find.package("frmtmb")), "version",
    format(packageVersion("frmtmb")), "RTMB", format(packageVersion("RTMB")),
    "TMB", format(packageVersion("TMB")), "\n")
ns <- asNamespace("frmtmb")
st <- new.env()
st$rid <- 0L
put <- function(name, f) {
  unlockBinding(name, ns)
  assign(name, f, envir = ns)
}
o_opt <- ns$optimize_obj
put("optimize_obj", function(obj, control, bounds = list(lower = -Inf,
                                                         upper = Inf),
                             par_units = NULL, verbose = 0L,
                             start_par = obj$par, tally = NULL) {
  if (isTRUE(st$on)) {
    st$rid <- st$rid + 1L
    st$hist[[st$rid]] <- list(f = list(), g = list())
  }
  o_opt(obj, control, bounds, par_units, verbose, start_par, tally)
})
o_esc <- ns$escape_stationary
put("escape_stationary", function(obj, opt, frame, control, bounds,
                                  par_units = NULL, verbose = 0L,
                                  tally = NULL) {
  starts <- ns$stationary_escapes(obj, opt, frame)
  rule <- !is.null(tally) && length(obj$env$random) == 0L
  ob <- obj
  if (rule && length(starts)) {
    st$obj <- obj
    st$hist <- list()
    st$on <- TRUE
    ob$fn <- function(x, ...) {
      v <- obj$fn(x, ...)
      if (st$rid >= 1L) {
        h <- st$hist[[st$rid]]
        h$f[[length(h$f) + 1L]] <- list(x = x, v = v)
        st$hist[[st$rid]] <- h
      }
      v
    }
    ob$gr <- function(x, ...) {
      g <- obj$gr(x, ...)
      if (st$rid >= 1L) {
        h <- st$hist[[st$rid]]
        h$g[[length(h$g) + 1L]] <- list(x = x, g = g)
        st$hist[[st$rid]] <- h
      }
      g
    }
  }
  o <- tryCatch(o_esc(ob, opt, frame, control, bounds, par_units, verbose,
                      tally), finally = st$on <- FALSE)
  o
})

sn <- function(n, alpha) {
  dl <- alpha / sqrt(1 + alpha^2)
  dl * (abs(stats::rnorm(n)) - sqrt(2 / pi)) +
    sqrt(1 - dl^2) * stats::rnorm(n)
}
set.seed(1011)
invisible(stats::rnorm(60 + 15 + 60))
invisible(stats::rnorm(50 + 10))
invisible(sn(50, 1.5))
d <- data.frame(x = stats::rnorm(80), g = factor(rep(1:16, each = 5)))
d$y <- 0.5 * d$x + stats::rnorm(16, 0, 1)[d$g] + sn(80, 2)
d$y <- d$y * (1 + 2e-15 * k)
fo <- bf(y ~ x + (1 | g), sigma ~ 1, alpha ~ x)
f <- suppressWarnings(frm(fo, family = skew_normal(), data = d,
                          quadrature = TRUE,
                          start = list(betad = c(0, 0, 0))))
cat(sprintf("FIT k %g ll %.10f alpha %s\n", k, as.numeric(logLik(f)),
            paste(format(f$opt$par[4:5], digits = 5), collapse = ",")))

h <- st$hist[[1]]
gfin <- vapply(h$g, function(e) all(is.finite(e$g)), NA)
xs <- h$g[[which(!gfin)[1]]]$x
fx <- function(x) {
  m <- Filter(function(q) identical(unname(as.numeric(q$x)),
                                    unname(as.numeric(x))), h$f)
  as.numeric(m[[length(m)]]$v)
}
bestv <- Inf
for (e in h$g[gfin]) {
  v <- fx(e$x)
  if (is.finite(v) && v < bestv) {
    bestv <- v
    xb <- e$x
  }
}
o <- st$obj
tpl <- o$env$parameters
cat("anchor template (the rule's):\n")
str(tpl)
show <- function(lbl, x) {
  v <- o$fn(x)
  g <- o$gr(x)
  cat(sprintf("RULE %s fn %.8f gr %s\n", lbl, v,
              paste(format(g, digits = 4), collapse = " ")))
}
show("xb", xb)
show("x*", xs)
cat("xb:", format(xb, digits = 8), "\nx*:", format(xs, digits = 8), "\n")

ldsn <- function(y, m, s, al) {
  dl <- al / sqrt(1 + al^2)
  om <- s / sqrt(1 - 2 * dl^2 / pi)
  xi <- m - om * dl * sqrt(2 / pi)
  z <- (y - xi) / om
  log(2) - log(om) + stats::dnorm(z, log = TRUE) +
    stats::pnorm(al * z, log.p = TRUE)
}
exact_g <- function(par, dg) {
  b <- par[names(par) == "beta"]
  bd <- par[names(par) == "betad"]
  sdg <- exp(par[names(par) == "theta"][1])
  mu <- b[1] + b[2] * dg$x
  sig <- exp(bd[1])
  al <- bd[2] + bd[3] * dg$x
  L <- 12 * sdg
  lf <- function(u) {
    out <- stats::dnorm(u, 0, sdg, log = TRUE)
    for (j in seq_len(nrow(dg))) out <- out + ldsn(dg$y[j], mu[j] + u, sig,
                                                   al[j])
    out
  }
  dl <- al / sqrt(1 + al^2)
  om <- sig / sqrt(1 - 2 * dl^2 / pi)
  edge <- dg$y - mu + om * dl * sqrt(2 / pi)
  edge <- edge[edge > -L & edge < L]
  br <- sort(unique(c(-L, edge, L)))
  grid <- sort(c(seq(-L, L, length.out = 200001), edge))
  lv <- lf(grid)
  mx <- max(lv[is.finite(lv)])
  s <- 0
  for (j in seq_len(length(br) - 1L)) {
    s <- s + stats::integrate(function(u) exp(lf(u) - mx), br[j], br[j + 1],
                              rel.tol = 1e-11, abs.tol = 0,
                              subdivisions = 2000L)$value
  }
  list(ll = mx + log(s), width = if (length(edge)) diff(range(edge)) else NA,
       support = {
         lo <- edge[al[edge > -L & edge < L] < 0]
         hi <- edge[al[edge > -L & edge < L] > 0]
         c(lo = if (length(lo)) max(lo) else -Inf,
           hi = if (length(hi)) min(hi) else Inf)
       })
}

gk <- ns$quad_gk_integrate()
replica <- function(dg, tplg) {
  dat <- dg
  dat$gi <- as.integer(droplevels(dat$g))
  nll <- function(p) {
    RTMB::getAll(p)
    mu <- beta[1] + beta[2] * dat$x + b[dat$gi]
    al <- betad[2] + betad[3] * dat$x
    -(sum(RTMB::dnorm(b, 0, exp(theta), log = TRUE)) +
        sum(RTMBdist::dskewnorm2(dat$y, mu, exp(betad[1]), al, log = TRUE)))
  }
  RTMB::MakeADFun(nll, tplg, random = "b", integrate = gk, silent = TRUE)
}
cat("\nREPLICA, all groups, at the rule's anchor:\n")
rp <- replica(d, list(beta = unname(tpl$beta), betad = unname(tpl$betad),
                      b = tpl$b, theta = tpl$theta))
for (lbl in c("xb", "x*")) {
  x <- if (lbl == "xb") xb else xs
  cat(sprintf("REP %s fn %.8f gr %s\n", lbl, rp$fn(unname(x)),
              paste(format(rp$gr(unname(x)), digits = 4), collapse = " ")))
}
cat("\nPER GROUP at x* and xb: rule value, gradient finite, exact, diff\n")
for (lbl in c("xb", "x*")) {
  x <- if (lbl == "xb") xb else xs
  tot_r <- 0
  tot_e <- 0
  for (gl in levels(d$g)) {
    ig <- which(d$g == gl)
    dg <- d[ig, ]
    rg <- replica(dg, list(beta = unname(tpl$beta),
                           betad = unname(tpl$betad),
                           b = tpl$b[as.integer(gl)], theta = tpl$theta))
    v <- rg$fn(unname(x))
    g <- rg$gr(unname(x))
    ex <- exact_g(x, dg)
    tot_r <- tot_r + v
    tot_e <- tot_e - ex$ll
    cat(sprintf(paste0("  %s g %2s rule_nll %12.4f gr_finite %-5s",
                       " exact_nll %10.4f rule-exact %10.4f",
                       " support [%.4g, %.4g] width %.3g\n"),
                lbl, gl, v, all(is.finite(g)), -ex$ll, v + ex$ll,
                ex$support[["lo"]], ex$support[["hi"]],
                ex$support[["hi"]] - ex$support[["lo"]]))
  }
  cat(sprintf("  %s sum rule %.6f sum exact %.6f\n", lbl, tot_r, tot_e))
}
cat("END\n")
