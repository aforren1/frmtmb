# Seeds for the tiny-curvature cliff: does the fit ever stop against
# hscale = exp(rate' dtheta) at a non-stationary point without a warning?
# One row per data seed: the gap is bisected so the first group's likelihood
# curvature at the Laplace optimum is just above zero, then the fit is compared
# with the exact optimum (integrate() on base R densities).
# Usage: Rscript dev/review-impsamp-tiny3.R <from> <to>
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs()
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
args <- as.integer(commandArgs(TRUE))
ff <- bf(y ~ x + (1 | g), nu = 3)
one <- function(ds) {
  set.seed(ds)
  G <- 25L
  g <- rep(seq_len(G), each = 2)
  x <- rnorm(length(g))
  u <- rnorm(G, 0, 0.7)
  base <- 0.3 + 0.5 * x + u[g] + 1.0 * rt(length(g), 6)
  mk <- function(gap) {
    y <- base
    y[1:2] <- c(0.3 + 0.5 * x[1] - gap, 0.3 + 0.5 * x[2] + gap)
    data.frame(y = y, x = x, g = factor(g))
  }
  state <- function(gap) {
    fr <- frm(ff, family = student(), data = mk(gap), dry_run = "frame")
    nll <- ns$build_objective(fr)
    lap <- RTMB::MakeADFun(nll, fr[["par_template"]], random = "b",
                           map = fr[["map"]], silent = TRUE)
    opt <- suppressWarnings(nlminb(lap$par, lap$fn, lap$gr))
    lay <- ns$imp_layout(fr[["re_blocks"]])
    plan <- ns$imp_plan(lap, fr, lay, setNames(opt$par, names(lap$par)), 8L, 1L)
    list(hl = plan$hlik[1, 1, 1], rate = plan$hlik_rate[1, ], opt = opt)
  }
  gaps <- seq(0.5, 5, by = 0.5)
  hl <- vapply(gaps, function(gp) state(gp)$hl, 0)
  if (!any(hl > 0) || !any(hl == 0)) return(NULL)
  lo <- max(gaps[hl > 0])
  hi <- min(gaps[gaps > lo & hl == 0])
  for (i in 1:40) {
    mid <- (lo + hi) / 2
    if (state(mid)$hl > 0) lo <- mid else hi <- mid
    if (hi - lo < 1e-9) break
  }
  s <- state(lo)
  d <- mk(lo)
  w <- character(0)
  fit <- withCallingHandlers(
    tryCatch(frm(ff, family = student(), data = d, importance = 200L),
             error = function(e) e),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  if (inherits(fit, "error")) {
    return(data.frame(ds = ds, gap = lo, hl = s$hl, maxrate = max(abs(s$rate)),
                      rounds = NA, m1 = NA, grad = NA, zmax = NA, zlap = NA,
                      loss = NA, nw = 0,
                      msg = substr(conditionMessage(fit), 1, 150)))
  }
  dl <- split(d, d$g)
  ex <- function(par) {
    sd <- exp(par[4])
    tot <- 0
    for (dg in dl) {
      f <- function(u) {
        exp(vapply(u, function(ui) {
          sg <- exp(par[3])
          sum(dt((dg$y - par[1] - par[2] * dg$x - ui) / sg, 3, log = TRUE) -
                log(sg))
        }, 0) + dnorm(u, 0, sd, log = TRUE))
      }
      br <- seq(-10 * sd - 8, 10 * sd + 8, length.out = 41L)
      v <- 0
      for (i in seq_len(length(br) - 1L)) {
        v <- v + integrate(f, br[i], br[i + 1L], rel.tol = 1e-8, abs.tol = 0,
                           stop.on.error = FALSE)$value
      }
      tot <- tot - log(v)
    }
    tot
  }
  o <- nlminb(unname(fit$opt$par), ex, control = list(iter.max = 40))
  H <- numDeriv::hessian(ex, o$par)
  se <- sqrt(abs(diag(solve(H))))
  data.frame(ds = ds, gap = lo, hl = s$hl, maxrate = max(abs(s$rate)),
             rounds = fit$importance$rounds, m1 = fit$importance$moves[[1]],
             grad = fit$importance$grad,
             zmax = max(abs(unname(fit$opt$par) - o$par) / se),
             zlap = max(abs(s$opt$par - o$par) / se),
             loss = ex(unname(fit$opt$par)) - o$objective, nw = length(w),
             msg = substr(paste(w, collapse = " | "), 1, 150))
}
rows <- list()
for (ds in seq(args[1], args[2])) {
  r <- tryCatch(one(ds), error = function(e) {
    cat("seed", ds, "scan error", conditionMessage(e), "\n")
    NULL
  })
  if (!is.null(r)) {
    rows[[length(rows) + 1L]] <- r
    print(format(r, digits = 4), row.names = FALSE)
  }
}
out <- do.call(rbind, rows)
write.table(out, sprintf(
  "C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-out/tiny3-%d.tsv",
  args[1]), sep = "\t", quote = FALSE, row.names = FALSE)
