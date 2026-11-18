# Lane slopetrig, item 5058: a correlation run to its limit at a short
# stop is called a boundary. The design of dev/slopescale-verdict.R:
# gaussian and poisson, y ~ 1 + (1 + x | g), x spread 1e-4 or 1e-5,
# 20 groups of 15, seeds 1 to 6, autoscale = FALSE. Per fit: the
# correlation parameter t (rho = t / sqrt(1 + t^2)), the shortfall
# against the same model fitted on the standardized column, the
# conditions, and the change of the objective (negative log-likelihood)
# when t alone is moved back to the clamp, to 1 and to 0 (a negative
# change is a gain). Then the counts.
# Usage: Rscript slopetrig-cor.R base|lane
arm <- commandArgs(TRUE)[1]
stopifnot(arm %in% c("base", "lane"))
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-slopetrig-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb", as.character(packageVersion("frmtmb")),
    "lib:", dirname(find.package("frmtmb")), "\n")
run <- function(expr) {
  txt <- character(0)
  cl <- character(0)
  f <- withCallingHandlers(expr,
    warning = function(w) {
      cl <<- c(cl, "W")
      txt <<- c(txt, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      cl <<- c(cl, "M")
      txt <<- c(txt, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  list(fit = f, cl = cl, txt = txt)
}
tc <- 0.9999 / sqrt(1 - 0.9999^2)
rows <- list()
for (fam in c("gau", "pois")) for (s in c(1e-4, 1e-5)) for (seed in 1:6) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  x0 <- rnorm(300)
  eta <- 0.5 + rnorm(20, 0, 0.5)[g] + rnorm(20, 0, 0.5)[g] * x0
  y <- if (fam == "gau") eta + rnorm(300) else rpois(300, exp(eta))
  d <- data.frame(g = g, x = x0 * s, y = y)
  fm <- if (fam == "gau") gaussian() else poisson()
  fo <- y ~ 1 + (1 + x | g)
  r <- run(frm(fo, family = fm, data = d,
               control = frmtmb_control(autoscale = FALSE)))
  dc <- d
  dc$x <- (d$x - mean(d$x)) / sd(d$x)
  ref <- run(frm(fo, family = fm, data = dc))
  f <- r$fit
  p <- f$opt$par
  jt <- grep("theta", names(p))
  j <- jt[3]
  t0 <- p[j]
  f0 <- f$obj$fn(p)
  dq <- function(v) {
    q <- p
    q[j] <- v
    f$obj$fn(q) - f0
  }
  pr <- c(clamp = dq(sign(t0) * tc), one = dq(sign(t0)), zero = dq(0))
  f$obj$fn(p)
  bnd <- any(grepl("at the limit of its range", r$txt, fixed = TRUE))
  rows[[length(rows) + 1L]] <- data.frame(
    fam = fam, s = s, seed = seed, conv = f$opt$convergence,
    t = signif(t0, 4),
    short = signif(as.numeric(logLik(ref$fit)) - as.numeric(logLik(f)), 3),
    n_warn = sum(r$cl == "W"), boundary_msg = bnd,
    d_clamp = signif(pr[["clamp"]], 3), d_one = signif(pr[["one"]], 3),
    d_zero = signif(pr[["zero"]], 3))
  if (bnd || abs(t0) > tc) {
    cat("\n==", fam, s, seed, "\n")
    cat(paste0("  ", substr(r$txt, 1, 400)), sep = "\n")
  }
}
r <- do.call(rbind, rows)
options(width = 200)
print(r, row.names = FALSE)
cat(sprintf(
  "\nTOTAL fits %d ; t past clamp %d ; boundary message %d ; short>1e-3 %d ; short and boundary message %d ; short and silent %d ; back move gains>1e-6 among past clamp %d\n",
  nrow(r), sum(abs(r$t) > tc), sum(r$boundary_msg), sum(r$short > 1e-3),
  sum(r$short > 1e-3 & r$boundary_msg),
  sum(r$short > 1e-3 & r$n_warn == 0),
  sum(abs(r$t) > tc & pmin(r$d_clamp, r$d_one, r$d_zero) < -1e-6)))
