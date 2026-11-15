# Item 5408: under quadrature the stationary escape ran on the Laplace
# objective (quad_fit()) and again on the rule (fit_assembled()). How
# often does the second one run, what does it cost, and does it ever
# gain? escape_stationary() is wrapped to log each call that has starts
# (rule calls carry a tally, Laplace ones do not), as
# dev/rel082-quad-escape.R does.
#   Rscript dev/quadseam2-escdup.R <libpath|base>
a <- commandArgs(TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
calls <- new.env()
orig <- get("escape_stationary", ns)
unlockBinding("escape_stationary", ns)
assign("escape_stationary", function(obj, opt, frame, control, bounds,
                                     par_units = NULL, verbose = 0L,
                                     tally = NULL) {
  where <- if (is.null(tally)) "laplace" else "rule"
  n_st <- length(ns$stationary_escapes(obj, opt, frame))
  t0 <- proc.time()[["elapsed"]]
  o <- orig(obj, opt, frame, control, bounds, par_units, verbose, tally)
  el <- proc.time()[["elapsed"]] - t0
  if (n_st) {
    calls$log <- rbind(calls$log, data.frame(
      where = where, starts = o$stationary_escape[["starts"]],
      gain = o$stationary_escape[["gain"]], sec = el))
  }
  o
}, envir = ns)
fam_forced <- custom_family(
  "qesc_norm", dpars = c("mu", "sigma"),
  links = list(mu = "identity", sigma = "log"),
  lpdf = function(y, dpars, aterms) {
    -0.5 * ((y - dpars[["mu"]]) / dpars[["sigma"]])^2 -
      log(dpars[["sigma"]]) - 0.5 * log(2 * pi)
  },
  init_dpars = list(mu = function(y, aterms) mean(y),
                    sigma = function(y, aterms) stats::sd(y)),
  post = list(mean_fn = function(dpars, aterms) dpars[["mu"]],
              stationary = list(sigma = list(at = 0, tol = 1e6,
                                             from = c(1, -1)))))
`%||%` <- function(a, b) if (is.null(a)) b else a
tot <- data.frame()
run <- function(lbl, s, expr) {
  calls$log <- NULL
  t0 <- proc.time()[["elapsed"]]
  f <- suppressWarnings(suppressMessages(expr))
  el <- proc.time()[["elapsed"]] - t0
  lg <- calls$log %||% data.frame(where = character(), starts = numeric(),
                                  gain = numeric(), sec = numeric())
  la <- lg[lg$where == "laplace", , drop = FALSE]
  ru <- lg[lg$where == "rule", , drop = FALSE]
  cat(sprintf(paste0("%s %2d fit %.2f s | laplace %s gain %s | rule %s ",
                     "%.2f s gain %s | ll %.9f\n"), lbl, s, el,
              if (nrow(la)) "ran" else "no",
              if (nrow(la)) format(la$gain[1], digits = 3) else "-",
              if (nrow(ru)) "ran" else "no", sum(ru$sec),
              if (nrow(ru)) format(ru$gain[1], digits = 3) else "-",
              as.numeric(logLik(f))))
  tot <<- rbind(tot, data.frame(lbl = lbl, fit = el, lap = nrow(la),
                                rule = nrow(ru), rsec = sum(ru$sec),
                                rgain = if (nrow(ru)) ru$gain[1] else NA))
}
for (s in 1:20) {
  set.seed(s)
  d <- data.frame(x = rnorm(120), g = factor(rep(1:20, 6)))
  d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.6)[d$g] + rnorm(120)
  run("FORCED", s, frm(bf(y ~ x + (1 | g)), family = fam_forced, data = d,
                       quadrature = TRUE))
}
for (s in 1:40) {
  set.seed(s)
  n <- 200
  d <- data.frame(x = rnorm(n), g = factor(rep(1:20, length.out = n)))
  d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.5)[d$g] + rnorm(n)
  run("SYM", s, frm(bf(y ~ x + (1 | g)), family = skew_normal(), data = d,
                    quadrature = TRUE))
}
for (s in 1:10) {
  set.seed(s)
  n <- 200
  xs <- -abs(rnorm(n)) * 3
  y <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
  dd <- data.frame(y = y, xs = xs, g = factor(rep(1:20, length.out = n)))
  m3 <- mean((y - mean(y))^3) / stats::sd(y)^3
  bad <- list(betad = c(log(stats::sd(y)), 2 * sign(m3) + 0.5 * m3))
  run("STALL", s, frm(bf(y ~ xs + (1 | g), sigma ~ 1, alpha ~ 1),
                      family = skew_normal(), data = dd, quadrature = TRUE,
                      start = bad))
}
for (l in unique(tot$lbl)) {
  t <- tot[tot$lbl == l, ]
  cat(sprintf(paste0("SUM %s fits %d laplace %d rule %d rule secs %.2f of ",
                     "%.2f max rule gain %s\n"), l, nrow(t), sum(t$lap),
              sum(t$rule), sum(t$rsec), sum(t$fit),
              if (any(!is.na(t$rgain))) format(max(t$rgain, na.rm = TRUE),
                                               digits = 3) else "-"))
}
cat("END\n")
