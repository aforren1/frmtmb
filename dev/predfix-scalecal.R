# Item 6 calibration, part 2, on the released build: where does the
# default's shortfall begin as n and the effect size change, and does it
# reach REML, profile and a mixed fit, where the mu coefficients are not
# (all) outer parameters?
#   PREDFIX_ARM=base Rscript dev/predfix-scalecal.R > dev/predfix-log/scalecal-base.txt
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
ll <- function(f) if (inherits(f, "try-error")) NA_real_ else
  as.numeric(logLik(f))
fitll <- function(...) ll(try(quiet(frm(...)), silent = TRUE))
options(width = 200)

cat("\n== A. onset scale by n and slope, poisson and gaussian, 0 + xs and 1 + xs\n")
rows <- list()
for (n in c(20, 50, 250, 5000)) for (slope in c(0.05, 0.4)) {
  for (seed in 1:3) {
    set.seed(seed)
    x <- rnorm(n)
    d <- data.frame(x = x, pois = rpois(n, exp(0.5 + slope * x)),
                    gau = 1 + slope * x + rnorm(n))
    for (y in c("pois", "gau")) for (r in c("0 + xs", "1 + xs")) {
      fam <- if (y == "pois") poisson() else gaussian()
      fo <- bf(stats::as.formula(paste(y, "~", r)))
      d$xs <- d$x
      ref <- fitll(fo, family = fam, data = d)
      for (s in 10^-(2:7)) {
        d$xs <- d$x * s
        rows[[length(rows) + 1L]] <- data.frame(
          n = n, slope = slope, seed = seed, y = y, rhs = r, scale = s,
          short_default = ref - fitll(fo, family = fam, data = d),
          short_auto = ref - fitll(fo, family = fam, data = d,
                                   control = frmtmb_control(autoscale = TRUE)))
      }
    }
  }
}
A <- do.call(rbind, rows)
A$def_bad <- A$short_default > 1e-6
A$auto_bad <- A$short_auto > 1e-6
cat("default short (> 1e-6) by n and log10 scale:\n")
print(with(A, tapply(def_bad, list(n = n, log10scale = log10(scale)), sum)))
cat("out of", nrow(A) / (4 * 6), "fits per cell\n")
cat("autoscale short, total:", sum(A$auto_bad, na.rm = TRUE), "of", nrow(A),
    "; NA:", sum(is.na(A$auto_bad)), "\n")
cat("largest scale at which the default is short:",
    max(A$scale[A$def_bad %in% TRUE]), "\n")
saveRDS(A, "C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-log/scalecal-A.rds")

cat("\n== B. REML, profile and (1 | g) at scale 1e-6, against the scale-1 fit\n")
rows <- list()
for (seed in 1:3) {
  set.seed(seed)
  n <- 250
  d <- data.frame(x = rnorm(n), g = factor(sample(1:10, n, TRUE)))
  u <- rnorm(10, 0, 0.5)
  d$pois <- rpois(n, exp(0.5 + 0.4 * d$x + u[d$g]))
  d$gau <- 1 + 0.4 * d$x + u[d$g] + rnorm(n)
  for (y in c("pois", "gau")) for (re in c("", " + (1 | g)")) {
    fam <- if (y == "pois") poisson() else gaussian()
    for (r in c("0 + xs", "1 + xs")) {
      fo <- bf(stats::as.formula(paste(y, "~", r, re)))
      for (mode in c("ML", "REML", "profile")) {
        args <- list(fo, family = fam, data = d)
        if (mode == "REML") args$REML <- TRUE
        ctl <- if (mode == "profile") list(profile = TRUE) else list()
        d$xs <- d$x
        args$data <- d
        ref <- ll(try(quiet(do.call(frm, c(args, list(
          control = do.call(frmtmb_control, ctl))))), silent = TRUE))
        d$xs <- d$x * 1e-6
        args$data <- d
        f0 <- ll(try(quiet(do.call(frm, c(args, list(
          control = do.call(frmtmb_control, ctl))))), silent = TRUE))
        f1 <- ll(try(quiet(do.call(frm, c(args, list(
          control = do.call(frmtmb_control, c(ctl, autoscale = TRUE)))))),
          silent = TRUE))
        rows[[length(rows) + 1L]] <- data.frame(
          seed = seed, y = y, rhs = paste0(r, re), mode = mode,
          short_default = ref - f0, short_auto = ref - f1)
      }
    }
  }
}
B <- do.call(rbind, rows)
print(B, digits = 6)
saveRDS(B, "C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-log/scalecal-B.rds")
