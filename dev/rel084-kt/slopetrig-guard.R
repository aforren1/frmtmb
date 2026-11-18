# rel084-kt copy of dev/slopetrig-guard.R: lib = rellib-r23 then user lib (lane lib and rellib-r22 removed). No other change.
# Lane slopetrig, item 5058: the guard's condition absent and present,
# over roundings. Absent: a true boundary correlation (no slope variance
# in the data, test-review-v29.R's "a lost correlation leaves the sd
# intervals of its block", seeds 2, 7, 8), which must keep its boundary
# message, and the dimension-3 fixture of the same file. Present: the
# 5058 fits (dev/slopetrig-cor.R; gaussian 1e-4 seeds 4 and 5, poisson
# 1e-5 seed 5), which must warn "stopped short" and give no boundary
# message. Each over y * (1 + 2e-15 k), k = 0..20.
# Prints per fit: the correlation parameter, se_cor_gain_back() of it,
# the lost reasons, whether the boundary message and a "stopped short"
# warning came. Then the counts per case.
# Usage: Rscript slopetrig-guard.R base|lane
arm <- commandArgs(TRUE)[1]
stopifnot(arm %in% c("base", "lane"))
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/rellib-r23",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb", as.character(packageVersion("frmtmb")),
    "lib:", dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
run <- function(expr) {
  w <- character(0)
  m <- character(0)
  f <- withCallingHandlers(expr,
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    },
    message = function(x) {
      m <<- c(m, paste(class(x)[1], conditionMessage(x)))
      invokeRestart("muffleMessage")
    })
  list(fit = f, w = w, m = m)
}
nobs_slope <- function(s) {
  set.seed(s)
  g <- factor(rep(seq_len(20), each = 10))
  d <- data.frame(g = g, x = stats::rnorm(200))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(20, 0, 0.7)[g] + stats::rnorm(200)
  d
}
slope_data <- function(seed, s, family = "gaussian") {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  x0 <- stats::rnorm(300)
  eta <- 0.5 + stats::rnorm(20, 0, 0.5)[g] +
    stats::rnorm(20, 0, 0.5)[g] * x0
  y <- if (family == "gaussian") eta + stats::rnorm(300) else
    stats::rpois(300, exp(eta))
  data.frame(g = g, x = x0 * s, y = y)
}
cases <- list(
  list(id = "bnd2", d = nobs_slope(2), fo = y ~ x + (1 + x | g),
       fam = gaussian(), ctl = frmtmb_control(), ulp = TRUE),
  list(id = "bnd7", d = nobs_slope(7), fo = y ~ x + (1 + x | g),
       fam = gaussian(), ctl = frmtmb_control(), ulp = TRUE),
  list(id = "bnd8", d = nobs_slope(8), fo = y ~ x + (1 + x | g),
       fam = gaussian(), ctl = frmtmb_control(), ulp = TRUE),
  list(id = "gau4", d = slope_data(4, 1e-4), fo = y ~ 1 + (1 + x | g),
       fam = gaussian(), ctl = frmtmb_control(autoscale = FALSE), ulp = TRUE),
  list(id = "gau5", d = slope_data(5, 1e-4), fo = y ~ 1 + (1 + x | g),
       fam = gaussian(), ctl = frmtmb_control(autoscale = FALSE), ulp = TRUE),
  list(id = "pois5", d = slope_data(5, 1e-5, "poisson"),
       fo = y ~ 1 + (1 + x | g), fam = poisson(),
       ctl = frmtmb_control(autoscale = FALSE), ulp = FALSE))
rows <- list()
for (cs in cases) for (k in if (cs$ulp) 0:20 else 0L) {
  d <- cs$d
  d$y <- d$y * (1 + 2e-15 * k)
  r <- run(frm(cs$fo, family = cs$fam, data = d, control = cs$ctl))
  f <- r$fit
  p <- f$opt$par
  j <- grep("theta", names(p))[3]
  lost <- frmtmb:::sdr_of(f)$se_lost
  gb <- if (exists("se_cor_gain_back", ns)) {
    ns$se_cor_gain_back(f, j)
  } else NA_real_
  rows[[length(rows) + 1L]] <- data.frame(
    id = cs$id, k = k, conv = f$opt$convergence, t = signif(p[j], 4),
    past = abs(p[j]) > ns$se_cor_clamp_theta(),
    gain_back = signif(gb, 3),
    lost = paste(paste0(names(lost), "=", lost), collapse = ","),
    bmsg = any(grepl("frmtmb_boundary_fit", r$m, fixed = TRUE)),
    short_w = any(grepl("stopped short of the maximum", r$w, fixed = TRUE)),
    n_w = length(r$w))
}
r <- do.call(rbind, rows)
options(width = 200)
print(r, row.names = FALSE)
for (id in unique(r$id)) {
  s <- r[r$id == id, ]
  cat(sprintf("%-6s fits %2d past clamp %2d boundary msg %2d short warn %2d any warn %2d gain_back>1e-6 %2d lost: %s\n",
              id, nrow(s), sum(s$past), sum(s$bmsg), sum(s$short_w),
              sum(s$n_w > 0), sum(s$gain_back > 1e-6, na.rm = TRUE),
              paste(names(table(s$lost)), table(s$lost), collapse = "; ")))
}
