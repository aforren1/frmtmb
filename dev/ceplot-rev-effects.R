# Reviewer check (lane ceplot): which named effects conditional_effects()
# takes, per model shape, against brms 2.23.0's rule
# (get_all_effects(brmsterms(), comb_all = TRUE), read without a fit).
# For every data variable v, and a pair of two model variables:
#   lane/base: "ok" (answered, no invalid-effect warning), "INVALID"
#   (dropped with the warning or the all-invalid error), or the error.
#   Rscript dev/ceplot-rev-effects.R lane|base > dev/ceplot-rev-log/effects-<arm>.txt
# Data seed 11.
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/wt-ceplot-lib",
          "C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressMessages(library(frmtmb))
has_brms <- requireNamespace("brms", quietly = TRUE)
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
set.seed(11)
n <- 240
d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n), extra = rnorm(n),
                f = factor(sample(c("a", "b"), n, TRUE)),
                g = factor(rep(1:12, 20)), h = factor(rep(1:8, 30)),
                o = sample(1:4, n, TRUE), sdx = runif(n, 0.1, 0.3),
                wt = runif(n, 0.5, 1.5), off = rnorm(n, 0, 0.1),
                t = rep(1:20, 12))
d$g1 <- factor(sample(1:12, n, TRUE), levels = 1:12)
d$g2 <- factor(sample(1:12, n, TRUE), levels = 1:12)
d$xm <- d$x + rnorm(n, 0, 0.2)
d$y <- rnorm(n, 1 + 0.5 * d$x + 0.3 * d$z + rnorm(12)[d$g])
d$yc <- rpois(n, exp(0.2 + 0.3 * d$x))
d$yo <- factor(cut(d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf)), ordered = TRUE)
d$yp <- abs(d$y) + 0.5
d$fg <- factor(ifelse(as.integer(d$g) <= 6, "a", "b"))
d$y2 <- rnorm(n, 0.4 * d$z)
d$xmi <- d$x
d$xmi[c(5, 17, 40)] <- NA
cases <- list(
  plain = list(bf(y ~ x + z), gaussian()),
  inter = list(bf(y ~ x * f), gaussian()),
  smooth = list(bf(y ~ s(x) + z), gaussian()),
  smooth_by = list(bf(y ~ s(x, by = f) + f), gaussian()),
  t2 = list(bf(y ~ t2(x, z)), gaussian()),
  gp = list(bf(y ~ gp(x) + z), gaussian()),
  gp_by = list(bf(y ~ gp(x, by = f) + f), gaussian()),
  mo = list(bf(y ~ mo(o) + x), gaussian()),
  me = list(bf(y ~ me(xm, sdx) + z), gaussian()),
  mi = list(bf(y ~ mi(xmi) + w) + bf(xmi | mi() ~ w), gaussian()),
  nl = list(bf(yp ~ a * exp(b * x), a ~ 1 + z, b ~ 1, nl = TRUE),
            gaussian()),
  dpar = list(bf(y ~ x, sigma ~ z), gaussian()),
  reslope = list(bf(y ~ x + (1 + w | g)), gaussian()),
  offset = list(bf(yc ~ x + offset(off)), poisson()),
  weights = list(bf(y | weights(wt) ~ x), gaussian()),
  cs = list(bf(yo ~ cs(x) + z), acat()),
  poly = list(bf(y ~ poly(x, 2) + log(abs(z) + 1)), gaussian()),
  grby = list(bf(y ~ x + (1 | gr(g, by = fg))), gaussian()),
  mm = list(bf(y ~ x + (1 | mm(g1, g2))), gaussian()),
  gh = list(bf(y ~ x + (1 | g:h)), gaussian()),
  mix = list(bf(y ~ 1, mu1 ~ x, mu2 ~ z),
             mixture(gaussian(), gaussian())),
  zi = list(bf(yc ~ x, zi ~ z), zero_inflated_poisson()),
  ar = list(bf(y ~ x + ar(time = t, gr = g)), gaussian()),
  fs = list(bf(y ~ x + s(z, g, bs = "fs", k = 4)), gaussian()),
  mv = list(bf(y ~ x) + bf(y2 ~ z), gaussian()),
  lf = list(bf(y ~ x) + lf(sigma ~ w), gaussian())
)
cand <- c("x", "z", "w", "extra", "f", "g", "h", "o", "sdx", "wt", "off",
          "t", "g1", "g2", "xm", "fg", "xmi")
try_eff <- function(fit, eff, resp) {
  warned <- character(0)
  r <- tryCatch(withCallingHandlers(
    conditional_effects(fit, effects = eff, resolution = 3, resp = resp),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(r, "error")) {
    m <- conditionMessage(r)
    if (grepl("invalid for this model", m)) return("INVALID")
    return(paste0("ERR[", substr(gsub("\n", " ", m), 1, 60), "]"))
  }
  if (any(grepl("invalid for this model", warned))) return("INVALID")
  "ok"
}
brms_valid <- function(bform, fam) {
  if (!has_brms) return(NULL)
  bt <- tryCatch({
    bf2 <- brms::bf(bform$formula %||% bform)
    NULL
  }, error = function(e) NULL)
  NULL
}
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fit <- tryCatch(suppressWarnings(frm(cs[[1]], family = cs[[2]], data = d)),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("%-10s FIT ERROR %s\n", nm, conditionMessage(fit)))
    next
  }
  resp <- names(fit$spec$responses)[1L]
  res <- vapply(cand, function(v) try_eff(fit, v, resp), "")
  cat(sprintf("%-10s %s\n", nm, paste0(cand, "=", res, collapse = " ")))
}
