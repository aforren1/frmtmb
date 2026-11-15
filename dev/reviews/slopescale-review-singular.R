# Review of lane slopescale, attack item 4: does counting the far point's
# gain (se_sd_gain_up0(), log sd + 10) give a false "stopped short"
# verdict on a fit whose sd truly is at zero? Genuinely singular data
# (true slope sd 0), gaussian, poisson and bernoulli, (0 + x | g) and
# (1 + x | g), seeds 1..20, at spread 1 under the default and at spread
# 1e-4 under autoscale = FALSE (where the far point reaches out); plus
# bernoulli (1 | g) with true sd 0 and groups of 3 (all-zero groups, the
# Laplace approximation at a huge sd). Base and lane in one process.
# Records every condition and, on the lane, the far point's gain.
# Usage: Rscript slopescale-review-singular.R [lane|both]
a <- commandArgs(TRUE)
mode <- if (length(a)) a[1] else "both"
LANE <- "C:/Users/adf44/source/r/wt-slopescale-lib"
BASE <- "C:/Users/adf44/source/r/rellib-r21"
.libPaths(c(BASE, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
use <- function(lib) {
  if ("frmtmb" %in% loadedNamespaces()) unloadNamespace("frmtmb")
  suppressPackageStartupMessages(library(frmtmb, lib.loc = lib))
}
cat("BLAS:", extSoftVersion()[["BLAS"]], La_library(), "\n")
mk <- function(seed, fam, s) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  z <- rnorm(300)
  eta <- 0.2 + rnorm(20, 0, 0.5)[g] + 0.3 * z
  y <- switch(fam, gaussian = eta + rnorm(300),
              poisson = rpois(300, exp(eta)),
              bernoulli = rbinom(300, 1, plogis(eta)))
  data.frame(g = g, x = s * z, y = y)
}
fam_of <- function(f) switch(f, gaussian = gaussian(), poisson = poisson(),
                             bernoulli = bernoulli())
cases <- list()
for (fam in c("gaussian", "poisson", "bernoulli")) for (form in c("0", "1"))
  for (sp in c("1", "1e-4")) for (seed in 1:20) {
    fo <- if (form == "0") y ~ 1 + (0 + x | g) else y ~ 1 + (1 + x | g)
    cases[[length(cases) + 1L]] <- list(
      nm = paste(fam, form, sp, seed), fam = fam, fo = fo,
      d = mk(seed, fam, as.numeric(sp)),
      ctl = if (sp == "1") list() else list(autoscale = FALSE))
  }
for (seed in 1:20) {
  set.seed(500 + seed)
  g <- factor(rep(1:60, each = 3))
  d <- data.frame(g = g, y = rbinom(180, 1, 0.2))
  cases[[length(cases) + 1L]] <- list(nm = paste("bern3 icpt", seed),
                                      fam = "bernoulli", fo = y ~ 1 + (1 | g),
                                      d = d, ctl = list())
}
run1 <- function(cs, far = FALSE) {
  cond <- character(0)
  f <- withCallingHandlers(
    tryCatch(frm(cs$fo, family = fam_of(cs$fam), data = cs$d,
                 control = do.call(frmtmb_control, cs$ctl)),
             error = function(e) e),
    warning = function(w) {
      cond <<- c(cond, paste0("W:", class(w)[1], ":", conditionMessage(w)))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      cond <<- c(cond, paste0("M:", class(m)[1], ":", conditionMessage(m)))
      invokeRestart("muffleMessage")
    })
  if (inherits(f, "error")) return(list(err = TRUE, cond = cond))
  th <- f$opt$par[names(f$opt$par) == "theta"]
  fg <- NA_real_
  if (far) {
    # the lane's answer of the probe for the last log sd (the slope's,
    # or the only one), recomputed from scratch
    k <- tryCatch(max(frmtmb:::log_sd_theta_index(f)), error = function(e) NA)
    if (is.finite(k)) {
      if (is.environment(f$cache)) f$cache$se_sd_up <- NULL
      g <- tryCatch(frmtmb:::se_sd_gain_up(f, k), error = function(e) NULL)
      fg <- if (is.null(g)) NA_real_ else as.numeric(g)
    }
  }
  list(err = FALSE, cond = cond, ll = as.numeric(logLik(f)), th = th,
       gain = fg)
}
res <- list()
if (mode == "both") {
  use(BASE)
  res$base <- lapply(cases, run1)
}
use(LANE)
res$lane <- lapply(cases, run1, far = TRUE)
short <- function(r) any(grepl("stopped short", r$cond))
zero <- function(r) any(grepl("is at zero", r$cond))
nw <- function(r) sum(startsWith(r$cond, "W:"))
tab <- do.call(rbind, lapply(seq_along(cases), function(i) {
  l <- res$lane[[i]]
  b <- if (!is.null(res$base)) res$base[[i]]
  data.frame(case = cases[[i]]$nm,
             same_cond = if (is.null(b)) NA else identical(b$cond, l$cond),
             same_ll = if (is.null(b)) NA else identical(b$ll, l$ll),
             short_l = short(l), zero_l = zero(l), nw_l = nw(l),
             short_b = if (is.null(b)) NA else short(b),
             zero_b = if (is.null(b)) NA else zero(b),
             nw_b = if (is.null(b)) NA else nw(b),
             last_th = signif(tail(l$th, 1), 4), gain = signif(l$gain, 3))
}))
options(width = 200)
print(tab, row.names = FALSE)
cat("\nfits:", nrow(tab), "; lane 'stopped short':", sum(tab$short_l),
    "; lane 'is at zero':", sum(tab$zero_l),
    "; lane probe gain > 1e-6:", sum(tab$gain > 1e-6, na.rm = TRUE),
    "; max lane gain", max(tab$gain, na.rm = TRUE), "\n")
if (!is.null(res$base)) {
  cat("base 'stopped short':", sum(tab$short_b), "; base 'is at zero':",
      sum(tab$zero_b), "; conditions identical:", sum(tab$same_cond),
      "of", nrow(tab), "; logLik identical:", sum(tab$same_ll), "\n")
  for (i in which(!tab$same_cond)) {
    cat("\n== differs:", tab$case[i], "\n base:\n")
    cat(paste0("  ", substr(res$base[[i]]$cond, 1, 250)), sep = "\n")
    cat(" lane:\n")
    cat(paste0("  ", substr(res$lane[[i]]$cond, 1, 250)), sep = "\n")
  }
}
