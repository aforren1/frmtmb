# Reviewer, short check of fix round 2.
.libPaths(c(
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
cat("mgcv after library(frmtmb):", "mgcv" %in% loadedNamespaces(), "\n")
suppressMessages(library(frmtmb.spatial))
cat("mgcv after library(frmtmb.spatial):", "mgcv" %in% loadedNamespaces(),
    "\n")
core <- names(frmtmb:::covstruct_core)
cat("built-ins:", core, "\n")
cat("built-ins in covstruct_formula_terms:",
    intersect(core, frmtmb:::covstruct_formula_terms), "\n")
cat("built-ins in covstruct_reserved_names():",
    intersect(core, frmtmb:::covstruct_reserved_names()), "\n")
cat("registry holds every built-in:", all(core %in% frmtmb:::cs_names()),
    "\n")
cat("cs_check_entry on every registered entry:",
    sum(vapply(frmtmb:::cs_names(), function(n)
      is.null(frmtmb:::cs_check_entry(n, frmtmb:::cs_entry(n))), TRUE)),
    "of", length(frmtmb:::cs_names()), "pass\n")
r <- function(lab, e) cat(sprintf("%-46s -> %s\n", lab, tryCatch({
  f <- suppressWarnings(e)
  paste("FIT", paste(vapply(f$frame$re_blocks, `[[`, "", "covstruct"),
                     collapse = ","), "logLik",
        format(as.numeric(logLik(f)), digits = 10), "b:",
        paste(names(f$estimates$beta), collapse = ","))
}, error = function(e) paste("ERR", substr(conditionMessage(e), 1, 160)))))
set.seed(3)
d <- data.frame(g = gl(15, 6), t = factor(rep(1:3, 30)), x = runif(90),
                z = runif(90, 1, 2))
d$y <- d$x + rnorm(15)[d$g] + rnorm(90)
r("y ~ x + exp(0 + t | g)", frm(y ~ x + exp(0 + t | g), data = d))
r("y ~ exp(x) + exp(0 + t | g)", frm(y ~ exp(x) + exp(0 + t | g), data = d))
r("y ~ exp(x) + (1 | g)", frm(y ~ exp(x) + (1 | g), data = d))
r("y ~ log(z) + sqrt(z) + (1 | g)", frm(y ~ log(z) + sqrt(z) + (1 | g),
                                       data = d))
r("y ~ poly(x, 2) + scale(z) + factor(t) + (1|g)",
  frm(y ~ poly(x, 2) + scale(z) + factor(t) + (1 | g), data = d))
r("y ~ x + gau(0 + t | g)", frm(y ~ x + gau(0 + t | g), data = d))
r("y ~ x + mat(0 + t | g)", frm(y ~ x + mat(0 + t | g), data = d))
r("y ~ x + ou(0 + t | g)", frm(y ~ x + ou(0 + t | g), data = d))
r("y ~ x + s(x) + (1 | g)", frm(y ~ s(z) + (1 | g), data = d))
# the same two fits on base, for the logLik
lp <- function(b, theta, blk) sum(RTMB::dnorm(b, 0, exp(theta[1]), log = TRUE))
vc <- function(theta, blk) matrix(exp(2 * theta[1]), 1, 1,
                                  dimnames = list(blk$cnms, blk$cnms))
for (nm in c("exp", "log", "sqrt", "toeplitz", "spline", "cov", "kernel",
             "power", "rw1"))
  cat(sprintf("  register %-9s %s\n", nm, tryCatch({
    frmtmb_register_covstruct(nm, lpdf = lp, vcov = vc, npar = 1L,
                              sd_idx = 1L); "ok"
  }, error = function(e) paste("refused:", substr(conditionMessage(e), 1, 90)))))
