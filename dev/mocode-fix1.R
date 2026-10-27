## Fix round 1 probes: fractional new data, conditional_effects() and
## emmeans reference values, the subset() grid, the gap standard error,
## new-data refusal causes, marginaleffects.
##   Rscript dev/mocode-fix1.R base|lane
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r13",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mocode-lib", base)
} else base)
suppressMessages(library(frmtmb))
cat("== frmtmb", format(packageVersion("frmtmb")), arm,
    dirname(system.file(package = "frmtmb")), "\n")
q <- function(e) suppressMessages(e)
show <- function(lab, e) {
  w <- character(0)
  r <- tryCatch(withCallingHandlers(q(e), warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  }), error = function(e) e)
  cat(sprintf("%-40s %s\n", lab, if (inherits(r, "error")) {
    paste("ERROR:", substr(conditionMessage(r), 1, 200))
  } else paste(signif(as.numeric(unlist(r)), 7), collapse = " ")))
  if (length(w)) cat(sprintf("%-40s warning: %s\n", "", substr(w, 1, 160)))
  invisible(r)
}
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$z <- rnorm(n)
d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] + 0.5 * d$z + rnorm(n)
fk <- q(frm(y ~ mo(k), data = d))
cat("-- 1. fractional and out-of-range new data, y ~ mo(k), k in 1..4\n")
for (v in c(0.6, 2.4, 2.5, 4.4, 0, 5, 1, 4)) {
  show(paste("fitted(newdata k =", v, ")"),
       fitted(fk, newdata = data.frame(k = v))[, 1])
}
cat("-- 8. new-data refusal causes\n")
show("k = NA", fitted(fk, newdata = data.frame(k = NA_real_)))
show("k = ordered factor", fitted(fk, newdata = data.frame(
  k = factor(2, levels = 1:4, ordered = TRUE))))
show("k = \"2\"", fitted(fk, newdata = data.frame(k = "2")))
show("y ~ mo(mi(k)) fit", frm(y ~ mo(mi(k)), data = d))
cat("-- 2. conditional_effects() and emmeans, y ~ mo(k) + z\n")
fz <- q(frm(y ~ mo(k) + z, data = d))
ce <- q(conditional_effects(fz, effects = "z"))[[1]]
cat("effects = z: k held at", unique(ce$k), "; first z", ce$z[1],
    "estimate", signif(ce$estimate__[1], 7), "\n")
cat("frmtmb fitted at k = 3, that z:",
    signif(fitted(fz, newdata = data.frame(k = 3, z = ce$z[1]))[, 1], 7),
    "\n")
ce2 <- q(conditional_effects(fz, effects = "z:k"))[[1]]
cat("effects = z:k: k values", unique(ce2$k), "; effect2__ labels",
    unique(as.character(ce2$effect2__)), "\n")
ce3 <- q(conditional_effects(fz, effects = "k"))[[1]]
cat("effects = k: grid", ce3$k, "; z held at", unique(ce3$z), "\n")
cat("   bands lower/upper:", signif(ce3$lower__, 6), "/",
    signif(ce3$upper__, 6), "\n")
if (requireNamespace("emmeans", quietly = TRUE)) {
  show("emmeans(fz, ~ z)", summary(emmeans::emmeans(fz, ~ z))$emmean)
  show("emmeans(fz, ~ z, at = list(k = 3))",
       summary(emmeans::emmeans(fz, ~ z, at = list(k = 3)))$emmean)
}
if (requireNamespace("marginaleffects", quietly = TRUE)) {
  show("avg_slopes(fz, variables = 'k')",
       marginaleffects::avg_slopes(fz, variables = "k")$estimate)
  show("avg_slopes(fz, variables = 'z')",
       marginaleffects::avg_slopes(fz, variables = "z")$estimate)
  show("avg_comparisons(fz, variables = list(k = 1))",
       marginaleffects::avg_comparisons(fz,
         variables = list(k = 1))$estimate)
}
cat("-- 3. subset() response\n")
ds <- d
ds$y2 <- ds$y + rnorm(n)
ds$s2 <- TRUE
fs <- q(frm(bf(y | subset(k > 1) ~ mo(k)) + bf(y2 | subset(s2) ~ mo(k)),
            data = ds))
show("ce(effects = k, resp = y) grid",
     q(conditional_effects(fs, effects = "k", resp = "y"))[[1]]$k)
show("ce(effects = k, resp = y) estimate",
     q(conditional_effects(fs, effects = "k", resp = "y"))[[1]]$estimate__)
show("ce(effects = k, resp = y2) grid",
     q(conditional_effects(fs, effects = "k", resp = "y2"))[[1]]$k)
cat("-- 6. empty interior category: x in 1, 2, 5, seed 7\n")
set.seed(7)
dg <- data.frame(x = sample(c(1L, 2L, 5L), 200, replace = TRUE))
dg$y <- dg$x + rnorm(200)
fg <- suppressWarnings(q(frm(y ~ mo(x), data = dg)))
show("fitted Estimate x = 1..5",
     fitted(fg, newdata = data.frame(x = 1:5))[, 1])
show("fitted Est.Error x = 1..5",
     fitted(fg, newdata = data.frame(x = 1:5))[, 2])
show("predict(se.fit) x = 1..5", frm_linpred(fg, newdata = data.frame(
  x = 1:5), se.fit = TRUE)$se.fit)
cg <- suppressWarnings(q(conditional_effects(fg, effects = "x")))[[1]]
show("ce bands se__ x = 1..5", cg$se__)
cat("-- SE with the simplex, seed 5 data (brms posterior sd of fitted at",
    "k = 1..4 is in dev/mocode-brms-fitted.R's log)\n")
show("fitted Est.Error k = 1..4",
     fitted(fk, newdata = data.frame(k = 1:4))[, 2])
cat("-- fitted Est.Error of y ~ mo(k) + z at k = 1..4, z = 0 (brms: see",
    "dev/mocode-out/brms-fix1.txt)\n")
show("fitted Estimate", fitted(fz, newdata = data.frame(k = 1:4, z = 0))[, 1])
show("fitted Est.Error", fitted(fz, newdata = data.frame(k = 1:4, z = 0))[, 2])
