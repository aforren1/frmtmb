# Every S3 method frmtmb registers for frmtmb_fit, called on a saved
# cox() and car() fit read where neither extension can be loaded (both
# installed as unloadable stubs), or where both are loaded. Each call is
# classified: refused naming both packages, refused naming one, answered
# (with a count of NA in the answer), or another error. Review finding
# 1 of dev/reviews/2026-10-24-smallfix.md.
#   Rscript dev/smallfix-methods-probe.R save <rds> <lib;lib>
#   Rscript dev/smallfix-methods-probe.R read <rds> <lib;lib> <none|all>
# Log: dev/smallfix-out/methods-probe-<arm>.txt
a <- commandArgs(trailingOnly = TRUE)
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
libs <- c(strsplit(a[3], ";", fixed = TRUE)[[1]], user)
if (a[1] == "save") {
  .libPaths(libs)
  suppressMessages({
    library(frmtmb); library(frmtmb.survival); library(frmtmb.spatial)
  })
  g <- expand.grid(r = 1:3, c = 1:3)
  n <- nrow(g)
  W <- 1 * (as.matrix(dist(g, method = "manhattan")) == 1)
  dimnames(W) <- list(paste0("L", 1:n), paste0("L", 1:n))
  set.seed(20261024)
  dd <- data.frame(loc = factor(rep(rownames(W), each = 12),
                                levels = rownames(W)), x = rnorm(12 * n))
  u <- rnorm(n, 0, 0.5)
  dd$time <- rexp(nrow(dd), exp(-0.5 + 0.6 * dd$x + u[dd$loc]))
  dd$out <- rbinom(nrow(dd), 1, 0.25)
  fit <- suppressWarnings(suppressMessages(frm(
    bf(time | cens(out) ~ x + car(W, gr = loc, type = "escar")),
    family = cox(df = 4), data = dd, data2 = list(W = W))))
  saveRDS(list(fit = fit, dd = dd, W = W), a[2])
  cat("saved\n")
  quit(save = "no")
}
need <- c("frmtmb.survival", "frmtmb.spatial")
if (identical(a[4], "none")) {
  broken <- tempfile("broken-lib")
  for (p in need) {
    dir.create(file.path(broken, p), recursive = TRUE)
    writeLines(c(paste("Package:", p), "Version: 0.0.0"),
               file.path(broken, p, "DESCRIPTION"))
  }
  libs <- c(broken, libs)
}
.libPaths(libs)
suppressMessages(library(frmtmb))
if (identical(a[4], "all")) {
  suppressMessages({library(frmtmb.survival); library(frmtmb.spatial)})
}
cat("frmtmb from", dirname(find.package("frmtmb")), "| loadable:",
    vapply(need, requireNamespace, NA, quietly = TRUE), "\n")
s <- suppressWarnings(readRDS(a[2]))
f <- s$fit
dd <- s$dd
ns <- readLines(file.path(find.package("frmtmb"), "NAMESPACE"))
m <- regmatches(ns, regexpr("^S3method\\(([^,]+),frmtmb_fit\\)", ns))
gens <- unique(sub("^.*::", "", sub("^S3method\\(([^,]+),frmtmb_fit\\)$",
                                     "\\1", m)))
gens <- sort(gens)
args <- list(
  hypothesis = list("x > 0"), refit = list(s$dd$time),
  update = list(. ~ . - x), drop1 = list(), simulate = list(nsim = 1),
  predict = list(newdata = dd[1:3, ]), anova = list(),
  loo_model_weights = NULL, kfold = NULL, reloo = NULL,
  loo_moment_match = NULL, loo_subsample = NULL, emm_basis = NULL,
  recover_data = NULL, set_coef = NULL, profile = NULL, plot = NULL,
  pp_check = NULL, expose_functions = NULL, influence = NULL)
cls <- function(e) {
  msg <- gsub("\\s+", " ", conditionMessage(e))
  if (inherits(e, "frmtmb_extension_missing")) {
    if (grepl("needs 2 extension packages", msg)) "REFUSED-BOTH" else
      "REFUSED-ONE"
  } else paste("ERROR", substr(msg, 1, 110))
}
for (gn in gens) {
  if (gn %in% names(args) && is.null(args[[gn]])) {
    cat(sprintf("%-26s SKIPPED (needs a setup the probe has not)\n", gn))
    next
  }
  fn <- tryCatch(getExportedValue("frmtmb", gn), error = function(e) {
    tryCatch(get(gn, mode = "function"), error = function(e) NULL)
  })
  if (is.null(fn)) {
    cat(sprintf("%-26s NO GENERIC\n", gn))
    next
  }
  r <- tryCatch({
    v <- suppressWarnings(suppressMessages(
      utils::capture.output(out <- do.call(fn, c(list(f), args[[gn]])))))
    nas <- tryCatch(sum(is.na(unlist(out))), error = function(e) NA)
    paste("OK", class(out)[1], "NA:", nas)
  }, error = cls)
  cat(sprintf("%-26s %s\n", gn, r))
}
