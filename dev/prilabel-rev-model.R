# Review of lane prilabel: does the nested-group expansion keep the
# MODEL? Fits the same nested models on rellib-r9 ("base") and on the
# lane library ("lane"); the lane arm reads the base arm's saved fits
# and compares logLik, fixef, sds, ranef by level and fitted values.
#   arg 1: "base" or "lane". Data dev/dblbar-data.R seed 20261013,
#   g3 seed 7, k1/k2 seed 9.
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-prilabel-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages(library(frmtmb))
lg <- file.path(wt, "dev/prilabel-rev-log")
sink(file.path(lg, paste0("model-", arm, ".txt")), split = TRUE)
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
dd <- dblbar_data()
set.seed(7)
dd$g3 <- factor(sample(c("r1", "r2"), nrow(dd), TRUE))
set.seed(9)
# levels that carry "_" and ":" themselves; a_b x c and a x b_c both
# paste to a_b_c
dd$k1 <- factor(sample(c("a_b", "a", "c:d"), nrow(dd), TRUE))
dd$k2 <- factor(sample(c("c", "b_c", "e"), nrow(dd), TRUE))
try_ <- function(expr) tryCatch(expr, error = function(e) {
  paste("ERROR:", conditionMessage(e))
})
forms <- list(
  M1 = bf(y ~ x + (1 + x | g1/g2)),
  M2 = bf(y ~ x + (1 + x | g1) + (1 + x | g1:g2)),
  M3 = bf(y ~ x + (1 | g1/g2/g3)),
  M4 = bf(y ~ x + (1 | g1/g2), sigma ~ (1 | g1/g2)),
  M5 = bf(y ~ x + (1 | k1/k2)),
  M6 = bf(y ~ x + (1 + x | g1) + (1 + x | g2:g1))
)
fits <- lapply(forms, function(f) try_(frm(f, data = dd)))
summ <- lapply(names(fits), function(nm) {
  ft <- fits[[nm]]
  if (is.character(ft)) return(list(err = ft))
  re <- ranef(ft)
  list(ll = as.numeric(logLik(ft)), fe = fixef(ft)[, 1],
       vc = lapply(VarCorr(ft), function(v) v$sd[, 1]),
       re = lapply(re, function(a) {
         if (length(dim(a)) == 3L) {
           m <- a[, 1L, , drop = FALSE]
           matrix(m, dim(a)[1L], dimnames = list(dimnames(a)[[1L]],
                                                 dimnames(a)[[3L]]))
         } else {
           as.matrix(a)
         }
       }),
       fit = fitted(ft)[, 1], conv = ft$opt$convergence,
       vars = grep("^(sd|cor)_", variables(ft), value = TRUE))
})
names(summ) <- names(fits)
saveRDS(summ, file.path(lg, paste0("model-", arm, ".rds")))
for (nm in names(summ)) {
  s <- summ[[nm]]
  cat("\n==", nm, deparse1(forms[[nm]]$formula), "\n")
  if (!is.null(s$err)) { cat(s$err, "\n"); next }
  cat(sprintf("logLik %.10f conv %s\n", s$ll, s$conv))
  cat("groups:", names(s$re), "\n")
  print(s$vars)
  for (g in names(s$re)) {
    rn <- if (is.matrix(s$re[[g]])) rownames(s$re[[g]]) else
      names(s$re[[g]])
    cat(" ", g, ": n =", length(rn), "dup =", sum(duplicated(rn)),
        "head:", head(rn, 4), "\n")
  }
}
if (arm == "lane") {
  base <- readRDS(file.path(lg, "model-base.rds"))
  rel <- function(a, b) max(abs(a - b) / pmax(1, abs(b)))
  cat("\n######## lane against base\n")
  for (nm in names(summ)) {
    a <- summ[[nm]]
    b <- base[[nm]]
    if (!is.null(a$err) || !is.null(b$err)) {
      cat(nm, "error on an arm\n"); next
    }
    cat(sprintf("%s dlogLik %.3e  max rel dfixef %.3e  max rel dfitted %.3e\n",
                nm, a$ll - b$ll, rel(a$fe, b$fe), rel(a$fit, b$fit)))
    # sds matched by group: base's g2:g1 is the lane's g1:g2
    flip <- function(g) {
      if (!grepl(":", g)) return(g)
      paste(rev(strsplit(g, ":", fixed = TRUE)[[1]]), collapse = ":")
    }
    for (g in setdiff(names(a$vc), "residual__")) {
      gb <- if (g %in% names(b$vc)) g else flip(g)
      if (!gb %in% names(b$vc)) { cat("  ", g, "no base group\n"); next }
      cat(sprintf("   sd %-10s (base %-10s) max rel diff %.3e\n", g, gb,
                  rel(a$vc[[g]], b$vc[[gb]])))
    }
  }
  cat("\n## lane M1 against base M6 (hand-written g2:g1) and lane M2\n")
  cat(sprintf("lane M1 - base M6 logLik %.3e\n", summ$M1$ll - base$M6$ll))
  cat(sprintf("lane M1 - lane M2 logLik %.3e, identical vars %s\n",
              summ$M1$ll - summ$M2$ll,
              identical(summ$M1$vars, summ$M2$vars)))
  cat(sprintf("lane M1 - lane M2 max rel dfitted %.3e\n",
              rel(summ$M1$fit, summ$M2$fit)))
  # ranef by level: base q1:p1 against lane p1_q1
  rb <- base$M1$re[["g2:g1"]]
  ra <- summ$M1$re[["g1:g2"]]
  key_b <- vapply(strsplit(rownames(rb), ":", fixed = TRUE),
                  function(p) paste(rev(p), collapse = "_"), "")
  cat(sprintf("ranef g1:g2 by level: matched %d of %d, max rel diff %.3e\n",
              sum(key_b %in% rownames(ra)), nrow(ra),
              rel(ra[key_b, , drop = FALSE], rb)))

  cat("\n## predictions and re_formula on lane M1\n")
  f1 <- fits$M1
  nd <- dd[1:6, ]
  full <- fitted(f1, newdata = nd)[, 1]
  cat(sprintf("newdata (old levels) vs fitted: max diff %.3e\n",
              max(abs(full - fitted(f1)[1:6, 1]))))
  for (rf in list(~ (1 + x | g1/g2), ~ (1 + x | g1) + (1 + x | g1:g2),
                  ~ (1 + x | g1), ~ (1 + x | g1:g2), ~ (1 + x | g2:g1))) {
    p <- try_(fitted(f1, newdata = nd, re_formula = rf)[, 1])
    cat(deparse1(rf), ":", if (is.character(p)) p else
      sprintf("max diff from full %.3e", max(abs(p - full))), "\n")
  }
  # the g1 part alone by hand: fixed + g1 intercept + slope
  re1 <- ranef(f1)$g1[, 1, ]
  fe <- fixef(f1)[, 1]
  hand <- fe[1] + fe[2] * nd$x + re1[as.character(nd$g1), 1] +
    re1[as.character(nd$g1), 2] * nd$x
  p1 <- fitted(f1, newdata = nd, re_formula = ~ (1 + x | g1))[, 1]
  cat(sprintf("re_formula ~ (1 + x | g1) against fixef + ranef g1: %.3e\n",
              max(abs(p1 - hand))))
  cat("\n## M5 levels with _ and :\n")
  f5 <- fits$M5
  if (!is.character(f5)) {
    r5 <- ranef(f5)
    print(lapply(r5, function(a) dimnames(a)[[1]]))
    print(try_(head(coef(f5)[[1]], 3)))
    nd5 <- dd[1:3, ]
    nd5$k2 <- factor("zz")
    cat(try_(fitted(f5, newdata = nd5)), "\n")
  }
}
sink()
