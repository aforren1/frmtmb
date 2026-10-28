# Reviewer: compare the saved draws of review-priordraws-yes.R and
# review-priordraws-only.R, brms against frmtmb.
#   Rscript dev/review-priordraws-compare.R
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
LOG <- "C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/review-priordraws-log"
cat("#### yes: names and two-sample KS (4000 a side)\n")
ps <- c()
for (nm in paste0("T", 1:11)) {
  fb <- file.path(LOG, sprintf("yes-brms-%s.rds", nm))
  ff <- file.path(LOG, sprintf("yes-frm-%s.rds", nm))
  if (!file.exists(fb) || !file.exists(ff)) {
    cat(nm, ": missing", file.exists(fb), file.exists(ff), "\n")
    next
  }
  b <- readRDS(fb); f <- readRDS(ff)
  same <- identical(names(b), names(f))
  cat(sprintf("%s names %s | brms: %s | frm: %s\n", nm,
              if (same) "SAME" else "DIFFER", paste(names(b), collapse = " "),
              paste(names(f), collapse = " ")))
  for (v in intersect(names(b), names(f))) {
    if (length(unique(b[[v]])) < 2) {
      cat(sprintf("   %-16s constant brms %g frm %g\n", v, b[[v]][1], f[[v]][1]))
      next
    }
    p <- suppressWarnings(stats::ks.test(b[[v]], f[[v]])$p.value)
    ps <- c(ps, setNames(p, paste(nm, v)))
    cat(sprintf("   %-16s KS p %.3f  range brms [%.3g, %.3g] frm [%.3g, %.3g]\n",
                v, p, min(b[[v]]), max(b[[v]]), min(f[[v]]), max(f[[v]])))
  }
}
cat(sprintf("SUMMARY yes: %d columns, KS p from %.3f to %.3f, below 0.05: %d, below 0.01: %d\n",
            length(ps), min(ps), max(ps), sum(ps < 0.05), sum(ps < 0.01)))

cat("\n#### only: z-scores of mean and 5/50/95%% quantiles, MCSE from posterior\n")
zz <- c()
for (nm in paste0("O", 1:6)) {
  fb <- file.path(LOG, sprintf("only-brms-%s.rds", nm))
  ff <- file.path(LOG, sprintf("only-frm-%s.rds", nm))
  if (!file.exists(fb) || !file.exists(ff)) {
    cat(nm, ": missing", file.exists(fb), file.exists(ff), "\n")
    next
  }
  b <- readRDS(fb); f <- readRDS(ff)
  cat(nm, "rhat brms", b$rhat, "frm", f$rhat, "\n")
  vs <- intersect(setdiff(names(b$draws), c("lp__", "lprior", ".chain",
                                             ".iteration", ".draw")),
                  names(f$draws))
  vs <- vs[!grepl("^r_|^z_|^Xme", vs) | vs == "r_g[1,Intercept]"]
  one <- function(x, ch) {
    m <- matrix(x, ncol = 4)
    c(mean = mean(x), se = posterior::mcse_mean(m),
      q05 = unname(quantile(x, 0.05)), s05 = unname(posterior::mcse_quantile(m, 0.05)),
      q50 = unname(quantile(x, 0.5)), s50 = unname(posterior::mcse_quantile(m, 0.5)),
      q95 = unname(quantile(x, 0.95)), s95 = unname(posterior::mcse_quantile(m, 0.95)))
  }
  add <- function(lab, xb, xf) {
    a <- one(xb); c <- one(xf)
    z <- c((a["mean"] - c["mean"]) / sqrt(a["se"]^2 + c["se"]^2),
           (a["q05"] - c["q05"]) / sqrt(a["s05"]^2 + c["s05"]^2),
           (a["q50"] - c["q50"]) / sqrt(a["s50"]^2 + c["s50"]^2),
           (a["q95"] - c["q95"]) / sqrt(a["s95"]^2 + c["s95"]^2))
    zz <<- c(zz, setNames(z, paste(nm, lab, c("mean", "q05", "q50", "q95"))))
    cat(sprintf("   %-18s brms mean %8.3f q95 %8.3f | frm mean %8.3f q95 %8.3f | z %s\n",
                lab, a["mean"], a["q95"], c["mean"], c["q95"],
                paste(sprintf("%5.2f", z), collapse = " ")))
  }
  for (v in vs) add(v, b$draws[[v]], f$draws[[v]])
  if (!is.null(b$yrep) && !is.null(f$yrep)) {
    cat("   yrep non-finite: brms", sum(!is.finite(b$yrep)), "frm", sum(!is.finite(f$yrep)), "of", length(f$yrep), "
")
    for (j in 1:3) if (all(is.finite(b$yrep[, j])) && all(is.finite(f$yrep[, j]))) add(paste0("yrep[,", j, "]"), b$yrep[, j], f$yrep[, j])
  } else cat("   yrep missing: brms", is.null(b$yrep), "frm", is.null(f$yrep), "\n")
}
cat(sprintf("SUMMARY only: %d z-scores, max |z| %.2f, |z| > 3: %d\n",
            length(zz), max(abs(zz)), sum(abs(zz) > 3)))
print(round(sort(abs(zz), decreasing = TRUE)[1:5], 2))
