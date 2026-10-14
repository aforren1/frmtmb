# Lane dblbar: the `||` fits against the base build, against the same
# model spelled (x | gr(g, cor = FALSE)) and diag(x | g), and against
# brms 2.23.0's posterior (dev/dblbar-brms-fit.R). Reads
# dev/dblbar-log/fits-{base,lane}.rds and brms-fit-<case>.rds.
# Output: dev/dblbar-log/fits-compare.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
lg <- function(f) file.path(wt, "dev/dblbar-log", f)
base <- readRDS(lg("fits-base.rds"))
lane <- readRDS(lg("fits-lane.rds"))
sink(lg("fits-compare.txt"), split = TRUE)
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))
ulps <- function(a, b) max(abs(a - b) / (.Machine$double.eps * abs(b)))
for (nm in names(lane)) {
  cat("\n######## ", nm, "\n", sep = "")
  L <- lane[[nm]]
  B <- base[[nm]]
  ref <- L[[1]]
  for (k in seq_along(L)[-1]) {
    o <- L[[k]]
    if (!is.null(o$error)) {
      cat(names(L)[k], ": ERROR ", o$error, "\n", sep = "")
      next
    }
    same <- names(o$vals)[names(o$vals) %in% names(ref$vals)]
    cat(sprintf("lane %s vs %s: identical vals %s, identical logLik %s, max rel %.3g\n",
                names(L)[1], names(L)[k],
                identical(unname(ref$vals[same]), unname(o$vals[same])),
                identical(ref$logLik, o$logLik),
                rel(ref$vals[same], o$vals[same])))
    for (p in c("p_gau", "p_pop", "p_int", "p_all")) {
      if (is.character(ref[[p]]) || is.character(o[[p]])) {
        cat("  ", p, ": ", if (is.character(ref[[p]])) ref[[p]] else
          "ok", " / ", if (is.character(o[[p]])) o[[p]] else "ok", "\n",
          sep = "")
      } else {
        cat(sprintf("  %s identical %s\n", p, identical(ref[[p]], o[[p]])))
      }
    }
  }
  b1 <- B[[1]]
  if (!is.null(b1$error)) {
    cat("base ", names(B)[1], ": ERROR ", b1$error, "\n", sep = "")
  } else {
    cat(sprintf("base vs lane %s: logLik base %.15g lane %.15g (diff %.3g)\n",
                names(L)[1], b1$logLik, ref$logLik, ref$logLik - b1$logLik))
    cat("  base names:", paste(names(b1$vals), collapse = " "), "\n")
    cat("  lane names:", paste(names(ref$vals), collapse = " "), "\n")
    same <- intersect(names(ref$vals), names(b1$vals))
    cat(sprintf("  shared values: max rel %.3g, max ulps %.3g\n",
                rel(ref$vals[same], b1$vals[same]),
                ulps(ref$vals[same], b1$vals[same])))
    for (p in c("p_gau", "p_pop", "p_int", "fitted_in")) {
      if (!is.character(b1[[p]]) && !is.character(ref[[p]])) {
        cat(sprintf("  %s: max rel %.3g\n", p, rel(ref[[p]], b1[[p]])))
      } else {
        cat("  ", p, ": base ", if (is.character(b1[[p]])) b1[[p]] else
          "ok", " / lane ", if (is.character(ref[[p]])) ref[[p]] else
          "ok", "\n", sep = "")
      }
    }
  }
  bf <- lg(sprintf("brms-fit-%s.rds", nm))
  if (file.exists(bf)) {
    bb <- readRDS(bf)
    v <- intersect(names(bb$median), names(ref$vals))
    z <- (ref$vals[v] - bb$median[v]) / bb$sd[v]
    cat("  brms parameters:", paste(names(bb$median), collapse = " "), "\n")
    cat("  (ML - brms median) / brms sd:\n")
    print(round(z, 3))
    cat("  brms rhat max:", format(max(bb$rhat), digits = 4), "\n")
    cat("  summary random rows equal brms's:",
        identical(unname(unlist(ref$summary_random)),
                  unname(unlist(bb$summary_random))), "\n")
    cat("  VarCorr g elements: frm", paste(names(ref$vc$g), collapse = ","),
        " brms", paste(names(bb$vc$g), collapse = ","), "\n")
    cat("  VarCorr g sd rownames equal:",
        identical(rownames(ref$vc$g$sd), rownames(bb$vc$g$sd)), "\n")
    # frmtmb's ranef() is a level-by-coefficient matrix and brms's an
    # array with a summary dimension between; levels and names compare
    dn <- function(x) list(x[[1L]], x[[length(x)]])
    cat("  ranef levels and coefficient names equal:",
        identical(dn(ref$ranef_dimnames$g), dn(bb$ranef_dimnames$g)),
        "  coef():", identical(dn(ref$coef_dimnames$g),
                                dn(bb$coef_dimnames$g)), "\n")
    cat("  new-level fitted (gaussian), frm Estimate and se vs brms mean",
        "and sd:\n")
    print(cbind(frm = ref$p_gau[, 1], brms = bb$newlev$gau_mean,
                frm_se = ref$p_gau[, 2], brms_sd = bb$newlev$gau_sd))
    cat("  re_formula = NA at the new rows, frm vs brms:\n")
    print(cbind(frm = ref$p_pop[, 1], brms = bb$newlev$pop_mean))
    cat("  re_formula = ~ (1 | g) on rows 1:6, frm vs brms mean:\n")
    print(cbind(frm = ref$p_int[, 1], brms = bb$re_formula_int))
  }
}
sink()
