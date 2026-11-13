# Lane names: frmtmb's names and layouts on the cases of
# dev/names-cases.R, ML fit and (with "sample") frm_sample() draws.
# Usage: Rscript dev/names-frm.R <base|lane> [sample] [case ...]
# Output: dev/names-out/frm-<arm>[-sample].txt
a <- commandArgs(TRUE)
arm <- a[1]
smp <- "sample" %in% a
sel <- setdiff(a[-1], "sample")
lib <- c(if (arm == "lane") "C:/Users/adf44/source/r/wt-names-lib",
         "C:/Users/adf44/source/r/rellib-r20",
         "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(lib)
suppressPackageStartupMessages(library(frmtmb))
if (smp) suppressPackageStartupMessages(library(frmtmb.sample))
source("dev/names-cases.R")
dir.create("dev/names-out", showWarnings = FALSE)
cases <- names_cases(frmtmb::bf, list(
  categorical = categorical(), skew_normal = skew_normal(),
  gaussian = gaussian(), poisson = poisson()))
if (length(sel)) cases <- cases[sel]
sink(sprintf("dev/names-out/frm-%s%s%s.txt", arm,
             if (smp) "-sample" else "",
             if (length(sel)) paste0("-", paste(sel, collapse = "-")) else
               ""),
     split = TRUE)
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
dshape <- function(x) {
  if (is.null(x)) return("NULL")
  if (is.null(dim(x))) return(paste("vector of", length(x)))
  paste0(paste(dim(x), collapse = " x "), "; dimnames: ",
         if (is.null(dimnames(x))) "NULL (none)" else
           paste(vapply(dimnames(x), function(d) {
             if (is.null(d)) "NULL" else paste(utils::head(d, 6),
                                                collapse = ",")
           }, ""), collapse = " | "))
}
res <- list()
for (nm in names(cases)) {
  cs <- cases[[nm]]
  cat("\n=====", nm, "=====\n")
  fit <- tryCatch(suppressWarnings(suppressMessages(frm(
    cs$f, family = cs$fam, data = cs$d))), error = function(e) e)
  if (inherits(fit, "error")) {
    cat("FIT ERROR:", conditionMessage(fit), "\n")
    next
  }
  obj <- fit
  if (smp) {
    obj <- tryCatch(suppressWarnings(suppressMessages(frm_sample(
      fit, chains = 1, iter = 60, warmup = 30, seed = 1, refresh = 0))),
      error = function(e) e)
    if (inherits(obj, "error")) {
      cat("SAMPLE ERROR:", conditionMessage(obj), "\n")
      next
    }
    cat("as_draws_df names identical to variables():",
        identical(setdiff(names(posterior::as_draws_df(obj)),
                          c(".chain", ".iteration", ".draw")),
                  variables(obj)), "\n")
  }
  v <- variables(obj)
  cat("variables (b[, r_, latent_state, Xme_, Yl_ dropped):\n")
  cat(v[!grepl("^(r_|Xme_|Yl_|latent_state|b\\[)", v)], fill = 78)
  for (p in c("^r_", "^b\\[", "^latent_state", "^zgp", "^simo")) {
    cat(" ", p, ":", head(grep(p, v, value = TRUE), 4),
        "count", sum(grepl(p, v)), "\n")
  }
  re <- tryCatch(ranef(obj), error = function(e) e)
  if (inherits(re, "error")) {
    cat("ranef ERROR:", conditionMessage(re), "\n")
  } else {
    cat("names(ranef):", names(re), "\n")
    for (g in seq_along(re)) cat(" ranef", names(re)[g], ":",
                                 dshape(re[[g]]), "\n")
  }
  co <- tryCatch(coef(obj), error = function(e) e)
  if (inherits(co, "error")) {
    cat("coef ERROR:", conditionMessage(co), "\n")
  } else {
    cat("names(coef):", names(co), "\n")
    for (g in seq_along(co)) {
      cat(" coef", names(co)[g], ":",
          if (is.list(co[[g]]) && !is.data.frame(co[[g]]))
            paste("list:", paste(names(co[[g]]), collapse = ","))
          else dshape(co[[g]]), "\n")
    }
  }
  if (nm %in% c("dist", "distpop")) {
    for (fn in c("posterior_predict", "predictive_error", "log_lik",
                 "posterior_epred", "posterior_linpred")) {
      r <- tryCatch(dshape(get(fn)(obj)), error = function(e)
        paste("ERROR", conditionMessage(e)))
      cat(fn, ":", r, "\n")
    }
    r <- tryCatch(dshape(predict(obj, summary = FALSE)),
                  error = function(e) paste("ERROR", conditionMessage(e)))
    cat("predict(summary = FALSE):", r, "\n")
    if (!smp) {
      cat("summary print, first lines:\n")
      out <- capture.output(print(summary(fit)))
      cat(head(out, 8), sep = "\n")
    }
  }
  if (nm == "dblbar" && !smp) {
    cv <- tryCatch(confint_varcorr(fit), error = function(e) e)
    if (!inherits(cv, "error")) {
      cat("confint_varcorr()$block:", unique(cv$block), "\n")
    }
    vc <- tryCatch(as.data.frame(VarCorr(fit)), error = function(e) e)
    if (!inherits(vc, "error")) print(utils::head(vc))
  }
  res[[nm]] <- list(
    variables = v,
    ranef = if (!inherits(re, "error")) {
      lapply(unclass(re), function(m) dimnames(m)[[length(dim(m))]])
    } else conditionMessage(re),
    coef = if (!inherits(co, "error")) {
      lapply(co, function(m) if (is.data.frame(m) || is.matrix(m))
        colnames(m) else if (is.array(m)) dimnames(m)[[3L]] else names(m))
    } else conditionMessage(co))
}
saveRDS(res, sprintf("dev/names-out/frm-%s%s.rds", arm,
                     if (smp) "-sample" else ""))
sink()
