# What brms 2.23.0 actually returns for the shapes item 2.6f moves.
# Measured against brms itself, not its documentation. The fixtures need
# rename_pars(), which is what brms's own test suite does: the stored
# example carries raw Stan names and no group-level draws until then.
#   Rscript dev/shapes-probe-brms.R > dev/shapes-log/probe-brms.txt 2>&1
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))

f1 <- rename_pars(brms:::brmsfit_example1)
f2 <- rename_pars(brms:::brmsfit_example2)
f4 <- rename_pars(brms:::brmsfit_example4)
f5 <- rename_pars(brms:::brmsfit_example5)
f6 <- rename_pars(brms:::brmsfit_example6)

show <- function(label, expr) {
  cat("\n== ", label, "\n", sep = "")
  r <- tryCatch(withCallingHandlers(expr,
         warning = function(w) {
           cat("  [warning] ", conditionMessage(w), "\n", sep = "")
           invokeRestart("muffleWarning") }),
    error = function(e) {
      cat("  ERROR: ", conditionMessage(e), "\n", sep = ""); NULL })
  invisible(r)
}
dd <- function(x) cat("dim:", paste(dim(x), collapse = " x "),
                      " class:", paste(class(x), collapse = ","), "\n")

show("fitted(f1)", { x <- fitted(f1); dd(x); print(dimnames(x)[[2]])
  print(head(x, 3)) })
show("fitted(f1, summary = FALSE)", dd(fitted(f1, summary = FALSE)))
show("fitted(f1, probs = c(.1,.9))",
     print(colnames(fitted(f1, probs = c(0.1, 0.9)))))
show("fitted(f1, robust = TRUE)", print(colnames(fitted(f1, robust = TRUE))))
show("fitted(f4) ordinal", { x <- fitted(f4); dd(x); print(dimnames(x)[-1]) })
show("fitted(f5) mixture", { x <- fitted(f5); dd(x); print(dimnames(x)[[2]]) })
show("fitted(f6) multivariate", { x <- fitted(f6); dd(x)
  print(dimnames(x)[-1]) })
show("fitted(f1, scale = 'linear')", { x <- fitted(f1, scale = "linear")
  dd(x); print(head(x, 2)) })
show("fitted(f1, dpar = 'sigma')", { x <- fitted(f1, dpar = "sigma")
  dd(x); print(head(x, 2)) })
show("residuals(f1)", { x <- residuals(f1); dd(x); print(colnames(x)) })
show("residuals(f1, summary = FALSE)", dd(residuals(f1, summary = FALSE)))
show("predict(f1)", { x <- predict(f1); dd(x); print(colnames(x))
  print(head(x, 3)) })
show("predict(f4) ordinal", { x <- predict(f4); dd(x); print(colnames(x))
  print(head(x, 3)) })
show("predict(f4, ntrys)", { x <- predict(f4, ntrys = 3); dd(x) })
show("predict(f1, summary = FALSE)", dd(predict(f1, summary = FALSE)))
show("predict(f5) mixture", { x <- predict(f5); dd(x); print(colnames(x)) })
show("predict(f6) multivariate", { x <- predict(f6); dd(x)
  print(dimnames(x)[-1]) })
show("formals predict.brmsfit", print(names(formals(brms:::predict.brmsfit))))
show("formals fitted.brmsfit", print(names(formals(brms:::fitted.brmsfit))))
show("formals residuals.brmsfit",
     print(names(formals(brms:::residuals.brmsfit))))
show("fixef(f1)", { x <- fixef(f1); dd(x); print(dimnames(x))
  print(round(x, 4)) })
show("formals fixef.brmsfit", print(names(formals(brms:::fixef.brmsfit))))
show("fixef(f1, pars = 'Trt1')", print(fixef(f1, pars = "Trt1")))
show("fixef(f1, summary = FALSE)", dd(fixef(f1, summary = FALSE)))
show("ngrps(f1)", str(ngrps(f1)))
show("ngrps(f6)", str(ngrps(f6)))
show("ngrps of a fit with no groups", str(ngrps(f2)))
show("vcov(f1)", { x <- vcov(f1); dd(x); print(rownames(x)) })
show("vcov(f1, correlation = TRUE)",
     print(round(diag(vcov(f1, correlation = TRUE)), 4)))
show("formals vcov.brmsfit", print(names(formals(brms:::vcov.brmsfit))))
show("summary(f1) names", { s <- summary(f1); print(names(s))
  cat("class:", class(s), "\n") })
show("summary(f1)$fixed", { s <- summary(f1); print(dimnames(s$fixed)) })
show("summary(f1)$random", str(summary(f1)$random))
show("summary(f4)$random", str(summary(f4)$random))
show("summary(f1) slot classes", { s <- summary(f1)
  for (n in names(s)) cat(" ", n, ": ", class(s[[n]])[1], " ",
    paste(dim(s[[n]]), collapse = "x"), "\n", sep = "") })
show("formals summary.brmsfit", print(names(formals(brms:::summary.brmsfit))))
show("nsamples(f1)", print(nsamples(f1)))
show("formals nsamples", print(names(formals(getS3method("nsamples",
  "brmsfit", envir = asNamespace("brms"))))))
show("posterior_samples(f1)", { x <- posterior_samples(f1); dd(x)
  print(head(names(x), 14)) })
show("formals posterior_samples", print(names(formals(getS3method(
  "posterior_samples", "brmsfit", envir = asNamespace("brms"))))))
show("posterior_samples(f1, pars = 'b_')",
     print(names(posterior_samples(f1, pars = "b_"))))
show("posterior_samples(f1, 'b_Trt1')",
     print(names(posterior_samples(f1, "b_Trt1"))))
show("posterior_samples(f1, fixed = TRUE)",
     print(names(posterior_samples(f1, "b_Trt1", fixed = TRUE))))
show("variables(f4)", print(variables(f4)))
show("variables(f1)", print(variables(f1)))
show("nobs(f1)", print(nobs(f1)))
