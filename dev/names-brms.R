# Lane names: brms 2.23.0's names and layouts on the cases of
# dev/names-cases.R. One small fit per case (chains = 1, iter = 60,
# seed 1); only names, dims and orders are read, never values.
# Output: dev/names-out/brms.txt and dev/names-out/brms.rds
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
source("dev/names-cases.R")
dir.create("dev/names-out", showWarnings = FALSE)
cases <- names_cases(brms::bf, list(
  categorical = categorical(), skew_normal = skew_normal(),
  gaussian = gaussian(), poisson = poisson()))
args <- commandArgs(TRUE)
if (length(args)) cases <- cases[args]
sink(sprintf("dev/names-out/brms-%s.txt",
             if (length(args)) paste(args, collapse = "-") else "all"),
     split = TRUE)
cat("brms", format(packageVersion("brms")), "\n")
res <- list()
dshape <- function(x) {
  if (is.null(x)) return("NULL")
  paste0(paste(dim(x), collapse = " x "), "; dimnames: ",
         paste(vapply(dimnames(x), function(d) {
           if (is.null(d)) "NULL" else paste(utils::head(d, 6),
                                              collapse = ",")
         }, ""), collapse = " | "))
}
for (nm in names(cases)) {
  cs <- cases[[nm]]
  cat("\n=====", nm, "=====\n")
  fit <- tryCatch(suppressWarnings(suppressMessages(brm(
    cs$f, family = cs$fam, data = cs$d, chains = 1, iter = 60,
    warmup = 30, seed = 1, refresh = 0, backend = "rstan"))),
    error = function(e) e)
  if (inherits(fit, "error")) {
    cat("FIT ERROR:", conditionMessage(fit), "\n")
    next
  }
  v <- variables(fit)
  v <- v[!grepl("^(z_|L_|Xme_|Yl_|err\\[|zerr\\[|lp__|lprior)", v) |
           grepl("^zgp", v)]
  cat("variables (z_, L_, Xme_, err, lp__, lprior dropped):\n")
  cat(v, fill = 78)
  vall <- variables(fit)
  cat("variables matching err/zerr/zgp/simo (first 4 each):\n")
  for (p in c("^err\\[", "^zerr\\[", "^zgp", "^simo")) {
    cat(" ", p, ":", head(grep(p, vall, value = TRUE), 4),
        "count", sum(grepl(p, vall)), "\n")
  }
  re <- tryCatch(ranef(fit), error = function(e) e)
  if (inherits(re, "error")) {
    cat("ranef ERROR:", conditionMessage(re), "\n")
  } else {
    cat("names(ranef):", names(re), "\n")
    for (g in names(re)) cat(" ranef", g, ":", dshape(re[[g]]), "\n")
  }
  co <- tryCatch(coef(fit), error = function(e) e)
  if (inherits(co, "error")) {
    cat("coef ERROR:", conditionMessage(co), "\n")
  } else {
    cat("names(coef):", names(co), "\n")
    for (g in names(co)) cat(" coef", g, ":", dshape(co[[g]]), "\n")
  }
  if (nm %in% c("dist", "distpop")) {
    cat("posterior_predict:", dshape(posterior_predict(fit)), "\n")
    cat("predictive_error:", dshape(predictive_error(fit)), "\n")
    cat("predict(summary = FALSE):",
        dshape(predict(fit, summary = FALSE)), "\n")
    cat("log_lik:", dshape(log_lik(fit)), "\n")
    cat("posterior_epred:", dshape(posterior_epred(fit)), "\n")
    cat("posterior_linpred:", dshape(posterior_linpred(fit)), "\n")
    cat("summary print, first lines:\n")
    out <- capture.output(print(summary(fit)))
    cat(head(out, 8), sep = "\n")
  }
  res[[nm]] <- list(variables = vall,
                    ranef = if (!inherits(re, "error")) lapply(re, dimnames)
                    else conditionMessage(re),
                    coef = if (!inherits(co, "error")) lapply(co, dimnames)
                    else conditionMessage(co))
}
saveRDS(res, sprintf("dev/names-out/brms-%s.rds",
                     if (length(args)) paste(args, collapse = "-") else
                       "all"))
sink()
