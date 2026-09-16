## Does brms carry the arguments the new refusal now rejects?
.libPaths(c("C:/Users/adf44/source/r/asrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
gs <- c("print", "summary", "model.frame", "nobs", "coef", "logLik",
        "formula", "vcov", "family", "terms", "update", "fixef", "ranef",
        "VarCorr", "residuals", "fitted", "predict")
for (g in gs) {
  m <- tryCatch(getFromNamespace(paste0(g, ".brmsfit"), "brms"),
                error = function(e) NULL)
  cat(sprintf("%-14s %s\n", g,
              if (is.null(m)) "<no brmsfit method>"
              else paste(names(formals(m)), collapse = ", ")))
}
for (nm in c("print.brmssummary", "print.brmsfit", "summary.brmsfit")) {
  m <- tryCatch(getFromNamespace(nm, "brms"), error = function(e) NULL)
  if (!is.null(m)) cat(sprintf("%-20s %s\n", nm,
                               paste(names(formals(m)), collapse = ", ")))
}
cat("\n-- who calls these generics with the refused argument? --\n")
probe <- list(
  c("nobs", "use.fallback"), c("coef", "complete"),
  c("print", "digits"), c("model.frame", "data"),
  c("summary", "digits"), c("logLik", "REML"))
pkgs <- rownames(utils::installed.packages())
pkgs <- intersect(pkgs, c("stats", "emmeans", "insight", "marginaleffects",
                          "broom", "broom.mixed", "lme4", "nlme", "car",
                          "sandwich", "lmtest", "MuMIn", "performance",
                          "parameters", "bayesplot", "loo", "posterior",
                          "DHARMa", "ggeffects", "modelbased", "effects",
                          "multcomp", "mgcv", "glmmTMB", "brms", "knitr",
                          "rmarkdown", "testthat", "pkgdown"))
for (pk in pkgs) {
  ok <- tryCatch({ loadNamespace(pk); TRUE }, error = function(e) FALSE)
  if (!ok) next
  ns <- asNamespace(pk)
  for (pr in probe) {
    pat <- paste0(pr[1], "\\(.*", pr[2], " *=")
    hits <- character()
    for (nm in ls(ns, all.names = TRUE)) {
      o <- tryCatch(get(nm, envir = ns), error = function(e) NULL)
      if (!is.function(o)) next
      src <- tryCatch(deparse(body(o)), error = function(e) character())
      h <- grep(pat, src, value = TRUE)
      if (length(h)) hits <- c(hits, sprintf("%s::%s | %s", pk, nm,
                                             trimws(h[1L])))
    }
    if (length(hits)) {
      cat(sprintf("\n[%s(%s=)] %d site(s)\n", pr[1], pr[2], length(hits)))
      cat(paste0("   ", utils::head(hits, 4)), sep = "\n")
    }
  }
}
cat("\nDONE\n")
