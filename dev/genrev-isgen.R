.libPaths(c("C:/Users/adf44/source/r/genrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(...) suppressMessages(suppressWarnings(...))
for (p in c("posterior","loo","rstantools","bayesplot","brms","lme4",
            "nlme","generics")) q(loadNamespace(p))
ns <- asNamespace("frmtmb")
tab <- get("frm_generic_owners", envir = ns)
say <- function(...) cat(sprintf(...))
bad <- 0L
for (g in names(tab)) {
  for (o in tab[[g]]) {
    f <- tryCatch(getExportedValue(o, g), error = function(e) NULL)
    if (!is.function(f)) { say("%-20s %-10s NOT A FUNCTION\n", g, o); bad <- bad+1L; next }
    src <- paste(deparse(body(f)), collapse = " ")
    isg <- grepl("UseMethod", src, fixed = TRUE) ||
           grepl("standardGeneric", src, fixed = TRUE)
    if (!isg) { say("%-20s %-10s NOT A GENERIC: %s\n", g, o, substr(src,1,60)); bad <- bad+1L }
  }
}
say("owner exports that are not S3/S4 generics: %d of %d checks\n", bad,
    sum(lengths(tab)))
cat("GENREVDONE\n")
