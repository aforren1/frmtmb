# What brms 2.23.0 does with an intercept in an ordinal disc formula.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
ns <- asNamespace("brms")
hits <- character()
for (nm in ls(ns, all.names = TRUE)) {
  f <- get(nm, ns); if (!is.function(f)) next
  src <- paste(deparse(f), collapse = "\n")
  if (grepl("disc", src) && grepl("[Ii]ntercept", src) && grepl("warn|stop", src)) hits <- c(hits, nm)
}
print(hits)
set.seed(1)
d <- data.frame(y = sample(0:4, 60, TRUE), z = rnorm(60), x = rnorm(60))
for (fam in list(cumulative(), hurdle_cumulative())) {
  cat("\n==", fam$family, "==\n")
  r <- withCallingHandlers(
    tryCatch(stancode(bf(y ~ x, disc ~ 1 + z), data = transform(d, y = if (fam$family == "cumulative") y + 1 else y), family = fam),
             error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(w) { cat("WARNING:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") },
    message = function(m) { cat("MESSAGE:", conditionMessage(m)); invokeRestart("muffleMessage") })
  if (is.character(r) && length(r) == 1 && startsWith(r, "ERROR")) cat(r, "\n") else
    cat(grep("disc", strsplit(as.character(r), "\n")[[1]], value = TRUE), sep = "\n")
}
