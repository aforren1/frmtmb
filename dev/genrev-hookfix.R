.libPaths(c("C:/Users/adf44/source/r/genrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(...) suppressMessages(suppressWarnings(...))
say <- function(...) cat(sprintf(...))
own <- function(f) environmentName(topenv(environment(f)))
q(library(loo)); q(library(frmtmb))
ns <- asNamespace("frmtmb")
# the one-line fix: set the hook for EVERY owner, not only the ones
# that are not loaded yet.  Simulated here by adding the missing one.
tab <- get("frm_generic_owners", envir = ns)
adopt <- get("frm_adopt", envir = ns)
gens <- names(tab)[vapply(tab, function(o) "loo" %in% o, NA)]
fn <- local({ g <- gens; function(...) {
  for (one in g) tryCatch(adopt(one, tab[[one]], "frmtmb"),
                          error = function(e) NULL)
  invisible() } })
setHook(packageEvent("loo", "onLoad"), fn)
detach("package:loo", unload = FALSE)
q(unloadNamespace("loo"))
q(loadNamespace("brms"))
f <- get("loo", envir = globalenv())
stale <- !identical(topenv(environment(f)), asNamespace("loo"))
tb <- get(".__S3MethodsTable__.", envir = environment(f))
say("WITH the hook set for an already-loaded owner:\n")
say("  loo generic from %s, STALE=%s, loo.brmsfit reachable=%s\n",
    own(f), stale, exists("loo.brmsfit", envir = tb, inherits = FALSE))
say("  identical to live loo::loo: %s\n", identical(f, loo::loo))
cat("GENREVDONE\n")
