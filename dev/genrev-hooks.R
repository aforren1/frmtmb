LIB <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(...) suppressMessages(suppressWarnings(...))
say <- function(...) cat(sprintf(...))
mine <- function(p) {
  h <- getHook(packageEvent(p, "onLoad"))
  n <- sum(vapply(h, function(f) {
    e <- environment(f)
    !is.null(e) && exists("g", envir = e, inherits = FALSE) &&
      exists("p", envir = e, inherits = FALSE)
  }, NA))
  c(total = length(h), frmtmb_shaped = n)
}
if (!requireNamespace("pkgload", quietly = TRUE)) { say("no pkgload\n"); quit(save="no") }
q(library(brms))
say("lme4 loaded with brms: %s\n", isNamespaceLoaded("lme4"))
for (p in c("lme4","bayesplot","posterior")) {
  v <- mine(p); say("before  %-10s total=%d frmtmb-shaped=%d\n", p, v[1], v[2])
}
for (i in 1:3) {
  q(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-generics",
                      export_all = FALSE, quiet = TRUE, helpers = FALSE))
  for (p in c("lme4","bayesplot","posterior")) {
    v <- mine(p); say("load_all%d %-10s total=%d frmtmb-shaped=%d\n", i, p, v[1], v[2])
  }
}
# does the CURRENT namespace still know about all of them?
ns <- asNamespace("frmtmb")
he <- get("frm_generic_hooks", envir = ns)
say("current instance tracks hooks for: %s\n", paste(ls(he), collapse=","))
q(pkgload::unload("frmtmb"))
for (p in c("lme4","bayesplot")) { v <- mine(p); say("after-unload %-10s total=%d frmtmb-shaped=%d\n", p, v[1], v[2]) }
# and does a LATER load of lme4 now run N stale hooks pointing into
# dead namespaces?
r <- tryCatch({ q(loadNamespace("lme4")); "ok" }, error = function(e) conditionMessage(e))
say("loadNamespace(lme4) after unload: %s\n", r)
say("GENREVDONE\n")
