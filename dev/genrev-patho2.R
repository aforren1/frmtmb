args <- commandArgs(trailingOnly = TRUE)
LIB <- args[1]; case <- args[2]
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(...) suppressMessages(suppressWarnings(...))
own <- function(f) { if (!is.function(f)) return("<absent>")
  e <- environment(f); if (is.null(e)) return("base")
  n <- environmentName(topenv(e)); if (nzchar(n)) n else "<anon>" }
reach <- function(g, cls) {
  f <- tryCatch(get(g, envir = globalenv()), error = function(e) NULL)
  if (!is.function(f)) return(c("<absent>", "<absent>", NA))
  genv <- environment(f)
  tbl <- tryCatch(get(".__S3MethodsTable__.", envir = genv), error=function(e) NULL)
  nm <- paste0(g, ".", cls)
  m <- if (!is.null(tbl) && exists(nm, envir = tbl, inherits = FALSE)) get(nm, envir=tbl) else NULL
  n <- environmentName(topenv(genv))
  st <- if (!nzchar(n) || !isNamespaceLoaded(n)) TRUE else !identical(topenv(genv), asNamespace(n))
  c(own(f), if (is.null(m)) "<none>" else own(m), st)
}
say <- function(...) cat(sprintf(...))

if (case == "stale") {
  # The owner is loaded BEFORE frmtmb, so frm_adopt_generics() leaves
  # NO hook for it.  Unload it, then let something else load it fresh.
  q(library(loo)); q(library(frmtmb))
  say("frmtmb ns loo -> %s\n", own(get("loo", envir = asNamespace("frmtmb"))))
  say("hook for loo: %d\n", length(getHook(packageEvent("loo","onLoad"))))
  detach("package:loo", unload = FALSE)
  r <- tryCatch({ unloadNamespace("loo"); "ok" }, error = function(e) conditionMessage(e))
  say("unload loo: %s\n", r)
  q(loadNamespace("brms"))
  say("loo re-loaded fresh: %s\n", isNamespaceLoaded("loo"))
  for (g in c("loo","waic","loo_compare")) {
    v <- reach(g, "brmsfit")
    say("  %-12s generic=%-10s method=%-10s STALE=%s\n", g, v[1], v[2], v[3])
  }
  say("frmtmb ns loo identical to live loo::loo: %s\n",
      identical(get("loo", envir = asNamespace("frmtmb")), loo::loo))
}

if (case == "sampleprop") {
  # The order dev/generics-propagate.R did not run: frmtmb.sample
  # attached, then brms LOADED but never attached.
  q(library(frmtmb.sample))
  gens <- c("loo","ngrps","posterior_summary","as_draws_df","pp_check",
            "VarCorr","fixef","variables","waic","bayes_R2","hypothesis",
            "conditional_effects","prior_summary","nchains","ndraws",
            "niterations","nvariables","as_draws","as_draws_array",
            "as_draws_list","as_draws_matrix","as_draws_rvars","LOO",
            "WAIC","expose_functions","ranef")
  q(loadNamespace("brms"))
  lost <- character()
  for (g in gens) { v <- reach(g, "brmsfit"); if (v[2] == "<none>") lost <- c(lost, g) }
  say("MODE frmtmb.sample attached, brms loaded only\n")
  say("re-exported-from-core lost %d of %d\n", length(lost), length(gens))
  say("lost: %s\n", paste(lost, collapse=", "))
  # and the same for the names frmtmb.sample defines itself
  say("own-defined log_lik generic owner: %s\n", reach("log_lik","brmsfit")[1])
}

if (case == "loadall") {
  if (!requireNamespace("pkgload", quietly = TRUE)) { say("pkgload absent\n"); quit(save="no") }
  q(library(brms))
  hookn <- function(p) length(getHook(packageEvent(p, "onLoad")))
  say("hooks-before lme4=%d bayesplot=%d\n", hookn("lme4"), hookn("bayesplot"))
  for (i in 1:3) {
    ok <- tryCatch({ q(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-generics",
            export_all = FALSE, quiet = TRUE, helpers = FALSE)); "ok" },
            error = function(e) conditionMessage(e))
    v <- reach("loo","brmsfit"); w <- reach("posterior_summary","brmsfit")
    say("load_all %d: %s hooks lme4=%d bayesplot=%d loo=%s/%s postsum=%s/%s\n",
        i, ok, hookn("lme4"), hookn("bayesplot"), v[1], v[2], w[1], w[2])
  }
}
say("GENREVDONE\n")
