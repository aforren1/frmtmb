# genrev: an INDEPENDENT reimplementation of the load-order probe.
#
# Differences from dev/generics-check.R on purpose:
#  * it reproduces UseMethod's lookup DIRECTLY (the .__S3MethodsTable__.
#    of environment(generic)) instead of inferring it from an error
#    message, so "dispatched to somebody else's .default" is visible as
#    a distinct outcome instead of counting as a success;
#  * it records WHICH package owns the method that would run;
#  * it covers load orders the lane did not run.
#
# Usage: Rscript genrev-probe.R <LIB> <mode>
args <- commandArgs(trailingOnly = TRUE)
LIB <- args[1]
mode <- args[2]
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

GENS <- c("as_draws", "as_draws_array", "as_draws_df", "as_draws_list",
          "as_draws_matrix", "as_draws_rvars", "nchains", "ndraws",
          "niterations", "nvariables", "variables", "loo", "loo_compare",
          "waic", "bayes_R2", "prior_summary", "pp_check",
          "conditional_effects", "expose_functions", "hypothesis",
          "posterior_summary", "LOO", "WAIC", "ngrps", "fixef", "ranef",
          "VarCorr", "refit")

owner_of <- function(f) {
  if (!is.function(f)) return("<none>")
  e <- environment(f)
  if (is.null(e)) return("base")
  n <- environmentName(topenv(e))
  if (nzchar(n)) n else "<anon>"
}

# Reproduce what UseMethod() does for a generic reached from the global
# environment: look in the caller's top environment, then in the S3
# methods table of the namespace where the generic was DEFINED.
lookup <- function(g, cls) {
  f <- tryCatch(get(g, envir = globalenv()), error = function(e) NULL)
  if (!is.function(f)) {
    return(list(gen = "<absent>", meth = "<absent>", via = "<absent>"))
  }
  genv <- environment(f)
  gen <- owner_of(f)
  tbl <- NULL
  if (!is.null(genv) && isNamespace(genv)) {
    tbl <- tryCatch(get(".__S3MethodsTable__.", envir = genv),
                    error = function(e) NULL)
  }
  pick <- function(nm) {
    v <- NULL
    if (!is.null(tbl) && exists(nm, envir = tbl, inherits = FALSE)) {
      v <- get(nm, envir = tbl)
    } else if (!is.null(genv) && exists(nm, envir = genv, inherits = FALSE)) {
      v <- get(nm, envir = genv)
    } else if (exists(nm, envir = globalenv())) {
      v <- get(nm, envir = globalenv())
    }
    v
  }
  m <- pick(paste0(g, ".", cls))
  if (!is.null(m)) return(list(gen = gen, meth = owner_of(m), via = "class"))
  d <- pick(paste0(g, ".default"))
  if (!is.null(d)) return(list(gen = gen, meth = owner_of(d), via = "default"))
  list(gen = gen, meth = "<none>", via = "none")
}

# getS3method() is the other column the lane reports; it is a different
# code path from UseMethod and is what tooling (broom, emmeans, etc.)
# uses to find a method.
gs3 <- function(g, cls) {
  m <- tryCatch(utils::getS3method(g, cls, optional = TRUE),
                error = function(e) NULL)
  if (is.null(m)) "<none>" else owner_of(m)
}

setup <- function(mode) {
  q <- function(...) suppressMessages(suppressWarnings(...))
  if (mode == "ctrlB")      { q(library(brms)) }
  else if (mode == "ctrlL") { q(library(lme4)) }
  else if (mode == "A")     { q(library(brms)); q(library(frmtmb)) }
  else if (mode == "B")     { q(library(frmtmb)); q(library(brms)) }
  else if (mode == "D")     { q(library(frmtmb)); q(loadNamespace("brms")) }
  else if (mode == "C")     { q(library(frmtmb)) }
  else if (mode == "L")     { q(library(lme4)); q(library(frmtmb)) }
  else if (mode == "THIRD") {
    # frmtmb LOADED but never ATTACHED, reached through a third package
    # that Imports it.  brms attached first, as in mode A.
    q(library(brms)); q(library(genrevdep))
  }
  else if (mode == "THIRD2") {
    # the same, with brms only loadNamespace()d
    q(library(genrevdep)); q(loadNamespace("brms"))
  }
  else if (mode == "UNLOADB") {
    # both attached, then brms taken away and put back
    q(library(brms)); q(library(frmtmb))
    detach("package:brms", unload = FALSE)
    q(unloadNamespace("brms"))
    q(loadNamespace("brms"))
  }
  else if (mode == "RELOADOWNER") {
    # frmtmb attached, an owner loaded (so no hook is left for it),
    # then that owner unloaded and RELOADED
    q(library(frmtmb)); q(loadNamespace("loo"))
    q(unloadNamespace("loo"))
    q(loadNamespace("loo"))
  }
  else if (mode == "REATTACH") {
    q(library(frmtmb)); detach("package:frmtmb", unload = FALSE)
    q(loadNamespace("brms")); q(library(frmtmb))
  }
  else if (mode == "NSCYCLE") {
    q(loadNamespace("frmtmb")); q(unloadNamespace("frmtmb"))
    q(library(brms)); q(library(frmtmb))
  }
  else stop("unknown mode ", mode)
}

setup(mode)
cls <- if (mode %in% c("ctrlL", "L")) "merMod" else "brmsfit"
cat(sprintf("GENREVMODE %s CLASS %s\n", mode, cls))
cat(sprintf("SEARCH %s\n", paste(search(), collapse = " ")))
for (g in GENS) {
  r <- lookup(g, cls)
  cat(sprintf("ROW\t%s\t%s\t%s\t%s\t%s\n", g, r$gen, r$meth, r$via,
              gs3(g, cls)))
}
cat("GENREVDONE\n")
