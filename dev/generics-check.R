# Phase 2.5a/b/c measurement.  One process per load order; a single
# process cannot see this defect, because the search path is built once.
#
# Usage: Rscript dev/generics-check.R <LIB> <mode>
#   mode X  library(brms) alone               the CONTROL
#   mode A  library(brms); library(frmtmb)    the reported order
#   mode B  library(frmtmb); library(brms)    the reverse
#   mode C  library(frmtmb) alone             every owner absent
#   mode D  library(frmtmb); loadNamespace(brms)  loaded, not attached
#   mode L  library(lme4); library(frmtmb)    the ngrps/refit residue
#   mode M  library(lme4) alone               control for L
#
# THE INSTRUMENT, which the first version of this script got wrong.
#
# The first version decided "did dispatch work" by matching the string
# "no applicable method" in the error.  Two things were wrong with it
# and review found both.
#
#  * It scored a generic as HEALTHY when the call fell silently into
#    frmtmb's own `.default`.  `posterior_summary` and `loo_compare`
#    are exactly that case on the base build, so the base damage was
#    reported as 25 of 27 when it is 27 of 27, and the two it missed
#    are the worst two, because a silent wrong answer beats an error.
#  * On a `frmtmb_fit` it scored a REGRESSION as an improvement: the
#    six `as_draws*` refusals went from R's "no applicable method" to
#    posterior's "All list elements must be lists themselves", and the
#    counter, which tested only for the first string, went 1 to 0.
#
# So this version does not read error messages at all.  It reproduces
# UseMethod()'s own lookup: the visible environment first, then the
# `.__S3MethodsTable__.` of `environment(generic)`, for the object's
# class and then for "default".  It reports WHICH package supplies the
# method that would run, so "answered by somebody else's default" is a
# distinct outcome from "answered".
args <- commandArgs(trailingOnly = TRUE)
LIB <- args[1]
mode <- args[2]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

gens <- c("as_draws", "as_draws_array", "as_draws_df", "as_draws_list",
          "as_draws_matrix", "as_draws_rvars", "nchains", "ndraws",
          "niterations", "nvariables", "variables", "loo", "loo_compare",
          "waic", "bayes_R2", "prior_summary", "pp_check",
          "conditional_effects", "expose_functions", "hypothesis",
          "posterior_summary", "LOO", "WAIC", "ngrps", "fixef", "ranef",
          "VarCorr", "refit")

cls <- "brmsfit"
if (mode == "X") {
  suppressMessages(library(brms))
} else if (mode == "A") {
  suppressMessages(library(brms)); suppressMessages(library(frmtmb))
} else if (mode == "B") {
  suppressMessages(library(frmtmb)); suppressMessages(library(brms))
} else if (mode == "C") {
  suppressMessages(library(frmtmb))
} else if (mode == "D") {
  suppressMessages(library(frmtmb)); suppressMessages(loadNamespace("brms"))
} else if (mode == "L") {
  suppressMessages(library(lme4)); suppressMessages(library(frmtmb))
  cls <- "merMod"
} else if (mode == "M") {
  suppressMessages(library(lme4)); cls <- "merMod"
}

where_from <- function(f) {
  if (!is.function(f)) return("-")
  e <- environment(f)
  if (is.null(e)) return("base")
  n <- environmentName(topenv(e))
  if (nzchar(n)) n else "?"
}

# UseMethod()'s lookup, reproduced.  Returns the package that supplies
# the method and whether it matched the CLASS or fell through to
# `default`.
resolve <- function(gen, k) {
  f <- tryCatch(get(gen), error = function(e) NULL)
  if (!is.function(f)) return(c(pkg = "absent", how = "absent"))
  env <- environment(f)
  if (is.null(env)) env <- baseenv()
  tab <- if (exists(".__S3MethodsTable__.", envir = env,
                    inherits = FALSE)) {
    get(".__S3MethodsTable__.", envir = env, inherits = FALSE)
  } else NULL
  for (want in c(k, "default")) {
    nm <- paste0(gen, ".", want)
    v <- tryCatch(get(nm, envir = globalenv()), error = function(e) NULL)
    if (is.function(v)) return(c(pkg = where_from(v), how = want))
    if (!is.null(tab) && exists(nm, envir = tab, inherits = FALSE)) {
      return(c(pkg = where_from(get(nm, envir = tab)), how = want))
    }
  }
  c(pkg = "none", how = "none")
}

cat(sprintf("MODE %s CLASS %s\n", mode, cls))
cat(sprintf("%-20s %-12s %-12s %-10s %-6s\n", "generic",
            "generic-from", "method-from", "matched", "getS3"))
for (g in gens) {
  f <- tryCatch(get(g), error = function(e) NULL)
  r <- resolve(g, cls)
  gs <- !is.null(tryCatch(getS3method(g, cls, optional = TRUE),
                          error = function(e) NULL))
  cat(sprintf("%-20s %-12s %-12s %-10s %-6s\n", g, where_from(f),
              r[["pkg"]], r[["how"]], gs))
}

# frmtmb's OWN dispatch, on a real fit, measured the same way.  The
# question is not "did it error" but "does the method that would run
# belong to frmtmb".  A foreign `.default` walking into a frmtmb_fit is
# a FAILURE here, which is what the old counter could not see.
if ("frmtmb" %in% loadedNamespaces()) {
  set.seed(20260915)
  dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
  fit <- frmtmb::frm(frmtmb::bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  own <- c("fixef", "ranef", "VarCorr", "ngrps", "prior_summary",
           "loo", "waic", "bayes_R2", "LOO", "WAIC", "hypothesis",
           "conditional_effects", "expose_functions", "variables",
           "as_draws", "as_draws_df", "as_draws_array", "as_draws_list",
           "as_draws_matrix", "as_draws_rvars", "ndraws", "nchains",
           "niterations", "nvariables", "pp_check", "refit")
  bad <- character()
  cat("\n-- frmtmb_fit, method that would run --\n")
  for (g in own) {
    r <- resolve(g, "frmtmb_fit")
    ok <- r[["how"]] == "frmtmb_fit" &&
      r[["pkg"]] %in% c("frmtmb", "frmtmb.sample")
    if (!ok) bad <- c(bad, g)
    cat(sprintf("  %-20s %-12s %-12s %s\n", g, r[["pkg"]], r[["how"]],
                if (ok) "ok" else "NOT OURS"))
  }
  cat(sprintf("FITMETHOD_NOT_OURS %d of %d\n", length(bad), length(own)))
  if (length(bad)) {
    cat(sprintf("FITMETHOD_LIST %s\n", paste(bad, collapse = ",")))
  }
  # and the message a user actually sees on the six that regressed
  for (g in c("as_draws", "as_draws_df")) {
    m <- tryCatch({ do.call(g, list(fit)); "NO ERROR" },
                  error = function(e) conditionMessage(e))
    cat(sprintf("MSG %s: %s\n", g, substr(gsub("\n", " ", m), 1, 64)))
  }
  cat(sprintf("REFIT_OK %s\n",
              inherits(tryCatch(refit(fit, dd$y), error = function(e) e),
                       "frmtmb_fit")))
  cat(sprintf("NGRPS %s\n",
              paste(tryCatch(ngrps(fit), error = function(e) "ERR"),
                    collapse = ",")))
  ps <- tryCatch(posterior_summary(cbind(a = rnorm(50))),
                 error = function(e) e)
  cat(sprintf("POSTSUMMARY %s\n",
              if (inherits(ps, "error")) paste("ERR", conditionMessage(ps))
              else paste(colnames(ps), collapse = ",")))
}
cat("DONE\n")
