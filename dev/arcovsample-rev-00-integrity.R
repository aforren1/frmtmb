# REVIEW script 00: does the worker's installed library match the
# worktree sources, and are the run-tests.R assertions true?
#
#   Rscript dev/arcovsample-rev-00-integrity.R

LIB <- "C:/Users/adf44/source/r/wt-arcovsample-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
cat("tmbstan  ", format(utils::packageVersion("tmbstan")), "\n")
cat("makevars ", tools::makevars_user(), "\n")
cat("gnu++17  ", any(grepl("-std=gnu++17",
                           readLines(tools::makevars_user(), warn = FALSE),
                           fixed = TRUE)), "\n")
cat("rstan    ", format(utils::packageVersion("rstan")), "\n")
cat("StanHead ", format(utils::packageVersion("StanHeaders")), "\n")
cat("brms     ", format(utils::packageVersion("brms")), "\n")
cat("mvtnorm  ", format(utils::packageVersion("mvtnorm")), "\n")

suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
cat("\nlib paths:\n"); print(.libPaths())
cat("frmtmb        ", format(packageVersion("frmtmb")), " at ",
    dirname(system.file("DESCRIPTION", package = "frmtmb")), "\n")
cat("frmtmb.sample ", format(packageVersion("frmtmb.sample")), " at ",
    dirname(system.file("DESCRIPTION", package = "frmtmb.sample")), "\n")

# The installed object versus the worktree source, for every function
# the diff touched. deparse() of the installed closure against a
# deparse() of the source expression parsed out of the file: a stale
# install would show here rather than in a test that happens to pass.
src_fun <- function(path, name) {
  ex <- parse(file = path)
  for (e in ex) {
    if (is.call(e) && length(e) >= 3L &&
        as.character(e[[1L]])[1L] %in% c("<-", "=") &&
        identical(as.character(e[[2L]]), name)) {
      return(eval(e[[3L]], envir = new.env()))
    }
  }
  NULL
}

cmp <- function(pkg, path, name) {
  ns <- asNamespace(pkg)
  inst <- if (exists(name, envir = ns, inherits = FALSE)) {
    get(name, envir = ns)
  } else NULL
  s <- src_fun(path, name)
  if (is.null(inst)) { cat("  MISSING in install: ", pkg, "::", name, "\n",
                          sep = ""); return(invisible(FALSE)) }
  if (is.null(s)) { cat("  NOT FOUND in source: ", path, " ", name, "\n",
                        sep = ""); return(invisible(FALSE)) }
  ok <- identical(deparse(inst), deparse(s))
  cat("  ", if (ok) "same" else "DIFFERS", "  ", pkg, "::", name, "\n",
      sep = "")
  invisible(ok)
}

WT <- "C:/Users/adf44/source/r/frmtmb-wt-arcovsample"
cat("\ninstalled body vs worktree source:\n")
ok <- c(
  cmp("frmtmb", file.path(WT, "R/autocor.R"), "arma_cond_resp"),
  cmp("frmtmb", file.path(WT, "R/autocor.R"), "arma_cond_dpars"),
  cmp("frmtmb", file.path(WT, "R/predict-brms.R"), "rescor_row_loglik"),
  cmp("frmtmb.sample", file.path(WT, "extensions/frmtmb.sample/R/loo.R"),
      "draws_chain_id"),
  cmp("frmtmb.sample", file.path(WT, "extensions/frmtmb.sample/R/loo.R"),
      "draws_loglik_factors"),
  cmp("frmtmb.sample", file.path(WT, "extensions/frmtmb.sample/R/loo.R"),
      "draws_row_loglik"),
  cmp("frmtmb.sample",
      file.path(WT, "extensions/frmtmb.sample/R/methods-draws.R"),
      "posterior_predict.frmtmb_draws"))
cat("all same: ", all(ok), "\n")

cat("\nexported from frmtmb: ",
    all(c("arma_cond_resp", "arma_cond_dpars") %in%
          getNamespaceExports("frmtmb")), "\n")

# The reference build must NOT have them: the floor claim rests on it.
cat("\nreference build (rellib-r3):\n")
rl <- "C:/Users/adf44/source/r/rellib-r3"
nsr <- loadNamespace("frmtmb", lib.loc = rl)
cat("  frmtmb ", format(getNamespaceVersion(nsr)), "\n")
cat("  has arma_cond_resp : ",
    "arma_cond_resp" %in% getNamespaceExports(nsr), "\n")
cat("  has arma_cond_dpars: ",
    "arma_cond_dpars" %in% getNamespaceExports(nsr), "\n")
cat("DONE\n")
