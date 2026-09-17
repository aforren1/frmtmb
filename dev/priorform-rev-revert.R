# Reviewer, lane wt-priorform: put back ONE piece of the base behavior
# in a lane process and run the lane's test files, to see whether the
# tests that guard that piece fail without it.
#   Rscript dev/priorform-rev-revert.R <revert> <filter>
# revert: none dup cs nested bf bfnl lf asprior setprior print blockdpar
args <- commandArgs(trailingOnly = TRUE)
rev <- args[1]; filt <- args[2]
.libPaths(c("C:/Users/adf44/source/r/priorform-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
root <- "C:/Users/adf44/source/r/frmtmb-wt-priorform"
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
base <- readRDS(file.path(root, "dev", "priorform-rev-revert-base.rds"))
put <- function(nm, f) {
  environment(f) <- ns
  utils::assignInNamespace(nm, f, ns = "frmtmb")
  if (nm %in% getNamespaceExports("frmtmb")) {
    # the attached copy is what a test's bare call reaches
    env <- as.environment("package:frmtmb")
    unlockBinding(nm, env); assign(nm, f, envir = env); lockBinding(nm, env)
  }
}
switch(rev,
  none = NULL,
  dup = put("refuse_duplicated_re", function(cps) invisible(NULL)),
  cs = put("parse_linpred", base$parse_linpred),
  nested = put("refuse_nested_formula", function(f) invisible(f)),
  bf = put("bf", base$bf),
  bfnl = {
    txt <- deparse(get("bf", ns))
    hit <- grepl("bf(f1, ..., family = family, nl = nl)", txt, fixed = TRUE)
    stopifnot(sum(hit) == 1L)
    txt <- sub("bf(f1, ..., family = family, nl = nl)",
               "bf(f1, ..., family = family)", txt, fixed = TRUE)
    put("bf", eval(parse(text = txt)))
  },
  lf = put("lf", base$lf),
  asprior = put("as_priorlist", base$as_priorlist),
  setprior = put("set_prior", base$set_prior),
  print = {
    f <- base$print.frmtmb_priorlist; environment(f) <- ns
    registerS3method("print", "frmtmb_priorlist", f, envir = ns)
  },
  blockdpar = put("block_dpar", function(spec, frame, bk) ""),
  stop("unknown revert ", rev))
setwd(file.path(root, "tests", "testthat"))
res <- testthat::test_dir(".", filter = filt, package = "frmtmb",
                          reporter = "silent", stop_on_failure = FALSE)
d <- as.data.frame(res)
bad <- d[d$failed > 0 | d$error, c("file", "test")]
cat(sprintf("REVERT %s %s pass %d fail %d error %d\n", rev, filt,
            sum(d$passed), sum(d$failed), sum(d$error)))
if (nrow(bad)) for (i in seq_len(nrow(bad))) {
  cat("   ", bad$file[i], ":", bad$test[i], "\n")
}
