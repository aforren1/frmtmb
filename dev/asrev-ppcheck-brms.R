## SETTLE IT BY RUNNING IT: does brms's pp_check() accept and HONOR
## `re.form`, and does predictive_interval()?
##
## A real brmsfit is reachable without a fresh Stan compile: the method
## tier's helper builds one from a CACHED compiled program plus
## rstan::sampling(algorithm = "Fixed_param"), which the lane rules say
## is the safe half of the Stan surface on this machine.
##
## Usage: Rscript dev/asrev-ppcheck-brms.R
LIB <- "C:/Users/adf44/source/r/asrev-lib2"
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
            FRMTMB_STAN_CACHE = file.path(getwd(), "dev", "stan-cache"))
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
cat("StanHeaders:", as.character(packageVersion("StanHeaders")), "\n")
cat("rstan      :", as.character(packageVersion("rstan")), "\n")
cat("brms       :", as.character(packageVersion("brms")), "\n")

setwd("tests/testthat")
for (f in c("helper-brms.R", "helper-brms-methods.R")) sys.source(
  f, envir = globalenv())

s <- brms_shape("rC0")
bfit <- s$brmsfit
cat("\nbrmsfit built:", paste(class(bfit), collapse = "/"),
    " ndraws:", brms::ndraws(bfit), "\n")

show <- function(lab, expr) {
  r <- tryCatch(force(expr), error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("  %-40s ERROR: %s\n", lab,
                gsub("\n", " ", conditionMessage(r))))
    return(invisible(NULL))
  }
  cat(sprintf("  %-40s ok\n", lab))
  invisible(r)
}

cat("\n==== brms::pp_check(re.form = ) ====\n")
set.seed(101); p_def <- show("pp_check(ndraws = 5)",
                             brms::pp_check(bfit, ndraws = 5))
set.seed(101); p_ref <- show("pp_check(ndraws = 5, re_formula = NA)",
                             brms::pp_check(bfit, ndraws = 5,
                                            re_formula = NA))
set.seed(101); p_alias <- show("pp_check(ndraws = 5, re.form = NA)",
                               brms::pp_check(bfit, ndraws = 5,
                                              re.form = NA))
set.seed(101); p_typo <- show("pp_check(ndraws = 5, re_frmula = NA)",
                              brms::pp_check(bfit, ndraws = 5,
                                             re_frmula = NA))
if (!is.null(p_def) && !is.null(p_alias)) {
  cat("  re.form = NA  changes the answer vs default : ",
      !isTRUE(all.equal(p_def$data, p_alias$data)), "\n", sep = "")
}
if (!is.null(p_ref) && !is.null(p_alias)) {
  cat("  re.form = NA  equals re_formula = NA        : ",
      isTRUE(all.equal(p_ref$data, p_alias$data)), "\n", sep = "")
}
if (!is.null(p_def) && !is.null(p_typo)) {
  cat("  a TYPO (re_frmula) changes the answer       : ",
      !isTRUE(all.equal(p_def$data, p_typo$data)),
      "   <- FALSE means brms swallows unknown names\n", sep = "")
}

cat("\n==== brms::predictive_interval(re.form = ) ====\n")
set.seed(202); i_def <- show("predictive_interval()",
                             brms::predictive_interval(bfit))
set.seed(202); i_ref <- show("predictive_interval(re_formula = NA)",
                             brms::predictive_interval(bfit,
                                                       re_formula = NA))
set.seed(202); i_alias <- show("predictive_interval(re.form = NA)",
                               brms::predictive_interval(bfit,
                                                         re.form = NA))
if (!is.null(i_def) && !is.null(i_alias)) {
  cat("  re.form = NA  changes the answer vs default : ",
      !isTRUE(all.equal(unname(i_def), unname(i_alias))), "\n", sep = "")
}
if (!is.null(i_ref) && !is.null(i_alias)) {
  cat("  re.form = NA  equals re_formula = NA        : ",
      isTRUE(all.equal(unname(i_ref), unname(i_alias))), "\n", sep = "")
  cat("  max |re.form - re_formula|                  : ",
      format(max(abs(i_ref - i_alias)), digits = 12), "\n", sep = "")
}

cat("\n==== the static half, for the record ====\n")
b <- body(getFromNamespace("pp_check.brmsfit", "brms"))
d <- deparse(b)
cat(paste0("  ", grep("pred_args|do_call|yrep|method <-", d,
                      value = TRUE)[1:12]), sep = "\n")
cat("\nDONE\n")
