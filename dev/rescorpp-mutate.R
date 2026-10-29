## Lane rescorpp, fix round 1: run a test file with one function mutated
## in the loaded namespaces (no install), after the reviewer's
## dev/review-rescorpp-mutate.R.
##   Rscript dev/rescorpp-mutate.R <mutation> <pkg> <file>
a <- commandArgs(TRUE)
mutation <- a[1]; pkg <- a[2]; tf <- a[3]
.libPaths(c("C:/Users/adf44/source/r/wt-rescorpp-lib",
            "C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-rescorpp"
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = file.path(wt, "dev/stan-cache"))
suppressMessages(library(pkg, character.only = TRUE))
suppressMessages(library(frmtmb.sample))
patch <- function(name, from, old, new, owner = "frmtmb") {
  ns <- asNamespace(owner)
  txt <- deparse(get(name, ns))
  hit <- grep(old, txt, fixed = TRUE)
  if (length(hit) != 1L) stop("mutation ", mutation, " did not apply")
  txt[hit] <- sub(old, new, txt[hit], fixed = TRUE)
  fn <- eval(parse(text = txt))
  environment(fn) <- ns
  utils::assignInNamespace(name, fn, owner)
  if (owner == "frmtmb") {
    imp <- parent.env(asNamespace("frmtmb.sample"))
    if (exists(name, envir = imp, inherits = FALSE)) {
      unlockBinding(name, imp); assign(name, fn, envir = imp)
      lockBinding(name, imp)
    }
    e <- if ("package:frmtmb" %in% search()) {
      as.environment("package:frmtmb")
    } else emptyenv()
    if (exists(name, envir = e, inherits = FALSE)) {
      unlockBinding(name, e); assign(name, fn, envir = e); lockBinding(name, e)
    }
  } else {
    S3 <- get(".__S3MethodsTable__.", envir = ns)
    if (exists(name, envir = S3, inherits = FALSE)) assign(name, fn, envir = S3)
  }
}
switch(mutation,
  # a student rescor draw as if gaussian (the reviewer's "nonu")
  nonu = patch("predict_rescor_draw", NULL, "nu <- rescor_shared_nu(fs)",
               "nu <- NULL"),
  # the correlation of a resp subset in model order (the reviewer's sortC)
  sortC = patch("rescor_row_loglik", NULL, "C <- C[k, k, drop = FALSE]",
                "C <- C[sort(k), sort(k), drop = FALSE]"),
  # the draw's correlation in model order, whatever order rspecs has
  sortR = patch("predict_rescor_draw", NULL,
                "R <- rescor_matrix(fs)[names(rspecs), names(rspecs), drop = FALSE]",
                "R <- rescor_matrix(fs)[sort(names(rspecs)), sort(names(rspecs)), drop = FALSE]"),
  # dps paired by position, as before the fix
  pospair = patch("predict_rescor_draw", NULL, "dps <- dps[names(rspecs)]",
                  "dps <- dps"),
  # predict() without the category refusal, as the lane had it before
  nocatref = patch("predict.frmtmb_draws", NULL,
                   "if (summary && draws_resp_stacks(fit, resp)) {",
                   "if (FALSE) {", owner = "frmtmb.sample"),
  # a repeated name passed through, as the lane had it before
  nodup = patch("draws_resp_set", NULL, "if (anyDuplicated(hit)) {",
                "if (FALSE) {", owner = "frmtmb.sample"),
  stop("unknown mutation"))
cat("mutation", mutation, "applied\n")
tdir <- if (pkg == "frmtmb") file.path(wt, "tests/testthat") else
  file.path(wt, "extensions", pkg, "tests/testthat")
res <- as.data.frame(testthat::test_file(
  file.path(tdir, tf), package = pkg, env = testthat::test_env(pkg),
  reporter = "silent", stop_on_failure = FALSE))
print(res[res$failed > 0 | res$error, c("test", "failed", "error")])
cat("RESULT", mutation, tf, "pass", sum(res$passed), "fail", sum(res$failed),
    "error", sum(res$error), "skip", sum(res$skipped), "\n")
