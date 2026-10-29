## Review rescorpp: run a test file with one core function mutated in
## the loaded namespaces (no install).
##   Rscript dev/review-rescorpp-mutate.R <mutation> <pkg> <file>
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
ns <- asNamespace("frmtmb")
patch <- function(name, fn) {
  environment(fn) <- ns
  utils::assignInNamespace(name, fn, "frmtmb")
  for (p in c("frmtmb.sample")) {
    imp <- parent.env(asNamespace(p))
    if (exists(name, envir = imp, inherits = FALSE)) {
      unlockBinding(name, imp); assign(name, fn, envir = imp)
      lockBinding(name, imp)
    }
  }
  if ("package:frmtmb" %in% search() &&
      exists(name, envir = as.environment("package:frmtmb"), inherits = FALSE)) {
    e <- as.environment("package:frmtmb")
    unlockBinding(name, e); assign(name, fn, envir = e); lockBinding(name, e)
  }
}
src <- get("predict_rescor_draw", ns)
txt <- deparse(src)
if (mutation == "nochol") {
  # independent draws: drop the Cholesky factor
  txt <- sub("L <- chol(R)", "L <- diag(nrow(R))", txt, fixed = TRUE)
  patch("predict_rescor_draw", eval(parse(text = txt)))
} else if (mutation == "sepW") {
  # one chi-square per response instead of one per row
  txt <- sub("Z <- Z/sqrt(stats::rchisq(m, nu)/nu)",
             "Z <- Z/sqrt(matrix(stats::rchisq(length(Z), nu), m)/nu)",
             txt, fixed = TRUE)
  patch("predict_rescor_draw", eval(parse(text = txt)))
} else if (mutation == "sortC") {
  # the subset's correlation in model order, not the order resp names
  t2 <- deparse(get("rescor_row_loglik", ns))
  t2 <- sub("C <- C[k, k, drop = FALSE]", "C <- C[sort(k), sort(k), drop = FALSE]",
            t2, fixed = TRUE)
  if (!any(grepl("sort(k)", t2, fixed = TRUE))) stop("sortC did not apply")
  patch("rescor_row_loglik", eval(parse(text = t2)))
  txt <- character(0)
} else if (mutation == "nonu") {
  # a student rescor draw as if gaussian
  txt <- sub("nu <- rescor_shared_nu(fs)", "nu <- NULL", txt, fixed = TRUE)
  patch("predict_rescor_draw", eval(parse(text = txt)))
}
if (length(txt) && identical(deparse(get("predict_rescor_draw", ns)), deparse(src))) {
  stop("mutation ", mutation, " did not apply")
}
cat("mutation", mutation, "applied\n")
tdir <- if (pkg == "frmtmb") file.path(wt, "tests/testthat") else
  file.path(wt, "extensions", pkg, "tests/testthat")
res <- as.data.frame(testthat::test_file(
  file.path(tdir, tf), package = pkg, env = testthat::test_env(pkg),
  reporter = "check", stop_on_failure = FALSE))
print(res[res$failed > 0 | res$error, c("test", "failed", "error")])
cat("RESULT", mutation, tf, "pass", sum(res$passed), "fail", sum(res$failed),
    "error", sum(res$error), "skip", sum(res$skipped), "\n")
