# Reviewer of lane surface, claim 9: does sd2_floored() change any fit
# a test file makes? One process, two arms (pitfall 21: no cross-process
# ulp noise): run the test file with the floor, then with sd2_floored()
# replaced by the 0.68.1 expression exp(2 * log_sd), recording every
# fit_assembled() result, and compare the two arms with identical().
#   Rscript dev/surface-rev-floorab.R <test file> [core|sample]
a <- commandArgs(TRUE)
source("dev/surface-rev-env.R")
rev_env("lane")
Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")
pkg <- if (length(a) > 1 && a[2] == "sample") "frmtmb.sample" else "frmtmb"
suppressPackageStartupMessages({
  library(testthat); library(frmtmb)
  if (pkg == "frmtmb.sample") library(frmtmb.sample)
})
tdir <- if (pkg == "frmtmb") "tests/testthat" else
  "extensions/frmtmb.sample/tests/testthat"
ns <- asNamespace("frmtmb")
orig_fa <- get("fit_assembled", ns)
rec <- list()
wrap <- function(...) {
  out <- orig_fa(...)
  cs <- vapply(out$frame$re_blocks %||% list(), function(b) {
    b[["covstruct"]] %||% ""
  }, "")
  dims <- vapply(out$frame$re_blocks %||% list(), function(b) {
    as.integer(b[["dim"]] %||% NA)
  }, 1L)
  rec[[length(rec) + 1L]] <<- list(
    cs = cs, dims = dims, est = out$estimates,
    ll = tryCatch(as.numeric(logLik(out)), error = function(e) NA),
    v = tryCatch(vcov(out, full = TRUE), error = function(e) NULL))
  out
}
`%||%` <- function(x, y) if (is.null(x)) y else x
assignInNamespace("fit_assembled", wrap, "frmtmb")
run <- function() {
  rec <<- list()
  res <- as.data.frame(test_file(file.path(tdir, a[1]), package = pkg,
                                 env = test_env(pkg), reporter = "silent",
                                 stop_on_failure = FALSE))
  bad <- res[res$failed > 0 | res$error, "test"]
  if (length(bad)) cat("  failing:", paste(bad, collapse = " | "), "\n")
  list(rec = rec, fail = sum(res$failed), err = sum(res$error),
       pass = sum(res$passed))
}
set.seed(1)
A <- run()
assignInNamespace("sd2_floored", function(log_sd) exp(2 * log_sd), "frmtmb")
set.seed(1)
B <- run()
# a second run of a file can stop early on state the first left behind
# (test-data2.R); compare the fits both arms made, in order
nA <- length(A$rec); nB <- length(B$rec)
if (nA != nB) {
  # pair the fits by structure, parameter names and logLik to 1e-6
  # relative, rather than by position
  key <- function(r) paste(paste(r$cs, collapse = ","),
                           paste(names(unlist(r$est)), collapse = ","),
                           signif(r$ll, 7))
  kA <- vapply(A$rec, key, ""); kB <- vapply(B$rec, key, "")
  m <- match(kA, kB)
  cat("  arms differ in length: paired", sum(!is.na(m)), "of", nA,
      "fits by key\n")
  A$rec <- A$rec[!is.na(m)]; B$rec <- B$rec[m[!is.na(m)]]
}
one <- c("gp", "ou", "homcs", "homtoep", "exp", "gau", "mat", "gr_cov")
touch <- vapply(A$rec, function(r) {
  any(r$cs %in% setdiff(one, "gr_cov")) ||
    any(r$cs == "gr_cov" & r$dims == 1L)
}, NA)
same <- mapply(function(x, y) {
  identical(x$est, y$est) && identical(x$ll, y$ll) && identical(x$v, y$v)
}, A$rec, B$rec)
if (!length(same)) same <- logical()
cat(sprintf(paste("FLOORAB %s fits %d/%d | floor-reachable fits %d |",
                  "identical (est, logLik, vcov) %s | differ %s |",
                  "arms pass %d/%d fail %d/%d err %d/%d\n"),
            a[1], nA, nB, sum(touch),
            if (anyNA(same)) "NA" else sum(same[touch]),
            if (anyNA(same)) "NA" else sum(!same),
            A$pass, B$pass, A$fail, B$fail, A$err, B$err))
cat("structures seen:",
    paste(names(table(unlist(lapply(A$rec[touch], `[[`, "cs")))),
          table(unlist(lapply(A$rec[touch], `[[`, "cs"))), sep = "=",
          collapse = " "), "\n")
