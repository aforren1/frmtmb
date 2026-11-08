# Review of lane latperf: one test file in one process, optionally under an
# in-process mutant of R/joint-sparse.R.
# Usage: Rscript dev/review-latperf-test.R <test-file> [mutant]
# mutants:
#   none      the lane build as installed
#   novf      jcs_rowvar() drops the outer part g' Vf g
#   negj      jcs_joint_cov() stores -J in place of J
#   perm      the selected inverse in the factor's permuted order
#   norepair  jcs_joint_cov() ignores the SE check's propagation covariance
#   nojvj     `[` drops J Vf J' from the inner block
#   diagonly  jcs_rowvar_inner() reads only the diagonal pairs of W
local({
  .libPaths(c("C:/Users/adf44/source/r/wt-latperf-lib",
              "C:/Users/adf44/source/r/rellib-r18",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  Sys.setenv(NOT_CRAN = "true")
  a <- commandArgs(TRUE)
  tf <- a[1]
  mut <- if (length(a) > 1) a[2] else "none"
  suppressMessages(library(testthat))
  suppressMessages(library(frmtmb))
  cat("lib:", dirname(find.package("frmtmb")), "frmtmb",
      format(packageVersion("frmtmb")), "\n")
  cat("mutant:", mut, "file:", tf, "\n")
  ns <- asNamespace("frmtmb")
  put <- function(nm, fn) {
    environment(fn) <- ns
    utils::assignInNamespace(nm, fn, ns = "frmtmb")
  }
  src_sub <- function(nm, from, to) {
    fn <- get(nm, envir = ns)
    txt <- deparse(fn, control = "useSource")
    hit <- grepl(from, txt, fixed = TRUE)
    if (sum(hit) != 1L) stop("mutant site not unique in ", nm, ": ", from)
    txt[hit] <- sub(from, to, txt[hit], fixed = TRUE)
    put(nm, eval(parse(text = txt)))
  }
  switch(mut,
    none = NULL,
    novf = src_sub("jcs_rowvar",
                   "v <- rowSums((G %*% S[[\"Vf\"]]) * G)",
                   "v <- 0 * rowSums((G %*% S[[\"Vf\"]]) * G)"),
    negj = src_sub("jcs_joint_cov", "J <- -as.matrix(", "J <- as.matrix("),
    perm = put("jcs_selinv", function(S) {
      .Call("tmb_invQ", S[["L"]], PACKAGE = "TMB")
    }),
    norepair = src_sub("jcs_joint_cov", "repair <- !is.null(Vp) &&",
                       "repair <- FALSE && !is.null(Vp) &&"),
    nojvj = src_sub("[.frm_joint_sparse",
                    "out[a, b] <- jcs_W(x, ri[a], rj[b]) +",
                    "out[a, b] <- jcs_W(x, ri[a], rj[b]) + 0 *"),
    diagonly = src_sub("jcs_rowvar_inner",
                       "e2 <- sequence(mi, from = cumsum(c(1, m))[row])",
                       "e2 <- e1"),
    stop("unknown mutant ", mut))
  if (mut == "diagonly") {
    # e1 repeats each entry mi times; the diagonal alone is one each
    src_sub("jcs_rowvar_inner", "e1 <- rep(seq_along(row), times = mi)",
            "e1 <- seq_along(row)")
  }
  path <- file.path("C:/Users/adf44/source/r/frmtmb-wt-latperf/tests/testthat",
                    tf)
  t0 <- proc.time()[["elapsed"]]
  res <- testthat::test_file(path, package = "frmtmb",
                             env = testthat::test_env("frmtmb"),
                             reporter = "silent", stop_on_failure = FALSE)
  df <- as.data.frame(res)
  cat(sprintf(paste("RESULT mutant %s file %s tests %d pass %d fail %d",
                    "err %d skip %d warn %d time %.1f\n"),
              mut, tf, nrow(df), sum(df$passed), sum(df$failed),
              sum(df$error), sum(df$skipped), sum(df$warning),
              proc.time()[["elapsed"]] - t0))
  bad <- df[df$failed > 0 | df$error, c("test", "failed", "error")]
  if (nrow(bad)) {
    for (i in seq_len(nrow(bad))) {
      cat("  FAILED:", bad$test[i], "\n")
    }
  }
})
