# Fix round 1, item 6: mutants on the final source, written by me, not the lane's.
# make <m>: source copy under rev-impsamp-src/fr1-<m>. test <m> [desc]: run
# test-importance.R, or one named test, stopping if the mutant's library did not load.
WT <- "C:/Users/adf44/source/r/frmtmb-wt-impsamp"
SRC <- "C:/Users/adf44/source/r/rev-impsamp-src"
LIB <- "C:/Users/adf44/source/r/rev-impsamp-lib"
a <- commandArgs(TRUE)
cmd <- a[[1]]
mut <- a[[2]]
muts <- list(
  partial_rate = list(
    list("      up[rnd] <- up[rnd] - hsolve(d) * h\n", ""),
    list("      up[rnd] <- up[rnd] + du\n      dn[rnd] <- dn[rnd] - du\n", "")),
  half_jacobian = list(
    list("    list(u = u, logdet = logdet, prior = tm)\n",
         "    list(u = u, logdet = 0.5 * (logdet + plan[[\"logdet_l\"]]), prior = tm)\n")),
  capped_settles = list(
    list("settled = !rose && !capped", "settled = !rose")),
  abs_step = list(
    list("    h[pn == cp][known] <- imp_cross_step * s[known]\n", "")),
  no_units = list(
    list("    u <- if (key %in% names(units)) units[[key]] else 1\n", "    u <- 1\n")),
  prior_twice = list(
    list("if (is.null(neg_log_prior)) val else val + neg_log_prior(pars)",
         "if (is.null(neg_log_prior)) val else val + 2 * neg_log_prior(pars)")),
  # a proposed fix: an ordinal threshold moves the predictor one for one
  fix_thres = list(
    list("  h <- imp_cross_step * pmax(1, abs(unname(at)))\n",
         paste0("  h <- imp_cross_step * pmax(1, abs(unname(at)))\n",
                "  h[names(at) %in% c(\"thres\", \"tau\", \"thres_first\", \"thres_inc\", \"thres_par\")] <- imp_cross_step\n")))
)
if (cmd == "make") {
  dst <- file.path(SRC, paste0("fr1-", mut))
  unlink(dst, recursive = TRUE)
  dir.create(file.path(dst, "R"), recursive = TRUE)
  dir.create(file.path(LIB, mut), recursive = TRUE, showWarnings = FALSE)
  for (f in c("DESCRIPTION", "NAMESPACE")) file.copy(file.path(WT, f), dst)
  file.copy(list.files(file.path(WT, "R"), full.names = TRUE), file.path(dst, "R"))
  for (d in c("inst", "data", "man")) file.copy(file.path(WT, d), dst, recursive = TRUE)
  p <- file.path(dst, "R", "importance.R")
  s <- paste0(paste(readLines(p), collapse = "\n"), "\n")
  for (ed in muts[[mut]]) {
    n <- lengths(gregexpr(ed[[1]], s, fixed = TRUE))
    if (!grepl(ed[[1]], s, fixed = TRUE) || n != 1L) stop("site not once: ", ed[[1]])
    s <- sub(ed[[1]], ed[[2]], s, fixed = TRUE)
  }
  writeBin(charToRaw(s), p)
  cat("made", mut, "\n")
} else {
  Sys.setenv(NOT_CRAN = "true")
  .libPaths(c(file.path(LIB, mut), "C:/Users/adf44/source/r/wt-impsamp-lib",
              "C:/Users/adf44/source/r/rellib-r16",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  suppressPackageStartupMessages({
    library(frmtmb)
    library(testthat)
  })
  if (!startsWith(find.package("frmtmb"), LIB)) stop("not the mutant's library")
  f <- file.path(WT, "tests/testthat/test-importance.R")
  res <- if (length(a) >= 3L) {
    test_file(f, package = "frmtmb", env = test_env("frmtmb"), reporter = "silent",
              stop_on_failure = FALSE, desc = a[[3]])
  } else {
    test_file(f, package = "frmtmb", env = test_env("frmtmb"), reporter = "silent",
              stop_on_failure = FALSE)
  }
  r <- as.data.frame(res)
  cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n", mut, sum(r$passed),
              sum(r$failed), sum(r$error), sum(r$skipped), sum(r$warning)))
  bad <- r[r$failed > 0 | r$error, c("test", "failed", "error")]
  for (i in seq_len(nrow(bad))) {
    cat(sprintf("  FAILED: %s (fail %d, error %s)\n", bad$test[i], bad$failed[i], bad$error[i]))
  }
}
