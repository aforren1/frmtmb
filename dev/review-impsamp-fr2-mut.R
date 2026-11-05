# Fix round 2, item 7: mutants on the final source, sites read from the final code.
# make <m>: source copy. test <m> [desc]: test-importance.R, or one named test.
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
  capped_settles = list(list("settled = !rose && !capped", "settled = !rose")),
  no_refusal = list(list(
    paste0("  if (it$rose) {\n    imp_refuse_draws(it, n_draw, iterate, lay, frame[[\"n_obs\"]],\n",
           "                     control$importance_ess %||% imp_ess_floor)\n  }\n"), "")),
  no_units = list(list("    u <- if (key %in% names(units)) units[[key]] else 1\n", "    u <- 1\n")),
  # item 3: the old size on theta alone, everything else as built
  theta_size = list(list(
    "  h <- rep(imp_cross_step, length(at))\n",
    paste0("  h <- rep(imp_cross_step, length(at))\n",
           "  th_ <- names(at) == \"theta\"\n",
           "  h[th_] <- imp_cross_step * pmax(1, abs(unname(at)[th_]))\n"))),
  no_prior_scale = list(list(
    "  prior_t <- imp_prior_scales(lay, theta_a)\n",
    "  prior_t <- lapply(lay[[\"q\"]], function(q) rep(1, q))\n"))
)
if (cmd == "make") {
  dst <- file.path(SRC, paste0("fr2-", mut))
  unlink(dst, recursive = TRUE)
  dir.create(file.path(dst, "R"), recursive = TRUE)
  dir.create(file.path(LIB, mut), recursive = TRUE, showWarnings = FALSE)
  for (f in c("DESCRIPTION", "NAMESPACE")) file.copy(file.path(WT, f), dst)
  file.copy(list.files(file.path(WT, "R"), full.names = TRUE), file.path(dst, "R"))
  for (d in c("inst", "data", "man")) file.copy(file.path(WT, d), dst, recursive = TRUE)
  p <- file.path(dst, "R", "importance.R")
  s <- paste0(paste(readLines(p), collapse = "\n"), "\n")
  for (ed in muts[[mut]]) {
    if (!grepl(ed[[1]], s, fixed = TRUE) || lengths(gregexpr(ed[[1]], s, fixed = TRUE)) != 1L) {
      stop("site not once: ", ed[[1]])
    }
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
