# Reviewer of lane ordstarts, build round. Probe 3 (question 3): what the
# holds check (mixture_swap_holds()) rejects. A wrapper runs the same four
# jittered runs (same private stream, same sd), returns the same verdict,
# and records each run's end and its class (degen_at() / tier_of() of the
# calling mixture_starts() frame): 0 proper, 1 flat, 2 step function.
#   Rscript ordstarts-build-review-p3.R <re> <n> <G> <seeds> <out.tsv>
args <- commandArgs(trailingOnly = TRUE)
re <- args[1]
nobs <- as.integer(args[2])
G <- as.integer(args[3])
seeds <- eval(parse(text = args[4]))
out <- args[5]
.libPaths(c("C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
ns <- asNamespace("frmtmb")
gen <- function(seed) {
  set.seed(seed)
  x <- rnorm(nobs)
  z <- rnorm(nobs)
  cls <- rbinom(nobs, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(nobs)
  d <- data.frame(x = x, z = z)
  if (re != "none") {
    g <- rep_len(seq_len(G), nobs)
    u0 <- rnorm(G, 0, 0.5)
    lat <- lat + u0[g]
    d$g <- factor(g)
  }
  d$y <- 1L + (lat > -1 + 0.1 * z) + (lat > 0) + (lat > 1 - 0.1 * z)
  d
}
log <- new.env()
log$rows <- list()
wrap <- function(o, run_from, np, tol) {
  e <- environment(run_from)
  Z <- ns$mixture_starts_draws(np, ns$mixture_swap_jitter_runs,
                               seed = 20261117L)
  ends <- lapply(seq_len(ncol(Z)), function(j) {
    run_from(o$par + ns$mixture_swap_jitter_sd * Z[, j])
  })
  climb <- vapply(ends, function(r) {
    !is.null(r) && r$objective < o$objective - tol
  }, NA)
  tiers <- vapply(ends, function(r) {
    if (is.null(r)) return(NA_integer_)
    dg <- e$degen_at(r)
    if (is.null(dg)) NA_integer_ else as.integer(e$tier_of(dg))
  }, 0L)
  log$rows[[length(log$rows) + 1L]] <- list(
    A = e$fit$opt$objective, S = o$objective,
    ends = vapply(ends, function(r) if (is.null(r)) NA_real_ else
      r$objective, 0), tiers = tiers, holds = !any(climb))
  !any(climb)
}
utils::assignInNamespace("mixture_swap_holds", wrap, "frmtmb")
form <- if (re == "none") bf(y ~ x + z) else bf(y ~ x + z + (1 | g))
first <- !file.exists(out)
for (seed in seeds) {
  log$rows <- list()
  f <- tryCatch(suppressWarnings(suppressMessages(
    frm(form, family = mixture(cumulative(), sratio()), data = gen(seed)))),
    error = function(e) e)
  if (inherits(f, "error")) next
  ms <- f$opt$mixture_starts
  for (r in log$rows) {
    row <- data.frame(re = re, n = nobs, seed = seed,
                      gainS = r$A - r$S, holds = r$holds,
                      swapped = ms[["swapped"]],
                      llfit = as.numeric(logLik(f)),
                      ends_minus_S = paste(sprintf("%+.4g", r$S - r$ends),
                                           collapse = ","),
                      tiers = paste(r$tiers, collapse = ","))
    utils::write.table(row, out, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = first, append = !first)
    first <- FALSE
  }
  cat(re, nobs, "seed", seed, "holds calls", length(log$rows), "\n")
}
cat("DONE p3\n")
