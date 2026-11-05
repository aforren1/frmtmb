# Attack 10: mutants, applied in-process with assignInNamespace (no build, no
# library touched). Each mutant must fail the lane's tests.
# Usage: Rscript review-latres-mutants.R <mutant> <testfile> [pkg]
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
args <- commandArgs(TRUE)
mut <- args[[1]]
tf <- args[[2]]
pkg <- if (length(args) > 2) args[[3]] else "frmtmb"
Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true")
suppressPackageStartupMessages({library(testthat); library(pkg, character.only = TRUE)})
ns <- asNamespace("frmtmb")
setf <- function(nm, f) {
  environment(f) <- ns
  assignInNamespace(nm, f, ns = "frmtmb")
}
patch_body <- function(nm, from, to) {
  f <- get(nm, ns)
  txt <- paste(deparse(f, width.cutoff = 500L), collapse = "\n")
  stopifnot(grepl(from, txt, fixed = TRUE))
  txt <- sub(from, to, txt, fixed = TRUE)
  g <- eval(parse(text = txt), ns)
  setf(nm, g)
}
switch(mut,
  none = NULL,
  ma_flip = patch_body("latres_expand", "ma1 = s * (b[seq_len(n)] + cf$ma * b[a[[\"prev\"]]])", "ma1 = s * (b[seq_len(n)] - cf$ma * b[a[[\"prev\"]]])"),
  arma_flip = patch_body("latres_expand", "arma1 = s * (b[seq_len(n)] + cf$ma * b[a[[\"prev\"]]])", "arma1 = s * (b[seq_len(n)] - cf$ma * b[a[[\"prev\"]]])"),
  start_var = patch_body("latres_lpdf", "0.5 * sum(s0 * s0)/(1 + t1 * t1)", "0.5 * sum(s0 * s0)"),
  cosy_swap = patch_body("latres_expand", "s * (sqrt(1 - cf$rho) * b[seq_len(n)] + sqrt(cf$rho) * b[a[[\"u\"]]])", "s * (sqrt(cf$rho) * b[seq_len(n)] + sqrt(1 - cf$rho) * b[a[[\"u\"]]])"),
  condar_sign = patch_body("latres_expand", "s * (b + latres_lag_sum(b, cf$ar, a))", "s * (b - latres_lag_sum(b, cf$ar, a))"),
  no_on_response = setf("latres_on_response", function(fam, av) FALSE),
  all_on_response = setf("latres_on_response", function(fam, av) TRUE),
  cdf_ignored = patch_body("latres_on_response", "(!cdf || identical(fam[[\"family\"]], \"cox\"))", "TRUE"),
  no_guard = setf("latres_guard", function(mu, a) mu),
  methods_no_shift = setf("latres_shift_of", function(fit, bk) 0 * fit$estimates[["b"]][bk[["b_idx"]]][seq_len(bk$aux$n)]),
  new_rows_no_draw = setf("latres_new_draw", function(fs, lat) 0 * seq_len(lat$blk$aux$n)),
  refuse_off = setf("latres_refuse_newdata", function(fit, newdata) invisible(NULL)),
  sderr_scale2 = patch_body("latres_expand", "ar1 = s * b", "ar1 = 2 * s * b"),
  row_order = {
    # the chain follows the data's row order, not the time order
    g <- get("autocor_cond_positions", ns)
    setf("autocor_cond_positions", function(gidx, tidx, n, p) g(gidx, seq_len(n), n, p))
  },
  sparse_dense = patch_body("latres_aux", "aux$prev <- prev", "prev[is.na(prev)] <- n + 1L; aux$prev <- prev"),
  stop("unknown mutant"))
# the registry holds the closures it was built with
reg <- get("frmtmb_covstruct_registry", ns)
reg$reg$latres$expand$apply <- get("latres_expand", ns)
reg$reg$latres$expand$jacobian <- get("latres_jacobian", ns)
reg$reg$latres$nll <- get("latres_lpdf", ns)
reg$reg$latres$draw <- get("latres_draw", ns)
r <- test_file(tf, package = pkg, env = test_env(pkg), reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(r)
cat(sprintf("RESULT %-18s %-28s tests=%d pass=%d fail=%d err=%d skip=%d\n", mut, basename(tf), sum(df$nb), sum(df$nb) - sum(df$failed) - sum(df$error) - sum(df$skipped), sum(df$failed), sum(df$error), sum(df$skipped)))
bad <- df[df$failed > 0 | df$error, "test"]
if (length(bad)) cat("  failed:", paste(substr(unique(bad), 1, 70), collapse = " | "), "\n")
