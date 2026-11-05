# Short check, item 7: mutants of the fix round, in-process, on the lane's test
# file. Usage: Rscript review-latres-mut2.R <mutant> <testfile>
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
args <- commandArgs(TRUE)
mut <- args[[1]]
tf <- args[[2]]
Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true")
suppressPackageStartupMessages(library(testthat))
ns <- asNamespace("frmtmb")
setf <- function(nm, f) {
  environment(f) <- ns
  assignInNamespace(nm, f, ns = "frmtmb")
}
switch(mut,
  none = NULL,
  m1_no_warning = setf("latres_laplace_report", function(fit) invisible(NULL)),
  r1_no_box = setf("latres_disp_upper", function(frame, spec, upper, map) upper),
  r1_no_cap = setf("latres_disp_guard", function(v) v),
  m4_no_margin = setf("latres_refuse_margin", function(fit) invisible(NULL)),
  r3_floor0 = assignInNamespace("latres_cubic_floor", 0, ns = "frmtmb"),
  r2_cosy_any = {
    f <- get("latres_limit_report", ns)
    txt <- paste(deparse(f, width.cutoff = 500L), collapse = "
")
    new <- sub("chain && val[2L] > 0", "chain", txt, fixed = TRUE)
    stopifnot(!identical(new, txt))
    setf("latres_limit_report", eval(parse(text = new), ns))
  },
  r4_skip_new = setf("latres_skip_se", function(fit, newdata = NULL) nrow(fit$frame$data_frame %||% data.frame()) > 3000 || fit$frame$n_obs > getOption("frmtmb.latres_se_rows", 3000L)),
  m2_no_hold = setf("latres_held_pos", function(fit, map) integer(0)),
  m3_no_report = setf("latres_limit_report", function(fit, control) invisible(NULL)),
  m3_no_na = setf("latres_at_limit", function(theta, a) rep(FALSE, 1L + length(latres_natural_coefs(theta[-1L], a)))),
  m4_no_refusal = setf("latres_refuse_laplace", function(spec, frame) invisible(NULL)),
  m5_old_guard = {
    f <- get("check_autocor_latent", ns)
    txt <- paste(deparse(f, width.cutoff = 500L), collapse = "\n")
    new <- gsub('"mu" %in% resp$nl_dpars', 'length(resp$nlpars) || length(resp$nl_dpars)', txt, fixed = TRUE)
    stopifnot(!identical(new, txt))
    setf("check_autocor_latent", eval(parse(text = new), ns))
  },
  m6_at_zero = setf("latres_fitted_mean", function(fit, newdata, re_formula, resp, dpar, allow_new_levels, out) out),
  m6_half_var = {
    f <- get("response_mean_gh", ns)
    txt <- paste(deparse(f, width.cutoff = 500L), collapse = "\n")
    new <- sub("sqrt(rep_len(var, length(eta)))", "sqrt(rep_len(var / 2, length(eta)))", txt, fixed = TRUE)
    stopifnot(!identical(new, txt))
    setf("response_mean_gh", eval(parse(text = new), ns))
  },
  m7_no_skip = setf("latres_skip_se", function(fit) FALSE),
  stop("unknown"))
r <- test_file(tf, package = "frmtmb", env = test_env("frmtmb"), reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(r)
cat(sprintf("RESULT %-16s %-22s tests=%d fail=%d err=%d skip=%d\n", mut, basename(tf), sum(df$nb), sum(df$failed), sum(df$error), sum(df$skipped)))
bad <- df[df$failed > 0 | df$error, "test"]
if (length(bad)) cat("  failed:", paste(substr(unique(bad), 1, 80), collapse = " | "), "\n")
