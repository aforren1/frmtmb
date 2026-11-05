# Does test-gp-noncentered.R kill each mutant? One mutant per process.
#   Rscript dev/ncgp-mut.R <mutant>
m <- commandArgs(TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
ns <- asNamespace("frmtmb")
put <- function(nm, f) {
  unlockBinding(nm, ns)
  environment(f) <- ns
  assign(nm, f, envir = ns)
}
switch(m,
  control = NULL,
  # brms's absolute line: 1e-12 on the covariance whatever sd is
  absolute = put("gp_corr", function(theta, blk) {
    Q <- 0
    for (j in seq_along(blk[["aux_D2"]])) {
      rho <- if (isTRUE(blk[["gp_iso"]])) exp(theta[2]) else exp(theta[1 + j])
      Q <- Q + blk[["aux_D2"]][[j]] / (2 * rho^2)
    }
    exp(-Q) + diag(1, nrow(blk[["aux_D2"]][[1]])) * 1e-12 * exp(-2 * theta[1])
  }),
  nugget6 = put("gp_nugget", 1e-6),
  nonull = put("gp_null_at_field", function(fit, null, r, f) null),
  noscan = put("gp_start_scan", function(obj, template, frame, start, to_beta = identity) obj),
  nowarm = put("gp_marginal_max", function(frame, keep = 3L, distinct = 4L) NULL),
  nounion = put("gp_marginal_opt", local({
    f <- get("gp_marginal_opt", asNamespace("frmtmb"))
    function(s, keep = 3L, distinct = 4L) f(s, 3L, 0L)
  })),
  nopolish = put("gp_marginal_polish", function(s, out) out),
  nowork = put("gp_start_working", function(frame) NULL),
  quiet_all = put("gp_nugget_message", function(fit) NULL),
  # F3 by headroom alone, the clause the review rejected
  quiet_headroom = {
    s <- paste(deparse(ns$check_convergence), collapse = "\n")
    s2 <- sub("isTRUE\\(abs\\(fit\\$opt\\$objective[^<]*<=\\s*control\\$grad_tol\\)", "TRUE", s)
    stopifnot(!identical(s, s2))
    put("check_convergence", eval(parse(text = s2)))
  },
  floorless = put("se_empty_row", 0),
  norowcheck = put("se_row_curved", function(j, fit, p, u) FALSE),
  noorder = put("gp_tape_order", function(frame, tpl) tpl),
  sorted = put("gp_farthest_first", function(blk) seq_len(nrow(blk[["aux_D2"]][[1L]]))),
  nothcols = put("sm_delta_cols", function(object, A, coef_pos, Zc, bk, b_pos, th_pos, rrj, krig = NULL) {
    if (block_b_positionwise(bk)) {
      return(list(A = Matrix::cbind2(A, Zc), coef_pos = c(coef_pos, b_pos[bk[["b_idx"]]])))
    }
    rrj <- rrj %||% rr_jacobians(object)
    list(A = Matrix::cbind2(A, Zc %*% rrj$Jb[bk[["c_idx"]], bk[["b_idx"]], drop = FALSE]),
         coef_pos = c(coef_pos, b_pos[bk[["b_idx"]]]))
  }),
  fd_dtheta = put("gp_field_dtheta", function(z, theta, blk) {
    L <- gp_factor_ordered(theta, blk)
    out <- matrix(0, length(z), length(theta))
    for (k in seq_along(theta)) {
      e <- replace(numeric(length(theta)), k, 1e-5)
      df <- (gp_field(z, theta + e, blk) - gp_field(z, theta - e, blk)) / 2e-5
      out[, k] <- backsolve(L, df[blk[["gp_order"]]], upper.tri = FALSE) / sd_floored(theta[1])
    }
    out
  }),
  stop("unknown mutant"))
r <- test_dir("tests/testthat", filter = "gp-noncentered|gp-start", package = "frmtmb",
               env = test_env("frmtmb"), reporter = "silent", stop_on_failure = FALSE)
res <- as.data.frame(r)
cat(sprintf("NCGPM %s pass=%d fail=%d error=%d: %s\n", m, sum(res$passed), sum(res$failed),
            sum(as.integer(res$error)),
            paste(unique(res$test[res$failed > 0 | res$error]), collapse = " | ")))
