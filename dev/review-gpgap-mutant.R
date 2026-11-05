# Review of lane gpgap, item 1: does a wrong model still FAIL through the
# instrument's new branch? Each mutant lives in this process only
# (namespace or instrument environment); nothing is installed.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-mutant.R <log> <id> <mutant>
args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
mutant <- args[2]
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
requireNamespace("brms"); requireNamespace("rstan")
attach(list(me = brms::me), name = "equiv_brms_me", warn.conflicts = FALSE)
eq_env <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = eq_env)
}
ns <- asNamespace("frmtmb")
put <- function(name, val) {
  unlockBinding(name, ns)
  assign(name, val, envir = ns)
  lockBinding(name, ns)
}
kern_mut <- function(cc) {
  put("gp_corr", function(theta, blk) {
    Q <- 0
    for (j in seq_along(blk[["aux_D2"]])) {
      rho <- if (isTRUE(blk[["gp_iso"]])) exp(theta[2]) else exp(theta[1 + j])
      Q <- Q + blk[["aux_D2"]][[j]] / (2 * rho^2)
    }
    exp(-cc * Q) + diag(1e-6, nrow(blk[["aux_D2"]][[1]]))
  })
}
tr_mut <- function(fn) {
  orig <- eq_env$stan_pars_from_fit
  eq_env$stan_pars_from_fit <- function(...) {
    p <- orig(...)
    fn(p)
  }
}
switch(mutant,
  none = NULL,
  `kern_1.1` = kern_mut(1.1),
  `kern_1.01` = kern_mut(1.01),
  `kern_1.001` = kern_mut(1.001),
  kern_rq = put("gp_corr", function(theta, blk) {
    Q <- blk[["aux_D2"]][[1]] / (2 * exp(2 * theta[2]))
    (1 + Q / 2)^(-2) + diag(1e-6, nrow(Q))
  }),
  `nug_1e-5` = put("gp_nugget", 1e-5),
  `nug_2e-6` = put("gp_nugget", 2e-6),
  `nug_1.1e-6` = put("gp_nugget", 1.1e-6),
  `nug_1e-7` = put("gp_nugget", 1e-7),
  `tr_lscale_1.01` = tr_mut(function(p) {
    n <- grep("^lscale_", names(p), value = TRUE)
    for (q in n) p[[q]] <- p[[q]] * 1.01
    p
  }),
  `tr_lscale_1.0001` = tr_mut(function(p) {
    n <- grep("^lscale_", names(p), value = TRUE)
    for (q in n) p[[q]] <- p[[q]] * 1.0001
    p
  }),
  `tr_sdgp_1.01` = tr_mut(function(p) {
    n <- grep("^sdgp_", names(p), value = TRUE)
    for (q in n) p[[q]] <- p[[q]] * 1.01
    p
  }),
  `tr_sdgp_1.0001` = tr_mut(function(p) {
    n <- grep("^sdgp_", names(p), value = TRUE)
    for (q in n) p[[q]] <- p[[q]] * 1.0001
    p
  }),
  tr_noorder = {
    eq_env$brms_gp_order <- function(sdat, zname) {
      seq_len(NROW(sdat[[sub("^zgp_", "Xgp_", zname)]]))
    }
  },
  stop("unknown mutant"))

smp <- read.delim("dev/equiv-data/sample.tsv")
row <- smp[smp$id == id, ]
key <- paste(row$feature_a, row$feature_b, sep = " x ")
rec <- eq_env$eq_measure(key, row)
cat(sprintf(paste0("MUT %s %s: verdict %s, level %s, l1_ok %s (val %s grad %s ",
                   "pts %s), const %.6g, extra %s, grad %.3g, tol %.3g, ",
                   "own line: %s\n"),
            id, mutant, rec$verdict, paste(rec$level, collapse = "+"),
            rec$l1_ok, rec$l1_val_ok, rec$l1_grad_ok, rec$l1_pts_ok,
            rec$l1_const %||% NA, paste(signif(rec$l1_extra_const, 4),
                                        collapse = ","),
            rec$l1_grad %||% NA, 1e-6 * max(1, abs(rec$l1_ours %||% NA)),
            rec$l1_by_decision %||% "(none)"))
