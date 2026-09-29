# Reviewer re-check, items (a) and (b).
#
# (a) Falsify the whole-block batch's correctness condition. For each
#     block the script reports, from the DESIGN and not from
#     re_row_support(), the largest number of the block's kept columns
#     any row loads. Where that is > 1 the whole-block batch MUST be
#     refused and the split must run. Every batch the shipped code
#     returns is then re-checked against the same independent support,
#     and the shipped Est.Error is compared with the UNBATCHED dense
#     route by identical().
#
# (b) Evaluation COUNTS, not clocks, for the worker's cells A, B and D,
#     on the shipped batcher and on an emulated PRE-NITS batcher (the
#     per-column-position split alone) in the same process. No pre-nits
#     library survives: wt-resmooth-lib was reinstalled at 01:25 and
#     wt-resmooth-lib2 already carries the fix, so the pre-nits arm is
#     reconstructed from its algorithm rather than measured. The 0.64.0
#     arm runs with REVLIB=base.
#
#   Rscript dev/resmooth-rev2-batch.R > dev/resmooth-rev2-batch-lane.txt
#   REVLIB=base Rscript dev/resmooth-rev2-batch.R > dev/resmooth-rev2-batch-base.txt
ARM <- if (identical(Sys.getenv("REVLIB"), "base")) "base" else "lane"
LIB <- if (ARM == "base") "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-resmooth-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("ARM:", ARM, "| frmtmb:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
gt <- function(nm) get(nm, envir = ns)
has_try <- exists("re_batch_try", envir = ns, inherits = FALSE)
cat("re_batch_try present:", has_try, "\n\n")

## The support matrix, rebuilt here from each linear predictor's Z, so
## that it does not inherit re_row_support()'s own reading of the design.
own_support <- function(fit) {
  n_c <- fit$frame[["n_c"]] %||% length(fit$estimates[["b"]])
  out <- NULL
  for (lp in fit$frame[["linpreds"]]) {
    Z <- lp[["Z"]]
    if (is.null(Z)) next
    M <- as.matrix(abs(as.matrix(Z)) > 0) * 1
    out <- if (is.null(out)) M else out + M
  }
  out
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## The PRE-NITS batcher: the per-column-position split alone, which is
## what re_b_batches() did before re_batch_try() was pulled out.
prenits_batches <- function(fit, b_idx, S) {
  if (!length(b_idx)) return(list())
  out <- list()
  for (bk in fit$frame[["re_blocks"]]) {
    ci <- bk[["c_idx"]]; bi <- bk[["b_idx"]]
    if (!any(bi %in% b_idx)) next
    if (!gt("block_b_positionwise")(bk)) return(NULL)
    if (length(ci) != length(bi)) return(NULL)
    D <- max(1L, bk[["dim"]])
    for (k in seq_len(D)) {
      at <- seq.int(k, length(ci), by = D)
      keep <- at[bi[at] %in% b_idx]
      if (!length(keep)) next
      Sb <- S[, ci[keep], drop = FALSE]
      if (max(rowSums(as.matrix(Sb))) > 1) return(NULL)
      owner <- rep(NA_integer_, nrow(Sb))
      nz <- which(as.matrix(Sb) != 0, arr.ind = TRUE)
      owner[nz[, 1L]] <- nz[, 2L]
      out[[length(out) + 1L]] <- list(idx = bi[keep], owner = owner)
    }
  }
  out
}

## The shipped configuration, replicated from fitted_point_se()
ship_cfg <- function(fit, newdata, re_formula, resp = NULL) {
  want <- gt("smooth_b_idx")(fit)
  if (gt("re_form_keeps")(re_formula)) {
    want <- sort(unique(c(want, gt("re_governed_b")(fit))))
  }
  b_idx <- if (!length(want)) NULL else {
    used <- gt("re_used_b")(fit, newdata, resp, FALSE)
    if (is.null(used)) want else intersect(want, used)
  }
  bt <- if (length(b_idx)) {
    gt("re_b_batches")(fit, newdata, resp, FALSE, b_idx)
  }
  list(b_idx = b_idx, bt = bt)
}

counted <- function(fit, f, b_idx, bt) {
  n <- 0L
  g <- function(x) { n <<- n + 1L; f(x) }
  se <- gt("fit_fd_se")(fit, g, b_idx = b_idx, b_batch = bt)
  list(se = as.vector(se), evals = n)
}

report <- function(lab, fit, nd, re_formula) {
  rf <- if (identical(re_formula, "NA")) NA else NULL
  lbl <- if (is.null(rf)) "NULL" else "NA"
  f <- function(x) gt("fitted_point")(x, nd, rf, "response", NULL, NULL,
                                      FALSE)
  cfg <- ship_cfg(fit, nd, rf)
  S <- own_support(fit)
  cat(sprintf("== %s | newdata %s | re_formula = %s\n", lab,
              if (is.null(nd)) "NULL (in sample)" else nrow(nd), lbl))
  bl <- fit$frame[["re_blocks"]]
  for (i in seq_along(bl)) {
    bk <- bl[[i]]
    keep <- which(bk[["b_idx"]] %in% (cfg$b_idx %||% integer(0)))
    if (!length(keep)) {
      cat(sprintf("   block %d %-7s not in b_idx\n", i, bk[["covstruct"]]))
      next
    }
    mx <- max(rowSums(S[, bk[["c_idx"]][keep], drop = FALSE]))
    cat(sprintf(paste0("   block %d %-7s dim %-4d kept coefs %-4d | ",
                       "MY max row load over kept cols: %d -> whole-block ",
                       "batch %s\n"),
                i, bk[["covstruct"]], max(1L, bk[["dim"]]), length(keep),
                mx, if (mx > 1) "MUST be refused" else "is allowed"))
  }
  # every batch the shipped code returned, re-checked against MY support
  nb <- if (is.null(cfg$bt)) NA_integer_ else length(cfg$bt)
  bad <- 0L
  if (!is.null(cfg$bt)) {
    for (b in cfg$bt) {
      cols <- integer(0)
      for (bk in bl) {
        at <- match(b$idx, bk[["b_idx"]])
        if (!anyNA(at)) cols <- bk[["c_idx"]][at]
      }
      if (!length(cols)) { bad <- bad + 1L; next }
      if (max(rowSums(S[, cols, drop = FALSE])) > 1) bad <- bad + 1L
      # and the owner map must agree with my own support
      Sb <- S[, cols, drop = FALSE]
      own <- rep(NA_integer_, nrow(Sb))
      nz <- which(as.matrix(Sb) != 0, arr.ind = TRUE)
      own[nz[, 1L]] <- nz[, 2L]
      if (!identical(own, b$owner)) bad <- bad + 1L
    }
  }
  cat(sprintf("   shipped batches: %s | batches failing MY re-check: %d\n",
              if (is.na(nb)) "NULL (no batching)" else as.character(nb), bad))
  # shipped (batched) against the UNBATCHED dense route
  sh <- counted(fit, f, cfg$b_idx, cfg$bt)
  de <- counted(fit, f, cfg$b_idx, NULL)
  pn <- counted(fit, f, cfg$b_idx,
                if (length(cfg$b_idx)) prenits_batches(fit, cfg$b_idx, S))
  cat(sprintf(paste0("   evals: shipped %-5d pre-nits split %-5d ",
                     "unbatched dense %-5d | b %d\n"),
              sh$evals, pn$evals, de$evals, length(cfg$b_idx %||% 0L)))
  cat(sprintf("   Est.Error shipped vs UNBATCHED: identical %-5s max|d| %.3e\n",
              identical(sh$se, de$se), max(abs(sh$se - de$se))))
  cat(sprintf("   Est.Error shipped vs PRE-NITS : identical %-5s max|d| %.3e\n",
              identical(sh$se, pn$se), max(abs(sh$se - pn$se))))
  # and the shipped value really is what fitted() reports
  ft <- as.vector(fitted(fit, newdata = nd, re_formula = rf)[, "Est.Error", ])
  cat(sprintf("   equals fitted()[, Est.Error]:   identical %-5s\n",
              identical(sh$se, ft)))
  cat("\n")
}

ord <- function(v) {
  factor(cut(v, c(-Inf, -0.8, 0.8, Inf), labels = FALSE), ordered = TRUE)
}

## ---------------- (a) the five falsification blocks ----------------
cat("#### (a) blocks where a row loads TWO of the block's columns ####\n\n")

# 1. (1 + x | g): dim 2, a row loads its level's intercept AND slope
set.seed(71)
ngJ <- 10; perJ <- 24
dJ <- data.frame(g = factor(rep(seq_len(ngJ), each = perJ)),
                 x = stats::runif(ngJ * perJ, -2, 2))
dJ$y <- ord(sin(1.5 * dJ$x) + stats::rnorm(ngJ, 0, 0.7)[dJ$g] +
              stats::rnorm(ngJ, 0, 0.5)[dJ$g] * dJ$x + stats::rlogis(nrow(dJ)))
fJ <- suppressWarnings(frm(bf(y ~ s(x, k = 6) + (1 + x | g)),
                          family = cumulative(), data = dJ))
report("J cumulative s(x,k=6) + (1 + x | g)", fJ, NULL, "NULL")
report("J cumulative s(x,k=6) + (1 + x | g)", fJ,
       dJ[c(2, 50, 200), c("x", "g")], "NULL")

# 2. multi-membership: a row loads TWO LEVELS at the SAME position
set.seed(73)
nM <- 220
dM <- data.frame(g1 = factor(sample(letters[1:7], nM, TRUE)),
                 g2 = factor(sample(letters[1:7], nM, TRUE)),
                 x = stats::runif(nM, -2, 2))
uM <- stats::rnorm(7, 0, 0.7)
dM$y <- ord(sin(1.5 * dM$x) + uM[dM$g1] + uM[dM$g2] + stats::rlogis(nM))
fM <- suppressWarnings(frm(bf(y ~ s(x, k = 6) + (1 | mm(g1, g2))),
                          family = cumulative(), data = dM))
report("M cumulative s(x,k=6) + (1 | mm(g1,g2))", fM, NULL, "NULL")

# 3. an fs factor smooth
set.seed(21)
ngF <- 6
dF <- data.frame(g = factor(rep(seq_len(ngF), each = 25)),
                 x = stats::runif(ngF * 25, -2, 2))
dF$y <- ord(stats::rnorm(ngF, 0, 1)[dF$g] * sin(dF$x) + 0.5 * dF$x +
              stats::rlogis(nrow(dF)))
fF <- suppressWarnings(frm(bf(y ~ s(x, g, bs = "fs", k = 5)),
                          family = cumulative(), data = dF))
report("FS cumulative s(x, g, bs = 'fs', k = 5)", fF, NULL, "NA")
report("FS cumulative s(x, g, bs = 'fs', k = 5)", fF,
       dF[c(3, 40, 90), c("x", "g")], "NA")

# 4. s(x, by = g): a by-factor smooth, one block per level
set.seed(75)
ngB <- 5
dB <- data.frame(g = factor(rep(seq_len(ngB), each = 44)),
                 x = stats::runif(ngB * 44, -2, 2))
dB$y <- ord(sin(1.5 * dB$x) * stats::rnorm(ngB, 1, 0.5)[dB$g] +
              stats::rlogis(nrow(dB)))
fB <- suppressWarnings(frm(bf(y ~ g + s(x, by = g, k = 5)),
                          family = cumulative(), data = dB))
report("BY cumulative g + s(x, by = g, k = 5)", fB, NULL, "NA")

# 5. a t2() with an re margin
set.seed(77)
ngT <- 8
dT <- data.frame(g = factor(rep(seq_len(ngT), each = 28)),
                 x = stats::runif(ngT * 28, -2, 2))
dT$y <- ord(sin(1.5 * dT$x) + stats::rnorm(ngT, 0, 0.8)[dT$g] * dT$x +
              stats::rlogis(nrow(dT)))
fT <- suppressWarnings(frm(bf(y ~ t2(x, g, bs = c("cr", "re"), k = 5)),
                          family = cumulative(), data = dT))
report("T2 cumulative t2(x, g, c('cr','re'), k = 5)", fT, NULL, "NA")
report("T2 cumulative t2(x, g, c('cr','re'), k = 5)", fT,
       dT[c(4, 60, 150), c("x", "g")], "NA")

## ---------------- (b) the worker's cells A, B, D ----------------
cat("#### (b) the worker's cells, by evaluation COUNT ####\n\n")
set.seed(41)
nD <- 160
dG <- data.frame(x = stats::runif(nD, -2, 2))
dG$y <- ord(1.1 * sin(2 * dG$x) + stats::rlogis(nD))
fG <- suppressWarnings(frm(bf(y ~ gp(x)), family = cumulative(), data = dG))
cat("gp(x) n_b =", length(fG$estimates[["b"]]), "\n")
# A: newdata at three positions the fit never saw
report("A gp(x), newdata 3 NEW positions", fG,
       data.frame(x = c(-1.234, 0.777, 1.611)), "NA")
# B: in sample
report("B gp(x), in sample", fG, NULL, "NA")
# D: s(x, k = 8) + (1 | g), 40 levels, in sample, NULL
set.seed(31)
ngD <- 40; perD <- 20
dD <- data.frame(g = factor(rep(seq_len(ngD), each = perD)),
                 x = stats::runif(ngD * perD, -2, 2))
dD$y <- ord(sin(1.5 * dD$x) + stats::rnorm(ngD, 0, 0.8)[dD$g] +
              stats::rlogis(nrow(dD)))
fD <- suppressWarnings(frm(bf(y ~ s(x, k = 8) + (1 | g)),
                          family = cumulative(), data = dD))
report("D s(x,k=8)+(1|g) 40 levels, in sample", fD, NULL, "NULL")
report("E s(x,k=8)+(1|g) 40 levels, in sample", fD, NULL, "NA")
cat("DONE\n")
