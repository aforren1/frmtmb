# Reviewer re-check, items (a) and (b), second pass.
#
# The first pass (dev/resmooth-rev2-batch-lane.txt) reported my
# independent re-check FAILING on every newdata cell. That was MY
# instrument, not the code: I built the support from each linear
# predictor's in-sample `Z`, which has n_obs rows, while a batch on the
# newdata route carries one owner per NEWDATA row. The comparison was
# between vectors of different length. The same flaw sent my pre-nits
# emulation down the dense branch on those cells, because fit_fd_se()
# skips a batch whose owner length does not match nrow(m0).
#
# So this pass restricts the two checks that need a design I built
# myself to the IN-SAMPLE route, where `Z` IS the design, and says so.
# On the newdata route the design can only be rebuilt by the package's
# own pred_design(), so an independent support is not available there;
# what IS available there, and is used, is the behavioural comparison
# against the unbatched dense route, which shares no batching code.
#
#   Rscript dev/resmooth-rev2-batch2.R > dev/resmooth-rev2-batch2-lane.txt
#   REVLIB=base Rscript dev/resmooth-rev2-batch2.R > dev/resmooth-rev2-batch2-base.txt
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
`%||%` <- function(a, b) if (is.null(a)) b else a
cat("re_batch_try present:",
    exists("re_batch_try", envir = ns, inherits = FALSE), "\n\n")

own_support <- function(fit) {
  out <- NULL
  for (lp in fit$frame[["linpreds"]]) {
    Z <- lp[["Z"]]
    if (is.null(Z)) next
    M <- (abs(as.matrix(Z)) > 0) * 1
    out <- if (is.null(out)) M else out + M
  }
  out
}

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
      if (max(rowSums(Sb)) > 1) return(NULL)
      owner <- rep(NA_integer_, nrow(Sb))
      nz <- which(Sb != 0, arr.ind = TRUE)
      owner[nz[, 1L]] <- nz[, 2L]
      out[[length(out) + 1L]] <- list(idx = bi[keep], owner = owner)
    }
  }
  out
}

ship_cfg <- function(fit, newdata, rf, resp = NULL) {
  want <- gt("smooth_b_idx")(fit)
  if (gt("re_form_keeps")(rf)) {
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
  list(se = as.vector(gt("fit_fd_se")(fit, g, b_idx = b_idx, b_batch = bt)),
       evals = n)
}
ulp <- function(a, b) {
  d <- abs(a - b) / (.Machine$double.eps * abs(a))
  max(d[is.finite(d)])
}

report <- function(lab, fit, nd, rf) {
  lbl <- if (is.null(rf)) "NULL" else "NA"
  insample <- is.null(nd)
  f <- function(x) gt("fitted_point")(x, nd, rf, "response", NULL, NULL,
                                      FALSE)
  cfg <- ship_cfg(fit, nd, rf)
  cat(sprintf("== %s | %s | re_formula = %s\n", lab,
              if (insample) "IN SAMPLE" else paste("newdata", nrow(nd)), lbl))
  nb <- if (is.null(cfg$bt)) NA_integer_ else length(cfg$bt)
  if (insample) {
    S <- own_support(fit)
    for (i in seq_along(fit$frame[["re_blocks"]])) {
      bk <- fit$frame[["re_blocks"]][[i]]
      keep <- which(bk[["b_idx"]] %in% (cfg$b_idx %||% integer(0)))
      if (!length(keep)) {
        cat(sprintf("   block %d %-7s not differenced\n", i,
                    bk[["covstruct"]])); next
      }
      mx <- max(rowSums(S[, bk[["c_idx"]][keep], drop = FALSE]))
      cat(sprintf(paste0("   block %d %-7s dim %-4d coefs %-4d | MY max row",
                         " load %d -> whole-block batch %s\n"),
                  i, bk[["covstruct"]], max(1L, bk[["dim"]]), length(keep),
                  mx, if (mx > 1) "MUST be refused" else "allowed"))
    }
    bad <- 0L
    for (b in cfg$bt %||% list()) {
      cols <- integer(0)
      for (bk in fit$frame[["re_blocks"]]) {
        at <- match(b$idx, bk[["b_idx"]])
        if (!anyNA(at)) cols <- bk[["c_idx"]][at]
      }
      Sb <- S[, cols, drop = FALSE]
      own <- rep(NA_integer_, nrow(Sb))
      nz <- which(Sb != 0, arr.ind = TRUE)
      own[nz[, 1L]] <- nz[, 2L]
      if (!length(cols) || max(rowSums(Sb)) > 1 ||
            !identical(own, b$owner)) bad <- bad + 1L
    }
    cat(sprintf(paste0("   shipped batches %s | failing MY independent",
                       " re-check: %d\n"),
                if (is.na(nb)) "NULL (none)" else as.character(nb), bad))
  } else {
    cat(sprintf(paste0("   shipped batches %s | no independent support on",
                       " the newdata route\n"),
                if (is.na(nb)) "NULL (none)" else as.character(nb)))
  }
  sh <- counted(fit, f, cfg$b_idx, cfg$bt)
  sh2 <- counted(fit, f, cfg$b_idx, cfg$bt)     # determinism control
  de <- counted(fit, f, cfg$b_idx, NULL)
  cat(sprintf("   evals shipped %-5d unbatched dense %-5d | b %d\n",
              sh$evals, de$evals, length(cfg$b_idx %||% integer(0))))
  cat(sprintf("   CONTROL same call twice: identical %s\n",
              identical(sh$se, sh2$se)))
  cat(sprintf(paste0("   shipped vs UNBATCHED: identical %-5s max|d| %.3e",
                     " | %.1f ulp\n"),
              identical(sh$se, de$se), max(abs(sh$se - de$se)),
              ulp(sh$se, de$se)))
  if (insample) {
    pn <- counted(fit, f, cfg$b_idx,
                  if (length(cfg$b_idx)) prenits_batches(fit, cfg$b_idx,
                                                         own_support(fit)))
    cat(sprintf(paste0("   evals PRE-NITS split %-5d | shipped vs PRE-NITS:",
                       " identical %-5s max|d| %.3e\n"),
                pn$evals, identical(sh$se, pn$se), max(abs(sh$se - pn$se))))
  }
  ft <- as.vector(fitted(fit, newdata = nd, re_formula = rf)[, "Est.Error", ])
  cat(sprintf("   equals fitted()[, Est.Error]: identical %s\n\n",
              identical(sh$se, ft)))
}

ord <- function(v) {
  factor(cut(v, c(-Inf, -0.8, 0.8, Inf), labels = FALSE), ordered = TRUE)
}

cat("#### (a) blocks where a row loads TWO of the block's columns ####\n\n")
set.seed(71)
ngJ <- 10; perJ <- 24
dJ <- data.frame(g = factor(rep(seq_len(ngJ), each = perJ)),
                 x = stats::runif(ngJ * perJ, -2, 2))
dJ$y <- ord(sin(1.5 * dJ$x) + stats::rnorm(ngJ, 0, 0.7)[dJ$g] +
              stats::rnorm(ngJ, 0, 0.5)[dJ$g] * dJ$x + stats::rlogis(nrow(dJ)))
fJ <- suppressWarnings(frm(bf(y ~ s(x, k = 6) + (1 + x | g)),
                          family = cumulative(), data = dJ))
report("J s(x,k=6) + (1 + x | g)", fJ, NULL, NULL)
report("J s(x,k=6) + (1 + x | g)", fJ, dJ[c(2, 50, 200), c("x", "g")], NULL)

set.seed(73)
nM <- 220
dM <- data.frame(g1 = factor(sample(letters[1:7], nM, TRUE)),
                 g2 = factor(sample(letters[1:7], nM, TRUE)),
                 x = stats::runif(nM, -2, 2))
uM <- stats::rnorm(7, 0, 0.7)
dM$y <- ord(sin(1.5 * dM$x) + uM[dM$g1] + uM[dM$g2] + stats::rlogis(nM))
fM <- suppressWarnings(frm(bf(y ~ s(x, k = 6) + (1 | mm(g1, g2))),
                          family = cumulative(), data = dM))
report("M s(x,k=6) + (1 | mm(g1,g2))", fM, NULL, NULL)

set.seed(21)
ngF <- 6
dF <- data.frame(g = factor(rep(seq_len(ngF), each = 25)),
                 x = stats::runif(ngF * 25, -2, 2))
dF$y <- ord(stats::rnorm(ngF, 0, 1)[dF$g] * sin(dF$x) + 0.5 * dF$x +
              stats::rlogis(nrow(dF)))
fF <- suppressWarnings(frm(bf(y ~ s(x, g, bs = "fs", k = 5)),
                          family = cumulative(), data = dF))
report("FS s(x, g, bs = 'fs', k = 5)", fF, NULL, NA)
report("FS s(x, g, bs = 'fs', k = 5)", fF, dF[c(3, 40, 90), c("x", "g")], NA)

set.seed(75)
ngB <- 5
dB <- data.frame(g = factor(rep(seq_len(ngB), each = 44)),
                 x = stats::runif(ngB * 44, -2, 2))
dB$y <- ord(sin(1.5 * dB$x) * stats::rnorm(ngB, 1, 0.5)[dB$g] +
              stats::rlogis(nrow(dB)))
fB <- suppressWarnings(frm(bf(y ~ g + s(x, by = g, k = 5)),
                          family = cumulative(), data = dB))
report("BY g + s(x, by = g, k = 5)", fB, NULL, NA)

set.seed(77)
ngT <- 8
dT <- data.frame(g = factor(rep(seq_len(ngT), each = 28)),
                 x = stats::runif(ngT * 28, -2, 2))
dT$y <- ord(sin(1.5 * dT$x) + stats::rnorm(ngT, 0, 0.8)[dT$g] * dT$x +
              stats::rlogis(nrow(dT)))
fT <- suppressWarnings(frm(bf(y ~ t2(x, g, bs = c("cr", "re"), k = 5)),
                          family = cumulative(), data = dT))
report("T2 t2(x, g, c('cr','re'), k = 5)", fT, NULL, NA)
report("T2 t2(x, g, c('cr','re'), k = 5)", fT, dT[c(4, 60, 150), c("x", "g")],
       NA)

cat("#### (b) the worker's cells A, B, D, E by evaluation COUNT ####\n\n")
set.seed(41)
nD <- 160
dG <- data.frame(x = stats::runif(nD, -2, 2))
dG$y <- ord(1.1 * sin(2 * dG$x) + stats::rlogis(nD))
fG <- suppressWarnings(frm(bf(y ~ gp(x)), family = cumulative(), data = dG))
cat("gp(x) n_b =", length(fG$estimates[["b"]]), "\n")
report("A gp(x), 3 NEW positions", fG,
       data.frame(x = c(-1.234, 0.777, 1.611)), NA)
report("B gp(x)", fG, NULL, NA)
set.seed(31)
ngD <- 40; perD <- 20
dD <- data.frame(g = factor(rep(seq_len(ngD), each = perD)),
                 x = stats::runif(ngD * perD, -2, 2))
dD$y <- ord(sin(1.5 * dD$x) + stats::rnorm(ngD, 0, 0.8)[dD$g] +
              stats::rlogis(nrow(dD)))
fD <- suppressWarnings(frm(bf(y ~ s(x, k = 8) + (1 | g)),
                          family = cumulative(), data = dD))
report("D s(x,k=8)+(1|g) 40 levels", fD, NULL, NULL)
report("E s(x,k=8)+(1|g) 40 levels", fD, NULL, NA)
cat("DONE\n")
