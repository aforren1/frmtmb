# Item 1: is the estimate unbiased for m_g(theta) at theta AWAY from the anchor,
# with the anchor held fixed so that the seeds vary the draws and nothing else?
# Three independent parts:
#  (a) recover m and L from the draws themselves by regressing u on z, and
#      compare log|L| with the log-determinant the tape uses;
#  (b) the joint density the weights use against the plain joint objective at
#      moved u;
#  (c) mean over seeds of m_g_hat / m_g_exact, per group, where m_g_exact is
#      integrate() on base R densities. Should be 1 within its error.
# Usage: Rscript dev/review-impsamp-unbiased.R <design> [seeds] [draws]
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs()
suppressPackageStartupMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-designs.R")
args <- commandArgs(TRUE)
nm <- args[[1]]
S <- if (length(args) >= 2) as.integer(args[[2]]) else 200L
N <- if (length(args) >= 3) as.integer(args[[3]]) else 200L
de <- design_list[[nm]](7L)
ns <- asNamespace("frmtmb")

fr <- frm(de$formula, family = de$family, data = de$data, dry_run = "frame")
nll <- ns$build_objective(fr)
tpl <- fr[["par_template"]]
lap <- RTMB::MakeADFun(nll, tpl, random = "b", map = fr[["map"]], silent = TRUE)
opt <- suppressWarnings(nlminb(lap$par, lap$fn, lap$gr))
lay <- ns$imp_layout(fr[["re_blocks"]])
gmap <- ns$imp_group_map(fr, lay)
anchor <- setNames(opt$par, names(lap$par))
V <- tryCatch(solve(numDeriv::jacobian(lap$gr, opt$par)),
              error = function(e) NULL)
se <- if (is.null(V)) rep(0.2, length(anchor)) else sqrt(pmax(diag(V), 1e-6))
cat("design", nm, "anchor", format(anchor, digits = 4), "\n")
cat("se", format(se, digits = 3), "\n")

pts <- list(anchor = anchor)
k <- length(anchor)
pts$sd_up <- replace(anchor, k, anchor[k] + 0.7)
pts$sd_dn <- replace(anchor, k, anchor[k] - 0.7)
pts$b1_3se <- replace(anchor, 1L, anchor[1L] + 3 * se[1L])
pts$b1_m6se <- replace(anchor, 1L, anchor[1L] - 6 * se[1L])
if (k >= 2) pts$b2_6se <- replace(anchor, 2L, anchor[2L] + 6 * se[2L])
pts$joint <- anchor + 4 * se * rep_len(c(1, -1), k)

gr <- split(seq_len(nrow(de$data)), de$data$g)
dl <- lapply(gr, function(i) de$data[i, , drop = FALSE])
exact_g <- function(p) {
  vapply(dl, function(d) {
    exact_group_1d(function(u) de$ll(unname(p), d, u), de$sdof(unname(p)))
  }, 0)
}

# one plan per seed; (a) and (b) at seed 1 only
res <- list()
for (pn in names(pts)) {
  p <- pts[[pn]]
  eg <- exact_g(p)
  ratios <- matrix(NA_real_, length(eg), S)
  for (s in seq_len(S)) {
    plan <- ns$imp_plan(lap, fr, lay, anchor, N, s)
    io <- ns$build_importance_objective(fr, lay, gmap, plan)
    pl <- ns$imp_par_list(tpl, p, fr[["map"]])
    pe <- environment(io$fn)
    a <- as.matrix(io$amat(pl))
    pr <- pe$proposal(pl)
    ld <- as.numeric(pr$logdet)
    lm <- apply(a, 1L, log_mean_exp) + 0.5 * plan$q * log(2 * pi) + ld
    ratios[, s] <- exp(lm - eg)
    if (s == 1L) {
      # (a) m and L read from the draws
      z <- plan$z
      ng <- plan$n_group
      nd <- plan$n_draw
      um <- do.call(rbind, lapply(seq_len(ng), function(g) {
        vapply(pr$u, function(uj) {
          as.numeric(uj)[(seq_len(nd) - 1L) * ng + g]
        }, numeric(nd))
      }))
      if (plan$q == 1L) {
        slope <- vapply(seq_len(ng), function(g) {
          ug <- as.numeric(pr$u[[1L]])[(seq_len(nd) - 1L) * ng + g]
          zg <- z[plan$idx[1L, g], ]
          unname(coef(lm(ug ~ zg))[2L])
        }, 0)
        err_ld <- max(abs(log(abs(slope)) - ld))
      } else err_ld <- NA
      # (b) the weights' joint density at the moved u, against the plain
      # joint objective. Sum over groups at two draw columns.
      err_joint <- 0
      for (j in c(1L, nd)) {
        joint <- pl
        b <- numeric(length(plan$uhat))
        for (qq in seq_len(plan$q)) {
          b[plan$idx[qq, ]] <- as.numeric(pr$u[[qq]])[(j - 1L) * ng + seq_len(ng)]
        }
        joint[["b"]] <- b
        ours <- sum(a[, j] - plan$zsq[, j])
        err_joint <- max(err_joint, abs(ours + nll(joint)))
      }
    }
  }
  mr <- rowMeans(ratios)
  se_r <- apply(ratios, 1L, sd) / sqrt(S)
  z_r <- (mr - 1) / se_r
  tot <- apply(ratios, 2L, prod)
  res[[pn]] <- data.frame(
    design = nm, pt = pn, S = S, N = N, err_logdet = err_ld,
    err_joint = err_joint,
    mean_ratio_min = min(mr), mean_ratio_max = max(mr),
    max_abs_z = max(abs(z_r)), mean_z = mean(z_r),
    total_mean_ratio = mean(tot), total_se = sd(tot) / sqrt(S),
    sum_group_var = sum(apply(log(ratios), 1L, var)),
    nll_exact = -sum(eg))
  print(format(res[[pn]][-(1:1)], digits = 4))
}
res <- do.call(rbind, res)
out <- "C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-out/unbiased"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.table(res, file.path(out, paste0(nm, ".tsv")), sep = "\t", quote = FALSE,
            row.names = FALSE)
cat("rate range:", format(range(plan$hlik_rate), digits = 4),
    " hlik range:", format(range(plan$hlik), digits = 4), "\n")
