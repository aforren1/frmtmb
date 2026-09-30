# Reviewer attack (lane ceplot) on the renamed-level copy of
# R/ce-levels.R, on draws, with no brms: per posterior draw, the curve
# at the grid row minus every OBSERVED part of the predictor (read from
# the same draw) is the new level's effect. For an intercept-only new
# block it is constant in x and, whitened by that draw's own block sd,
# standard normal; for (1 + x | g:h) it is linear in x and whitened by
# that draw's 2x2 covariance. Shapes: dev/ceplot-rev-shapes.R.
#   Rscript dev/ceplot-rev-crossed.R > dev/ceplot-rev-log/crossed.txt
# Seeds: data per shape (49, 51, 52, 53, 54, 55), draws 1, curves 7.
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
source("C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-rev-shapes.R")
cat("frmtmb from", find.package("frmtmb"), "\n")
say <- function(...) cat(sprintf(...), "\n", sep = "")
f4 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
N <- 400
S <- shape_fits()
# the observed part of each shape at its condition rows, per draw, and
# the new blocks (by group_name) whose effects the residual holds
obs_part <- list(
  A = list(c("r_g[1,Intercept]", "r_h[1,Intercept]")),
  B = list(c("r_g[1,Intercept]", "r_h[1,Intercept]"),
           c("r_g[2,Intercept]", "r_h[3,Intercept]"),
           c("r_g[3,Intercept]", "r_h[2,Intercept]", "r_g:h[3_2,Intercept]")),
  D = list(c("r_g[1,Intercept]", "r_h[1,Intercept]")),
  E = list(NULL),
  F = list(c("r_h[1,Intercept]")),
  G = list(c("r_g[1,Intercept]", "r_h[1,Intercept]", "r_k[1,Intercept]")),
  H = list(c("r_g[1,Intercept]", "r_h[1,Intercept]")))
new_blocks <- list(A = list("g:h"), B = list("g:h", "g:h", character(0)),
                   D = list("h:g"), E = list("g:h"), F = list("g1:h"),
                   G = list(c("g:h", "g:k")), H = list("g:h"))
block_cov <- function(fit, th, gname) {
  S <- 0
  for (bk in fit$frame$re_blocks) {
    if (!identical(bk$group_name, gname)) next
    S <- S + as.matrix(frmtmb:::covstruct_registry[[bk$covstruct]]$vcov(
      th[bk$theta_idx], bk))
  }
  S
}
for (nm in names(S)) {
  fit <- S[[nm]]$fit
  cat("\n== shape", nm, ":", vapply(fit$frame$re_blocks, `[[`, "",
                                    "term_label"), "\n")
  before <- serialize(list(fit$frame, fit$estimates), NULL)
  ds <- hand(fit, N)
  M <- ds$draws
  if (nm %in% c("E", "F")) {
    cat("labels:", grep("^r_", colnames(M), value = TRUE)[1:8], "\n")
  }
  if (nm == "E") {
    lab <- grep("^r_g\\[1,", colnames(M), value = TRUE)
    obs_part$E <- list(c(lab, "r_h[1,Intercept]"))
  }
  if (nm == "F") {
    lab <- grep("^r_mm", colnames(M), value = TRUE)[1:2]
    cat("mm labels used (members 1 and 2):", lab, "\n")
  }
  cond <- S[[nm]]$cond
  ce <- tryCatch(conditional_effects(ds, "x", resolution = 3,
                                     re_formula = NULL, conditions = cond,
                                     spaghetti = TRUE, seed = 7),
                 error = function(e) e)
  if (inherits(ce, "error")) {
    say("  ERROR %s", conditionMessage(ce))
    next
  }
  df <- ce$x
  sp <- attr(df, "spaghetti")
  conds <- levels(factor(sp$cond__))
  th_cols <- grep("^theta_", colnames(M), value = TRUE)
  zs <- list()
  for (ci in seq_along(conds)) {
    s <- sp[sp$cond__ == conds[ci], ]
    s <- s[order(s$sample__, s$x), ]
    xs <- sort(unique(s$x))
    E <- matrix(s$estimate__, ncol = length(xs), byrow = TRUE)
    stopifnot(nrow(E) == N)
    known <- M[, "b_Intercept"] + outer(M[, "b_x"], xs)
    op <- obs_part[[nm]][[min(ci, length(obs_part[[nm]]))]]
    for (l in op) known <- known + M[, l]
    if (nm == "F") {
      known <- known + 0.5 * (M[, "r_mmg1g2[1,Intercept]"] +
                                M[, "r_mmg1g2[2,Intercept]"])
    }
    R <- E - known
    nb <- new_blocks[[nm]][[min(ci, length(new_blocks[[nm]]))]]
    grow <- df[df$cond__ == conds[ci], ][1, ]
    if (!length(nb)) {
      say("  cond %s (observed): max |residual| %.3g", conds[ci],
          max(abs(R)))
      next
    }
    if (nm == "H") {
      u1 <- (R[, 3] - R[, 1]) / (xs[3] - xs[1])
      u0 <- R[, 1] - u1 * xs[1]
      lin <- max(abs(R[, 2] - (u0 + u1 * xs[2])))
      z <- t(vapply(seq_len(N), function(i) {
        Sg <- block_cov(fit, M[i, th_cols], nb)
        backsolve(chol(Sg), c(u0[i], u1[i]), transpose = TRUE)
      }, numeric(2)))
      say("  cond %s: linear-in-x residual, max dev %.3g", conds[ci], lin)
      say("  cond %s: whitened (u0, u1): means %s, var %s, cor %.3f",
          conds[ci], f4(colMeans(z)), f4(apply(z, 2, var)), cor(z)[1, 2])
      zs[[ci]] <- z[, 1]
      next
    }
    flat <- max(abs(R - R[, 1]))
    sdv <- vapply(seq_len(N), function(i) {
      sqrt(sum(vapply(nb, function(b) block_cov(fit, M[i, th_cols], b), 1)))
    }, 1)
    z <- R[, 1] / sdv
    zs[[ci]] <- z
    say("  cond %s (g/h as in grid: %s): constant in x to %.3g; z mean %.4f var %.4f (N = %d, var SE %.3f); KS p %.3f",
        conds[ci], paste(unlist(grow[intersect(names(grow), c("g", "h", "k", "g1", "g2", "f"))]), collapse = ","),
        flat, mean(z), var(z), N, sqrt(2 / (N - 1)),
        suppressWarnings(stats::ks.test(z, "pnorm")$p.value))
  }
  if (length(zs) >= 2 && !is.null(zs[[1]]) && !is.null(zs[[2]])) {
    say("  cor of the two new combinations' z: %.4f", cor(zs[[1]], zs[[2]]))
  }
  after <- serialize(list(fit$frame, fit$estimates), NULL)
  say("  fit frame and estimates unchanged after the call: %s",
      identical(before, after))
}
say("done")
