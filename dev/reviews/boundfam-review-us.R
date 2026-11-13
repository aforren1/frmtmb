# Review of lane boundfam: the us() row-edge NA rule (3098) in dimension
# 3 and 4, exactly / nearly rank deficient, on constructed thetas and on
# fits; base vs lane intervals on the same fits.
#
#   Rscript dev/reviews/boundfam-review-us.R <lane|base>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-boundfam-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
mcor <- function(C) {
  # multiple correlation of each variable with the ones before it
  d <- nrow(C)
  c(NA, vapply(2:d, function(i) {
    a <- C[i, 1:(i - 1)]
    B <- C[1:(i - 1), 1:(i - 1), drop = FALSE]
    sqrt(drop(t(a) %*% solve(B, a)))
  }, 0))
}
if (arm == "lane") {
  cat("\n== constructed thetas\n")
  mk <- function(d, i, R) {
    # variable i = combination of the earlier ones with multiple corr R
    set.seed(10 + d + i)
    A <- matrix(rnorm(d * d), d)
    A[i, ] <- 0
    w <- rnorm(i - 1)
    A[i, ] <- R * colSums(w * A[1:(i - 1), , drop = FALSE]) /
      sqrt(sum((colSums(w * A[1:(i - 1), , drop = FALSE]))^2)) *
      sqrt(sum(A[1, ]^2))
    S <- crossprod(t(A))
    S
  }
  for (d in 3:4) for (R in c(0.999, 0.99985, 0.99995, 1 - 1e-9)) {
    for (i in 3:d) {
      # build a correlation matrix whose row i has multiple corr R exactly
      set.seed(100 * d + i)
      C <- diag(d)
      L <- matrix(rnorm(d * d) * 0.4, d); L[upper.tri(L, TRUE)] <- 0
      diag(L) <- 1
      # make row i: off-diagonal part scaled to S = R^2 / (1 - R^2)
      off <- L[i, seq_len(i - 1)]
      off <- off / sqrt(sum(off^2)) * sqrt(R^2 / (1 - R^2))
      L[i, seq_len(i - 1)] <- off
      th <- c(rep(0, d), L[lower.tri(L)])
      Cm <- ns$us_chol_cor(L[lower.tri(L)], d)
      e <- ns$varcorr_us_row_edge(th, d)
      pr <- which(lower.tri(diag(d)), arr.ind = TRUE)
      cat(sprintf("d=%d row=%d R=%.10f mcor(row)=%.10f min eig=%.3g flagged=%s\n",
                  d, i, R, mcor(Cm)[i], min(eigen(Cm)$values),
                  paste(apply(pr[e, , drop = FALSE], 1, paste, collapse = ","), collapse = " ")))
    }
  }
}

cat("\n== fits\n")
one <- function(lab, d, form) {
  fit <- suppressWarnings(frm(form, data = d))
  th <- fit$estimates$theta
  tr <- suppressWarnings(ns$varcorr_trans_rows(fit))
  bk <- fit$frame$re_blocks[[1]]
  dd <- length(bk$cnms %||% bk[["cnms"]])
  C <- tryCatch(ns$us_chol_cor(th[-seq_len(dd)], dd), error = function(e) NULL)
  cat(sprintf("%s code=%s obj=%.6f mcor=%s\n", lab, fit$opt$convergence,
              fit$opt$objective,
              if (is.null(C)) "?" else paste(sprintf("%.7f", mcor(C)[-1]), collapse = ",")))
  for (i in seq_len(nrow(tr))) {
    cat(sprintf("   %-28s %-4s est_t=%12.6g se_t=%12.6g\n", tr$term[i],
                tr$type[i], tr$est_t[i] %||% NA, tr$se_t[i]))
  }
}
for (s in 1:6) {
  ng <- 40
  set.seed(s)
  d <- data.frame(x = rnorm(ng * 10), z = rnorm(ng * 10), w = rnorm(ng * 10),
                  g = factor(rep(seq_len(ng), each = 10)))
  u0 <- rnorm(ng, 0, 0.8); u1 <- rnorm(ng, 0, 0.5); u2 <- rnorm(ng, 0, 0.6)
  # dimension 4, rank 3: the w slope is a combination of the others
  u3 <- 0.5 * u0 - 0.4 * u1 + 0.3 * u2
  d$y4 <- 1 + 0.5 * d$x - 0.3 * d$z + 0.2 * d$w + u0[d$g] + u1[d$g] * d$x +
    u2[d$g] * d$z + u3[d$g] * d$w + rnorm(ng * 10, 0, 0.5)
  d$y <- d$y4
  one(paste0("rank3 d4 s", s), d, y ~ x + z + w + (1 + x + z + w | g))
  # dimension 3, nearly rank deficient: true multiple corr 0.999
  e <- rnorm(ng)
  u2n <- 0.4 * u0 - 0.6 * u1
  u2n <- u2n + e * sd(u2n) * sqrt(1 / 0.999^2 - 1)
  d$y <- 1 + 0.5 * d$x - 0.3 * d$z + u0[d$g] + u1[d$g] * d$x +
    u2n[d$g] * d$z + rnorm(ng * 10, 0, 0.5)
  one(paste0("near d3 s", s), d, y ~ x + z + (1 + x + z | g))
}
