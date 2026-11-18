# Lane drawdir, fix round 1: the facts for the user's decision on the
# reweighting where the per-direction reading fires on random-slope fits
# (review item 2). On the review's 8 random-slope fires (designs of
# dev/reviews/drawdir-review-survey.R) and the lane's 4 (S0.1 seeds 14,
# 16, 24, S0.2 seed 8; dev/drawdir-cases.R), on the fixed build:
# - every direction above the threshold, with its parameter, its
#   template component (theta = group-level sds, correlations, gp and
#   smooth scales; betad = a distributional parameter's coefficients),
#   and its excess;
# - in-sample MAE of the predictive mean against the true conditional
#   mean, and the new level's 90% and 98% coverage and width against
#   the known-parameter law, for the reweighted draws (lane) and the
#   normal law's (the verdict's fire set FALSE), 20000 draws, seed 1.
# No frm_sample() arm: the review's log has it for its 8.
# Usage: Rscript dev/drawdir-fires12.R <id> [...]
args <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-drawdir-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")), "\n")
ex <- parse("dev/reviews/drawdir-review-survey.R")
keep <- vapply(ex, function(e) {
  is.call(e) && identical(e[[1]], as.name("<-")) &&
    as.character(e[[2]]) %in% c("fq", "ilogit", "cum_mean", "cum_draw",
                                "design", "seeds_of", "said_of")
}, TRUE)
for (e in ex[keep]) eval(e, globalenv())
source("dev/drawdir-cases.R")
# the S designs' truth, rebuilt from their generator (case_fit())
s_truth <- function(kind, seed) {
  sds <- as.numeric(substring(kind, 2))
  set.seed(2000 + seed)
  G <- 15
  m <- 8
  g <- rep(seq_len(G), each = m)
  x <- stats::rnorm(G * m)
  ug <- stats::rnorm(G, 0, 0.8)
  sg <- stats::rnorm(G, 0, sds)
  list(mu = 1 + 0.5 * x + ug[g] + sg[g] * x,
       onew = function(n) stats::rnorm(n, 2, sqrt(0.8^2 + 4 * sds^2 + 1)))
}
# every direction's excess: draw_check_dirs() with its last steps
# replaced by a return of the whole vector
dirs_all <- local({
  f <- ns$draw_check_dirs
  b <- body(f)
  at <- which(vapply(as.list(b), function(e) {
    identical(deparse(e)[1], "i <- which.max(stat)")
  }, TRUE))
  b <- as.call(c(as.list(b)[seq_len(at - 1L)],
                 quote(list(stat = stat, R = R, k = k))))
  body(f) <- b
  f
})
real <- ns$draw_reweight
cap <- new.env()
grab <- function(object, dw) {
  dw$off <- length(dw$v0) - length(object$opt$par)
  cap$dw <- dw
  dw$D
}
cover <- function(pn, o) {
  q90 <- stats::quantile(pn, c(0.05, 0.95), names = FALSE)
  q98 <- stats::quantile(pn, c(0.01, 0.99), names = FALSE)
  sprintf("cov90 %.3f w90 %8.2f cov98 %.3f w98 %8.2f",
          mean(o >= q90[1] & o <= q90[2]), q90[2] - q90[1],
          mean(o >= q98[1] & o <= q98[2]), q98[2] - q98[1])
}
q <- function(expr) suppressMessages(suppressWarnings(expr))
for (id in args) {
  p <- strsplit(id, ":")[[1]]
  if (p[1] %in% c("S0.1", "S0.2")) {
    fit <- case_fit(id)$fit
    tr <- s_truth(p[1], as.integer(p[2]))
    r <- list(fit = fit, mu = tr$mu,
              nd = data.frame(x = 2, g = factor("new")), onew = tr$onew)
  } else {
    r <- design(p[1], as.integer(p[2]))
    fit <- r$fit
  }
  cat("==", id, "\n")
  cap$dw <- NULL
  assignInNamespace("draw_reweight", grab, ns = "frmtmb")
  set.seed(1)
  q(ns$predict_par_drawer(fit, TRUE, 2L))
  assignInNamespace("draw_reweight", real, ns = "frmtmb")
  dw <- cap$dw
  a <- dirs_all(fit, dw)
  comp <- ns$fit_draw_space(fit)$map$comp
  if (!is.null(a)) {
    j <- which(a$stat > ns$draw_dir_tol)
    for (i in j) {
      pj <- a$R[(i - 1L) %% a$k + 1L]
      cat(sprintf("  fires: %s%s (%s) %.3f\n", if (i > a$k) "-" else "+",
                  dw$labels[pj], comp[pj], a$stat[i]))
    }
    cat(sprintf("  directions read %d, finite %d\n", length(a$stat),
                sum(is.finite(a$stat))))
  }
  o <- NULL
  if (!is.null(r$onew)) {
    set.seed(99)
    o <- r$onew(2e5)
  }
  arm <- function(lab) {
    set.seed(1)
    s <- said_of(predict(fit, ndraws = 20000, summary = FALSE))
    chk <- fit$cache[["draw_check"]]
    line <- sprintf("  %-5s mae %.4f | %s", lab,
                    mean(abs(colMeans(s$val) - r$mu)), s$said)
    if (!is.null(o) && !is.null(r$nd)) {
      set.seed(1)
      pn <- q(predict(fit, newdata = r$nd, allow_new_levels = TRUE,
                      ndraws = 20000, summary = FALSE))[, 1]
      line <- paste(line, "|", cover(pn, o))
    }
    cat(line, sprintf("| fire %s dir %s (%s)\n", chk$fire,
                      format(chk$dir, digits = 4), chk$dir_par))
  }
  arm("lane")
  chk <- fit$cache[["draw_check"]]
  chk$fire <- FALSE
  assign("draw_check", chk, envir = fit$cache)
  arm("norm")
}
cat("DONE\n")
