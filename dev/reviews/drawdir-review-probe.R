# Reviewer of lane drawdir: each direction of the per-direction reading,
# one fit at a time. For a design id of dev/reviews/drawdir-review-survey.R:
# the law (captured from predict_par_drawer()), each parameter's direction
# at +-draw_edge_k through the side law: eA (the direction's excess), the
# slice's excess eB at draw_edge_k conditional sds, S2 = eA - max(eB, 0),
# the direction's move in each parameter, and the verdict. With a
# `sigid` design (identity-link sigma ~ x), the natural sigma at x = 0 and
# x = 1 at each direction's point (item 4: can the regression move leave
# the support?).
# Usage: Rscript dev/reviews/drawdir-review-probe.R <id> [...]
args <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-drawdir-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
ex <- parse("dev/reviews/drawdir-review-survey.R")
keep <- vapply(ex, function(e) {
  is.call(e) && identical(e[[1]], as.name("<-")) &&
    as.character(e[[2]]) %in% c("fq", "ilogit", "cum_mean", "cum_draw",
                                "design", "seeds_of")
}, TRUE)
for (e in ex[keep]) eval(e, globalenv())
ns <- asNamespace("frmtmb")
real <- ns$draw_reweight
cap <- new.env()
for (id in args) {
  p <- strsplit(id, ":")[[1]]
  r <- design(p[1], as.integer(p[2]))
  fit <- r$fit
  assignInNamespace("draw_reweight", function(object, dw) {
    dw$off <- length(dw$v0) - length(object$opt$par)
    cap$dw <- dw
    dw$D
  }, ns = "frmtmb")
  invisible(fq(ns$predict_par_drawer(fit, TRUE, 2L)))
  assignInNamespace("draw_reweight", real, ns = "frmtmb")
  dw <- cap$dw
  lab <- dw$labels
  cat("==", id, "par:", paste(sprintf("%s=%.4g", lab, dw$v0),
                              collapse = " "), "\n")
  cat("  side law:", if (is.null(dw$law)) "none" else paste(
    "edge", paste(lab[dw$law$edge], collapse = ","), "narrowed",
    paste(lab[dw$law$narrowed], collapse = ",")), "\n")
  off <- dw$off
  np <- length(fit$opt$par)
  pos <- off + seq_len(np)
  Lf <- dw$law$Lc %||% dw$L
  sv <- colSums(Lf^2)
  R <- setdiff(pos[sv[pos] > 0], c(dw$held, dw$law$edge))
  k <- length(R)
  B <- Lf[, R, drop = FALSE]
  U <- B / rep(sqrt(sv[R]), each = nrow(B))
  ck <- ns$draw_edge_k
  Z <- cbind(ck * U, -ck * U)
  D <- dw$v0 + if (is.null(dw$law)) crossprod(dw$L, Z) else
    ns$draw_edge_apply(dw$law, dw$L, Z)
  r2 <- ns$draw_logratio(fit, dw, D, Z)
  csd <- ns$draw_dirs_csd(B)
  obj <- fit$obj
  saved <- ns$obj_state_save(obj)
  f0 <- obj$fn(fit$opt$par)
  for (i in seq_len(2 * k)) {
    j <- (i - 1L) %% k + 1L
    q <- fit$opt$par
    q[R[j] - off] <- q[R[j] - off] + (if (i > k) -1 else 1) * ck * csd[j]
    ns$obj_state_restore(obj, saved)
    f1 <- tryCatch(suppressWarnings(obj$fn(q)), error = function(e) NA)
    eB <- f1 - f0 - ck^2 / 2
    mv <- D[, i] - dw$v0
    extra <- ""
    if (p[1] == "sigid") {
      si <- match("sigma_(Intercept)", lab)
      sx <- match("sigma_x", lab)
      extra <- sprintf(" | sigma(x=0) %.3f sigma(x=1) %.3f", D[si, i],
                       D[si, i] + D[sx, i])
    }
    cat(sprintf("  %-28s %s eA %9.3f eB %9.3f S2 %9.3f sratio %.3f | move %s%s\n",
                lab[R[j]], if (i > k) "-" else "+", r2$e[i], eB,
                r2$e[i] - max(eB, 0), sqrt(sv[R[j]]) / csd[j],
                paste(sprintf("%.3g", mv), collapse = " "), extra))
  }
  ns$obj_state_restore(obj, saved)
  set.seed(5)
  pr <- fq(predict(fit, ndraws = 200, summary = FALSE))
  chk <- fit$cache[["draw_check"]]
  cat(sprintf("  verdict fire %s excess %.3f rel %.3f dir %.3f (%s) n %d\n",
              chk$fire, chk$excess, chk$rel, chk$dir, chk$dir_par,
              chk$dir_n))
  if (p[1] == "sigid") {
    # how many of the normal law's draws put sigma at or below 0 inside
    # the data's x range, and what the check's own test draws read
    set.seed(1)
    Zn <- matrix(stats::rnorm(length(dw$v0) * 4000), length(dw$v0))
    Dn <- dw$v0 + crossprod(dw$L, Zn)
    si <- match("sigma_(Intercept)", lab)
    sx <- match("sigma_x", lab)
    lo <- pmin(Dn[si, ], Dn[si, ] + Dn[sx, ])
    cat(sprintf("  normal-law draws with sigma <= 0 on [0, 1]: %d of 4000\n",
                sum(lo <= 0)))
  }
}
