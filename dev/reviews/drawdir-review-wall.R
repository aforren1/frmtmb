# Reviewer of lane drawdir, check 4: where a direction of the reading
# reads an infinite excess (a "joint wall", S2 = Inf), is it the data's
# wall or the objective failing far out? Along each direction whose
# excess is not finite, the objective's excess at fractions t of the
# move (the side law applied to t times the deviates), and at the end
# point each moved parameter alone.
# Usage: Rscript dev/reviews/drawdir-review-wall.R <id> [...]
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
  fit <- design(p[1], as.integer(p[2]))$fit
  assignInNamespace("draw_reweight", function(object, dw) {
    dw$off <- length(dw$v0) - length(object$opt$par)
    cap$dw <- dw
    dw$D
  }, ns = "frmtmb")
  invisible(fq(ns$predict_par_drawer(fit, TRUE, 2L)))
  assignInNamespace("draw_reweight", real, ns = "frmtmb")
  dw <- cap$dw
  lab <- dw$labels
  off <- dw$off
  np <- length(fit$opt$par)
  pos <- off + seq_len(np)
  Lf <- dw$law$Lc %||% dw$L
  sv <- colSums(Lf^2)
  R <- setdiff(pos[sv[pos] > 0], c(dw$held, dw$law$edge))
  k <- length(R)
  U <- Lf[, R, drop = FALSE] / rep(sqrt(sv[R]), each = nrow(Lf))
  Z <- cbind(2 * U, -2 * U)
  D <- dw$v0 + if (is.null(dw$law)) crossprod(dw$L, Z) else
    ns$draw_edge_apply(dw$law, dw$L, Z)
  e <- ns$draw_logratio(fit, dw, D, Z)$e
  obj <- fit$obj
  saved <- ns$obj_state_save(obj)
  f0 <- obj$fn(fit$opt$par)
  at <- function(v) {
    ns$obj_state_restore(obj, saved)
    tryCatch(suppressWarnings(obj$fn(v[pos])), error = function(e) NA) - f0
  }
  cat("==", id, "estimates:", paste(sprintf("%s=%.4g", lab, dw$v0),
                                     collapse = " "), "\n")
  for (i in which(!is.finite(e))) {
    j <- (i - 1L) %% k + 1L
    cat(sprintf("  direction %s %s: end point %s\n", lab[R[j]],
                if (i > k) "-" else "+",
                paste(sprintf("%s=%.4g", lab, D[, i]), collapse = " ")))
    for (t in c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)) {
      Zt <- t * Z[, i, drop = FALSE]
      Dt <- dw$v0 + if (is.null(dw$law)) crossprod(dw$L, Zt) else
        ns$draw_edge_apply(dw$law, dw$L, Zt)
      ld <- sum(Zt^2) / 2
      cat(sprintf("    t %.2f drop %10.4g law %7.4f | %s\n", t, at(Dt), ld,
                  paste(sprintf("%.4g", Dt[, 1]), collapse = " ")))
    }
    for (m in which(abs(D[, i] - dw$v0) > 1e-6)) {
      v <- dw$v0
      v[m] <- D[m, i]
      cat(sprintf("    %s alone at %.4g: drop %.4g\n", lab[m], v[m], at(v)))
    }
  }
  ns$obj_state_restore(obj, saved)
}
