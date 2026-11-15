# Review of lane refits: exact references for influence() deletion refits
# of s(x, g, bs = "fs") and t2(x, v) (the fit's basis written out), the
# gp() deletion theta, and 3188 per-row predictions against nnet.
source("C:/Users/adf44/source/r/frmtmb-wt-refits/dev/reviews/refits-review-lib.R")
set.seed(8)
n <- 60
dl <- data.frame(x = runif(n, 0, 3), z = rnorm(n), v = runif(n),
                 f = factor(rep(c("a", "b"), n / 2)),
                 g = factor(rep(1:5, 12)), blk = factor(rep(1:6, 10)))
dl$y <- sin(dl$x) + 0.3 * dl$z + dl$v^2 + rnorm(n, 0, 0.3)

written <- function(fit, keep) {
  si <- fit$frame$linpreds[[1L]]$smooths[[1L]]
  M <- frmtmb:::smooth_basis_at(si, keep)
  at <- 0L
  terms <- character(0)
  for (b in seq_along(si$nr)) {
    nm <- paste0("r", b, "_", seq_len(si$nr[b]))
    for (k in seq_along(nm)) keep[[nm[k]]] <- M[, at + k]
    at <- at + si$nr[b]
    terms <- c(terms, paste0("homdiag(0 + ", paste(nm, collapse = " + "),
                             " | id)"))
  }
  fx <- character(0)
  if (si$nf > 0) for (k in seq_len(si$nf)) {
    keep[[paste0("fx", k)]] <- M[, at + k]; fx <- c(fx, paste0("fx", k))
  }
  keep$id <- factor(1)
  fo <- as.formula(paste("y ~ z", paste(c(fx, terms), collapse = " + "),
                         sep = " + "))
  cat("written formula:", deparse1(fo), "\n")
  cat("nr:", si$nr, " nf:", si$nf, "\n")
  suppressMessages(frm(fo, data = keep))
}
for (form in list(y ~ z + s(x, g, bs = "fs", k = 4), y ~ z + t2(x, v, k = 4),
                  y ~ z + s(x, k = 5))) {
  cat("\n====", deparse1(form), "====\n")
  fit <- suppressMessages(frm(form, data = dl))
  inf <- suppressWarnings(influence(fit, groups = "blk", data = dl))
  keep <- dl[dl$blk != "1", ]
  ref <- go("written refit", written(fit, keep))
  if (is.character(ref)) next
  cat("table fixed row 1:", sprintf("%.6f", inf$fixed["1", ]), "\n")
  cat("written fixef    :", sprintf("%.6f", fixef(ref, flatten = TRUE)), "\n")
  cat("table theta row 1:", sprintf("%.6f", inf$theta["1", ]), "\n")
  cat("written theta    :", sprintf("%.6f", ref$estimates$theta), "\n")
  cat("full theta       :", sprintf("%.6f", inf$theta_full), "\n")
  cat("written logLik   :", sprintf("%.8f", as.numeric(logLik(ref))), "\n")
}

cat("\n==== gp exact deletion of the max-x row: theta scale ====\n")
fg <- suppressMessages(frm(y ~ z + gp(x), data = dl))
imax <- which.max(dl$x)
dl$one <- factor(seq_len(n) == imax)
ig <- suppressWarnings(influence(fg, groups = "one", data = dl))
man <- suppressMessages(frm(y ~ z + gp(x), data = dl[-imax, ]))
cat("dmax full:", fg$frame$re_blocks[[1]]$gp_lscale_div %||% NA,
    " dmax manual:", man$frame$re_blocks[[1]]$gp_lscale_div %||% NA, "\n")
cat("table theta:", ig$theta["TRUE", ], "\n")
cat("manual theta:", man$estimates$theta, "\n")
cat("full theta:", ig$theta_full, "\n")
go("summary lscale manual vs full",
   list(manual = capture.output(summary(man))[grep("lscale", capture.output(summary(man)))],
        full = capture.output(summary(fg))[grep("lscale", capture.output(summary(fg)))]))

cat("\n==== 3188 per-row predictions against nnet::multinom ====\n")
set.seed(5)
dm <- data.frame(x = rnorm(200), N = 10L)
pr <- cbind(1, exp(0.3 + 0.5 * dm$x), exp(-0.2 * dm$x))
pr <- pr / rowSums(pr)
dm$Y <- t(sapply(1:200, function(i) rmultinom(1, 10, pr[i, ])))
colnames(dm$Y) <- c("a", "b", "c")
fm <- frm(Y | trials(N) ~ x, family = multinomial(), data = dm)
ds <- data.frame(x = dm$x, Ya = dm$Y[, 1], Yb = dm$Y[, 2], Yc = dm$Y[, 3])
mn <- nnet::multinom(cbind(Ya, Yb, Yc) ~ x, data = ds, trace = FALSE)
go("per row", {
  pd <- as.data.frame(marginaleffects::predictions(fm))
  pn <- as.data.frame(marginaleffects::predictions(mn, type = "probs"))
  pn$group <- c(Ya = "a", Yb = "b", Yc = "c")[as.character(pn$group)]
  of <- order(pd$group, pd$rowid); on <- order(pn$group, pn$rowid)
  c(rows = nrow(pd), rows_nn = nrow(pn),
    max_rel_est = max(abs(pd$estimate[of] / 10 - pn$estimate[on]) / pn$estimate[on]),
    max_rel_se = max(abs(pd$std.error[of] / 10 - pn$std.error[on]) / pn$std.error[on]))
})
go("confint_varcorr rownames of a plain (0 + f || g) fit",
   { dd <- dl; dd$f4 <- factor(rep(letters[1:4], 15))
     rownames(confint_varcorr(frm(y ~ z + (0 + f4 || g), data = dd))) })
