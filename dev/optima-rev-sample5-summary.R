# Reviewer of lane optima: summarize the simplex draws written by
# optima-rev-sample2-arm.R, -sample3-brms.R and -sample4-dirichlet.R.
#   Rscript dev/optima-rev-sample5-summary.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
out <- "C:/Users/adf44/source/r/frmtmb-wt-optima/dev/optima-rev-out"
D <- 3
# which sign sheet of the chart a draw is on: the signs of u
sheets <- function(Z, chain) {
  r2 <- rowSums(Z * Z)
  U <- (2 * Z %*% t(Bn) + (r2 - 1) * (-1 / sqrt(D))) / (r2 + 1)
  lab <- apply(sign(U), 1, function(s) paste(ifelse(s > 0, "+", "-"), collapse = ""))
  print(table(chain = chain, sheet = lab))
}
Bn <- local({ H <- stats::contr.helmert(D); sweep(H, 2, sqrt(colSums(H^2)), "/") })
chart <- function(Z) {
  r2 <- rowSums(Z * Z)
  U <- (2 * Z %*% t(Bn) + (r2 - 1) * (-1 / sqrt(D))) / (r2 + 1)
  W <- U * U; W / rowSums(W)
}
softmax <- function(Z) {
  A <- cbind(0, Z); A <- A - apply(A, 1, max)
  E <- exp(A); E / rowSums(E)
}
summ <- function(label, W, chain, Z = NULL, rhat = NULL, div = NULL) {
  mn <- apply(W, 1, min)
  db <- sqrt(rowSums((W - 1 / D)^2))
  q <- function(v) paste(sprintf("%.3g", quantile(v, c(.05, .5, .95))), collapse = "/")
  cat(sprintf("%-22s mean w %s | sd w %s\n", label,
              paste(sprintf("%.3f", colMeans(W)), collapse = " "),
              paste(sprintf("%.3f", apply(W, 2, sd)), collapse = " ")))
  cat(sprintf("%-22s min(w) q05/50/95 %s ; P(min<.01) %.3f | dist-bary q05/50/95 %s ; P(dist<.05) %.3f\n",
              "", q(mn), mean(mn < 0.01), q(db), mean(db < 0.05)))
  for (ch in sort(unique(chain))) {
    k <- chain == ch
    cat(sprintf("%-22s   chain %d: mean w %s ; median min(w) %.3g ; median dist %.3g\n",
                "", ch, paste(sprintf("%.3f", colMeans(W[k, ])), collapse = " "),
                median(mn[k]), median(db[k])))
  }
  if (!is.null(Z)) {
    cat(sprintf("%-22s zeta range [%.4g, %.4g], max|zeta| %.4g, median|zeta| (%s)\n",
                "", min(Z), max(Z), max(abs(Z)),
                paste(sprintf("%.3g", apply(abs(Z), 2, median)), collapse = ", ")))
  }
  nc <- length(unique(chain))
  rw <- sapply(seq_len(ncol(W)), function(j) {
    rstan::Rhat(matrix(W[, j], ncol = nc))
  })
  cat(sprintf("%-22s Rhat of w (rstan::Rhat) %s\n", "",
              paste(sprintf("%.3f", rw), collapse = " ")))
  if (!is.null(rhat)) cat(sprintf("%-22s Rhat zeta %s ; divergences %s\n", "",
                                  paste(sprintf("%.3f", rhat), collapse = " "),
                                  paste(div, collapse = "+")))
}
from_arr <- function(r, map) {
  a <- r$arr
  zc <- grep("^zeta", dimnames(a)[[3]])
  Z <- do.call(rbind, lapply(seq_len(dim(a)[2]), function(ch) a[, ch, zc]))
  ch <- rep(seq_len(dim(a)[2]), each = dim(a)[1])
  list(Z = Z, W = map(Z), chain = ch,
       rhat = r$summary[grep("^zeta", rownames(r$summary)), "Rhat"])
}
for (case in c("weak", "strong")) {
  cat("\n================", case, "================\n")
  for (arm in c("base", "lane")) {
    f <- file.path(out, paste0("draws-", arm, ".rds"))
    if (!file.exists(f)) next
    r <- readRDS(f)[[case]]
    x <- from_arr(r, if (arm == "lane") chart else softmax)
    cat(sprintf("ML zeta (%s): %s ; ML w %s\n", arm,
                paste(sprintf("%.4g", r$ml_zeta), collapse = ", "),
                paste(sprintf("%.4g", (if (arm == "lane") chart else softmax)(
                  matrix(r$ml_zeta, 1))), collapse = ", ")))
    summ(paste(arm, "frm_sample"), x$W, x$chain, x$Z, x$rhat, r$div)
    if (arm == "lane") sheets(x$Z, x$chain)
  }
  f <- file.path(out, "draws-brms.rds")
  if (file.exists(f)) {
    r <- readRDS(f)[[case]]
    sc <- grep("^simo_", colnames(r$draws))
    W <- r$draws[, sc]
    ch <- rep(1:2, each = nrow(W) / 2)
    rh <- r$summary[grep("^simo_", rownames(r$summary)), "Rhat"]
    cat("brms columns:", colnames(r$draws)[sc], "\n")
    summ("brms dirichlet(1)", W, ch, NULL, rh, r$div)
  }
  for (arm in c("base", "lane")) {
    f <- file.path(out, paste0("draws-dir-", arm, ".rds"))
    if (!file.exists(f)) next
    r <- readRDS(f)[[case]]
    x <- from_arr(r, if (arm == "lane") chart else softmax)
    summ(paste(arm, "+dir(1)+logJ"), x$W, x$chain, x$Z, x$rhat, r$div)
    if (arm == "lane") sheets(x$Z, x$chain)
  }
}
