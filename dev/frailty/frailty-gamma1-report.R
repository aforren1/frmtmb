source("frailty-common.R")
rd <- function(p) {
  x <- readLines(p)
  rows <- lapply(strsplit(x, "\t"), function(f) {
    kv <- strsplit(f, "=")
    v <- vapply(kv, function(z) trimws(z[2L]), character(1))
    names(v) <- vapply(kv, `[`, character(1), 1L)
    v
  })
  nm <- names(rows[[which.max(lengths(rows))]])
  d <- as.data.frame(do.call(rbind, lapply(rows, function(r) r[nm])),
                     stringsAsFactors = FALSE)
  for (j in nm) {
    if (j == "status") next
    d[[j]] <- suppressWarnings(as.numeric(d[[j]]))
  }
  d
}
d <- rd("frailty-gamma1.tsv")
cat("replicates", nrow(d), " ok", sum(d$status == "ok"), "\n")
d <- d[d$status == "ok", ]
n <- nrow(d)
mc <- function(k) sqrt(k / n * (1 - k / n) / n)
f3 <- function(v) formatC(v, digits = 4, format = "g")
cat("\ntruth: sd_u 0.2, gamma1 1.3, gamma0 -2.092269, beta 0.6\n\n")
cat(sprintf("%-18s %10s %10s %10s\n", "quantity", "mean", "sd", "cover"))
row <- function(lab, v, cov = NULL) {
  cat(sprintf("%-18s %10s %10s %10s\n", lab, f3(mean(v)), f3(sd(v)),
              if (is.null(cov)) "" else
                paste0(sum(cov), "/", n, " ", f3(mean(cov)))))
}
row("sd(gamma1|centre)", d$sd_hat, d$cov_sd)
row("gamma1", d$g1_on, d$cov_g1)
row("gamma0", d$g0_on, d$cov_g0)
row("beta(trt)", d$beta_on, d$cov_beta)
d$err_drop <- vapply(seq_len(n), function(i) {
  tr <- attr(g1_sim(d$seed[i], sd_u = 0.2, gamma1 = 1.3), "truth")
  mean(abs(d$g1_on[i] - tr$shape))
}, numeric(1))
cat("\nbinomial mcse at 0.95 with", n, "replicates:",
    f3(sqrt(0.95 * 0.05 / n)), "\n")
cat("\nper-centre shape error, mean |shape_hat - shape_true|:\n")
cat("  component ON  ", f3(mean(d$err_on)), " sd ", f3(sd(d$err_on)),
    "\n  component OFF ", f3(mean(d$err_off)), " sd ",
    f3(sd(d$err_off)), "\n  ratio off/on ",
    f3(mean(d$err_off) / mean(d$err_on)),
    "   ON better on ", sum(d$err_on < d$err_off), " of ", n, "\n",
    sep = "")
cat("  same fit with its deviations DROPPED ", f3(mean(d$err_drop)),
    "  ratio drop/on ", f3(mean(d$err_drop) / mean(d$err_on)),
    "  ON better on ", sum(d$err_on < d$err_drop), " of ", n, "
", sep = "")
cat("  correlation with the truth: mean ", f3(mean(d$cor_on)),
    " min ", f3(min(d$cor_on)), "\n", sep = "")
cat("\nlog likelihood ON - OFF: mean ", f3(mean(d$ll_on - d$ll_off)),
    " min ", f3(min(d$ll_on - d$ll_off)), "\n", sep = "")
cat("Laplace error against the exact integral: mean ",
    f3(mean(d$lap_err)), " range ",
    f3(min(d$lap_err)), " to ", f3(max(d$lap_err)), "\n", sep = "")
cat("\nmonotonicity: n_nonmono max ", max(d$n_nonmono),
    "; smallest fitted per-centre slope ", f3(min(d$min_slope)),
    " (mean ", f3(mean(d$min_slope)), ")\n", sep = "")
cat("margin to the floor, min slope / sd_hat: min ",
    f3(min(d$margin_sd)), " mean ", f3(mean(d$margin_sd)), "\n", sep = "")
cat("\nconvergence: conv=0 on ", sum(d$conv == 0), "/", n,
    ", pdHess on ", sum(d$pdHess == 1), "/", n,
    ", max maxgrad ", f3(max(d$maxgrad)), "\n", sep = "")
