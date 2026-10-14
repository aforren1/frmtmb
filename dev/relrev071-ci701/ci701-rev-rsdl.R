# (ci701-rev-rsdl.R) the log-likelihood change se_at_edge() reads on the
# rs20 fits (seeds 1..20): margin against se_edge_tol = 1e-6
.libPaths(c("C:/Users/adf44/source/r/rellib-r9", "C:/Users/adf44/source/r/rellib-r8", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
for (s in 1:20) {
  set.seed(s)
  g <- factor(rep(seq_len(20), each = 10))
  d <- data.frame(g = g, x = stats::rnorm(200))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(20, 0, 0.7)[g] + stats::rnorm(200)
  fit <- suppressMessages(suppressWarnings(frm(bf(y ~ x + (1 + x | g)), family = gaussian(), data = d)))
  p <- fit$opt$par; j <- 6
  q <- p; q[j] <- q[j] + 2 * sign(q[j])
  dl <- fit$obj$fn(p) - fit$obj$fn(q)
  cat(sprintf("RES seed=%d code=%d t=%.4g dl=%.3g edge=%s\n", s, fit$opt$convergence, p[j], dl,
              any(frmtmb:::se_edge_cor(fit, length(p)))))
}
