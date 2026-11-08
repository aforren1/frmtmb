# Lane latperf: a cumulative latent fit whose parameters keep their SEs.
.libPaths(c("C:/Users/adf44/source/r/wt-latperf-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
for (seed in 15:20) for (mult in c(2, 3)) {
  set.seed(seed)
  G <- 30L; K <- 12L
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  d <- d[order(d$g, d$t), ]
  n <- nrow(d)
  d$x <- stats::rnorm(n)
  e <- as.vector(replicate(G, stats::arima.sim(list(ar = 0.6), K, sd = 0.6)))
  d$ord <- cut(mult * (0.4 * d$x + e) + stats::rlogis(n),
               c(-Inf, -1, 0, 1, Inf), ordered_result = TRUE)
  w <- character(0)
  fit <- withCallingHandlers(frm(ord ~ x + ar(t, g, cov = TRUE),
                                 family = cumulative(), data = d),
    warning = function(x) {
      w <<- c(w, substr(conditionMessage(x), 1, 40))
      invokeRestart("muffleWarning")
    })
  se <- sqrt(diag(vcov(fit, full = TRUE)))
  cat(sprintf("seed %d mult %g finite %s warns %s\n", seed, mult,
              all(is.finite(se)), paste(w, collapse = " | ")))
}
