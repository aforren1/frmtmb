# Item 857: does the hmm-stationary test's assertion hold on every
# rounding? The response is perturbed by a few ulps, y * (1 + 2e-15 k),
# k = 1..20, and the test's three checks are recorded for each draw.
.libPaths(c("C:/Users/adf44/source/r/wt-silent-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.latent)})
cat("latent from", find.package("frmtmb.latent"), "\n")
sn_hmm_data <- function(seed = 1L, nseq = 20L, Tl = 40L, alpha = 5) {
  set.seed(seed)
  rsn <- function(n, mu, sigma, a) {
    d <- a / sqrt(1 + a^2)
    mu + sigma * (d * abs(stats::rnorm(n)) +
                    sqrt(1 - d^2) * stats::rnorm(n))
  }
  G <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE)
  do.call(rbind, lapply(seq_len(nseq), function(id) {
    s <- integer(Tl)
    s[1L] <- sample.int(2L, 1L)
    for (t in seq_len(Tl)[-1L]) {
      s[t] <- sample.int(2L, 1L, prob = G[s[t - 1L], ])
    }
    y <- ifelse(s == 1L, rsn(Tl, -3, 1, alpha), rsn(Tl, 3, 1, -alpha))
    data.frame(id = id, t = seq_len(Tl), y = y)
  }))
}
fam <- function() hmm(K = 2, skew_normal(), time = t, group = id)
d0 <- sn_hmm_data()
a1 <- "alpha1_(Intercept)"
rows <- NULL
for (k in 0:20) {
  d <- d0
  d$y <- d0$y * (1 + 2e-15 * k)
  fo <- bf(y ~ 1)
  u <- frm(fo, family = fam(), data = d, dry_run = "objective")
  bd <- u$estimates$betad
  stuck <- bd
  stuck[[a1]] <- 0
  fit <- suppressWarnings(frm(fo, family = fam(), data = d,
                              start = list(betad = stuck)))
  esc <- fit$opt[["stationary_escape"]]
  ll <- as.numeric(logLik(fit))
  ref <- max(vapply(c(2, -2), function(a) {
    s <- bd
    s[[a1]] <- a
    as.numeric(logLik(suppressWarnings(
      frm(fo, family = fam(), data = d, start = list(betad = s)))))
  }, 0))
  rows <- rbind(rows, data.frame(
    k = k, starts = if (is.null(esc)) NA else esc[["starts"]],
    gain_ratio = if (is.null(esc)) NA else esc[["gain"]] / abs(ll),
    alpha1 = fit$estimates$betad[[a1]],
    ref_minus_ll_ratio = (ref - ll) / abs(ref)))
}
print(format(rows, digits = 4), row.names = FALSE)
ok <- !is.na(rows$starts) & rows$starts == 2 & rows$gain_ratio > 1e-3 &
  abs(rows$alpha1) > 1 & rows$ref_minus_ll_ratio < 1e-8
cat(sprintf("SUMMARY all four checks hold on %d of %d draws\n",
            sum(ok), nrow(rows)))
