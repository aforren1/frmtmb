# Short check, item 5: sar()'s rho interval. y ~ x + sar(W) (lag),
# row-standardized rook W 5 x 5, true rho 0.9, 100 seeds; per-seed
# intervals saved, compared lane (clamped) against base (raw) in the base
# arm. Lane arm also: frm_multiple() pools then clamps.
source("dev/reviews/arspat2-short-common.R")
z <- qnorm(0.975)
rho <- 0.9
W <- rook(5)
W <- W / rowSums(W)
n <- nrow(W)
mkd <- function(seed) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  d$y <- as.vector(solve(diag(n) - rho * W, 1 + 0.5 * d$x + rnorm(n)))
  d
}
res <- do.call(rbind, lapply(seq_len(100), function(seed) {
  fit <- quiet(frm(y ~ x + sar(W), data = mkd(seed), data2 = list(W = W)))
  ac <- fit$frame$autocor[[1L]]
  cv <- quiet(confint_varcorr(fit))
  s <- quiet(summary(fit))$cor_pars
  data.frame(seed = seed, est = cv$estimate, lwr = cv$lwr, upr = cv$upr,
             type = cv$type, ee = s[1, "Est.Error"], s_lwr = s[1, 3],
             s_upr = s[1, 4], lo = ac$aux$lo, hi = ac$aux$hi,
             conv = fit$opt$convergence)
}))
saveRDS(res, file.path(outdir, sprintf("cover-%s.rds", arm)))
cat(sprintf("[%s] type %s; coverage %.3f; upr >= hi %d; lwr <= lo %d; NA %d; summary == confint %s\n",
            arm, paste(unique(res$type), collapse = ","),
            mean(res$lwr <= rho & rho <= res$upr), sum(res$upr >= res$hi),
            sum(res$lwr <= res$lo), sum(is.na(res$lwr)),
            isTRUE(all.equal(c(res$s_lwr, res$s_upr), c(res$lwr, res$upr)))))
if (arm == "base") {
  L <- readRDS(file.path(outdir, "cover-lane.rds"))
  want_l <- pmin(pmax(res$lwr, res$lo), res$hi)
  want_u <- pmin(pmax(res$upr, res$lo), res$hi)
  cat(sprintf("lane clamped vs clamp(base raw): max |diff| lwr %.2e upr %.2e; est %.2e; Est.Error %.2e\n",
              max(abs(L$lwr - want_l)), max(abs(L$upr - want_u)),
              max(abs(L$est - res$est)), max(abs(L$ee / res$ee - 1))))
  cat(sprintf("coverage lane %.3f base %.3f; seeds where base upr > hi %d (lane upr == hi there: %d)\n",
              mean(L$lwr <= rho & rho <= L$upr),
              mean(res$lwr <= rho & rho <= res$upr), sum(res$upr > res$hi),
              sum(L$upr[res$upr > res$hi] == L$hi[res$upr > res$hi])))
}
if (arm == "lane") {
  base <- mkd(3)
  set.seed(99)
  imps <- lapply(1:4, function(i) {
    d <- base
    d$y <- d$y + rnorm(n, 0, 0.1)
    d
  })
  mo <- quiet(frm_multiple(y ~ x + sar(W), data = imps,
                           data2 = list(W = W)))
  pv <- mo$pooled_varcorr
  print(pv)
  fits <- mo$fits
  ac <- fits[[1]]$frame$autocor[[1L]]
  est <- vapply(fits, function(f) quiet(confint_varcorr(f))$estimate, 0)
  se <- vapply(fits, function(f) quiet(summary(f))$cor_pars[1, "Est.Error"], 0)
  m <- length(fits)
  Tv <- mean(se^2) + (1 + 1 / m) * stats::var(est)
  tq <- qt(0.975, pv$df)
  raw <- mean(est) + c(-1, 1) * tq * sqrt(Tv)
  cat(sprintf("frm_multiple: type %s; pooled est %.10f vs mean %.10f; raw pooled interval (%.6f, %.6f), hi %.6f; reported (%.6f, %.6f); clamp(raw) diff %.2e\n",
              pv$type, pv$estimate, mean(est), raw[1], raw[2], ac$aux$hi,
              pv$lwr, pv$upr,
              max(abs(c(pv$lwr, pv$upr) - pmin(pmax(raw, ac$aux$lo), ac$aux$hi)))))
}
