# F1's order through every route of the exact gp(): the identity with the centered
# objective, the kriging, by, gr, ties, two terms, and the names in the draws.
# Rows on the lane's build (seed 20261103); the tests pin the pass/fail.
#   Rscript dev/gpgap-order-check.R
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
gsn <- getFromNamespace("gp_sample_nugget", "frmtmb.sample")
set.seed(20261103)
n <- 60
d <- data.frame(x = round(runif(n, 0, 6), 3), z = round(runif(n, 0, 6), 3),
                f = factor(sample(c("a", "b"), n, TRUE)),
                t = sample(seq(0, 6, by = 0.25), n, TRUE))
d$y <- sin(2 * d$x) + 0.5 * cos(d$z) + 0.3 * rnorm(n)
d$yf <- ifelse(d$f == "a", sin(2 * d$x), cos(2 * d$x)) + 0.3 * rnorm(n)
d$yt <- sin(d$t) + 0.3 * rnorm(n)
specs <- list(
  plain = quote(frm(bf(y ~ gp(x)), data = d)),
  two_terms = quote(frm(bf(y ~ gp(x) + gp(z)), data = d)),
  aniso = quote(frm(bf(y ~ gp(x, z, iso = FALSE)), data = d)),
  by_f = quote(frm(bf(yf ~ gp(x, by = f)), data = d)),
  rows_tied = quote(frm(bf(yt ~ gp(t)), data = d)),
  gr_false_tied = quote(frm(bf(yt ~ gp(t, gr = FALSE)), data = d)),
  poisson = quote(frm(bf(k ~ gp(x)), family = poisson(),
                      data = transform(d, k = rpois(n, exp(1 + sin(2 * x)))))),
  sigma_gp = quote(frm(bf(y ~ 1, sigma ~ gp(x)), data = d))
)
for (nm in names(specs)) {
  fit <- suppressWarnings(eval(specs[[nm]]))
  ib <- which(vapply(fit$frame$re_blocks, function(b) {
    identical(b$covstruct, "gp")
  }, NA))
  fs <- gsn(fit, ib)
  f_c <- build_objective(local({
    fr <- fs$frame
    fr
  }))
  fr_n <- fs$frame
  fr_n$ncp_blocks <- ib
  f_n <- build_objective(fr_n)
  p <- fit$estimates
  pz <- p
  jac <- 0
  ordok <- TRUE
  for (i in ib) {
    bk <- fs$frame$re_blocks[[i]]
    th <- p$theta[bk$theta_idx]
    L <- frmtmb:::ncp_block_chol(bk, th)
    # b is a field the model can draw at this nugget: a typical z through L,
    # which moves off frm()'s own b only by what the smaller nugget allows
    z <- ncp_unscale_b(bk, ncp_scale_b(bk, rnorm(bk$dim), th), th)
    pz$b[bk$b_idx] <- z
    p$b[bk$b_idx] <- as.numeric(L %*% z)
    jac <- jac + as.numeric(determinant(L, logarithm = TRUE)$modulus)
    ordok <- ordok && identical(sort(bk$gp_order), seq_len(bk$dim))
  }
  cl <- f_c(p)
  nl <- f_n(pz)
  cat(sprintf(paste0("GPORD %-14s blocks %d positions %s | centered - noncentered ",
                     "- log|L| %.3e of %.1f | order a permutation %s\n"),
              nm, length(ib),
              paste(vapply(fs$frame$re_blocks[ib], `[[`, 0L, "dim"), collapse = ","),
              cl - nl - jac, abs(jac), ordok))
}

# draws: names and the field's order on an ordinary run from a formula. The field
# must be in the sorted order of the positions, as frm()'s is.
set.seed(20261103)
ds <- suppressMessages(suppressWarnings(
  frm_sample(bf(y ~ gp(x)), data = d, family = gaussian(), chains = 2, iter = 600,
             warmup = 300, refresh = 0, seed = 7)))
dr <- as.matrix(ds$draws)
bc <- grep("^b[[]", colnames(dr), value = TRUE)
fit <- suppressWarnings(frm(bf(y ~ gp(x)), data = d))
cat(sprintf("GPORD draws: %d b columns, first %s last %s; cor with frm()'s field %.4f\n",
            length(bc), bc[1], bc[length(bc)],
            cor(colMeans(dr[, bc]), as.numeric(fit$estimates$b))))
ep <- posterior_epred(ds, newdata = data.frame(x = c(0.111, 3.333, 5.555)))
cat(sprintf("GPORD kriging at new positions: means %s against sin(2 x) %s; sd %s\n",
            paste(sprintf("%.3f", colMeans(ep)), collapse = " "),
            paste(sprintf("%.3f", sin(2 * c(0.111, 3.333, 5.555)) +
                            mean(0.5 * cos(d$z))), collapse = " "),
            paste(sprintf("%.3f", apply(ep, 2, sd)), collapse = " ")))
at <- posterior_epred(ds, newdata = d[1:3, "x", drop = FALSE])
fr1 <- posterior_epred(ds)[, 1:3]
cat(sprintf("GPORD epred at observed rows, newdata against none: max diff %.2e\n",
            max(abs(at - fr1))))
