# The comparison of dev/gpby-brms-epred.R on a brms fit with enough
# draws for the across-draw spread to mean something: brms's own
# y ~ gp(x, by = f) fitted to test-gp-by.R's data, 2 chains of 1500
# (750 warmup), seed 2026. Draw by draw, frmtmb's conditional law of the
# field at new positions against brms's, and the spread of the sampled
# field: brms's posterior_epred() against frmtmb's per-draw predictor
# with the kriging draw on and off.
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(brms); library(frmtmb)})
cat("lib:", find.package("frmtmb"), "| brms", format(packageVersion("brms")),
    "\n")
set.seed(5)
n <- 60
d <- data.frame(x = round(stats::runif(n, 0, 6), 1),
                f = factor(rep(c("a", "b", "c"), length.out = n)),
                w = stats::runif(n, 0.5, 2))
d$y <- 0.5 + ifelse(d$f == "a", sin(d$x),
                    ifelse(d$f == "b", cos(d$x), 0.2 * d$x)) +
  stats::rnorm(n, 0, 0.3)
rds <- "dev/gpby-brms-epred2-fit.rds"
bfit <- if (file.exists(rds)) readRDS(rds) else {
  b <- brm(y ~ gp(x, by = f), data = d, chains = 2, iter = 1500,
           warmup = 750, seed = 2026, refresh = 0, backend = "rstan")
  saveRDS(b, rds)
  b
}
dm <- as.matrix(bfit)
fr <- frm(frmtmb::bf(y ~ gp(x, by = f)), data = d)
gis <- fr$frame$linpreds[["y.mu"]]$gps
bks <- lapply(gis, function(g) fr$frame$re_blocks[[g$block_id]])
sdat <- standata(bfit, internal = TRUE)
lev <- levels(d$f)
nd <- data.frame(x = c(2.55, 6.5, 7.5, 2.55, 6.5, 7.5),
                 f = factor(c("a", "a", "a", "c", "c", "c"), levels = lev))
brms_cond <- function(s, j, xnew) {
  sd <- dm[s, paste0("sdgp_gpxf", lev[j])]
  ls <- dm[s, paste0("lscale_gpxf", lev[j])]
  dmax <- sdat[[paste0("dmax_1_", j)]]
  x <- as.numeric(sdat[[paste0("Xgp_1_", j)]])
  z <- dm[s, grep(paste0("^zgp_gpxf", lev[j], "\\["), colnames(dm))]
  k <- function(a, b) sd^2 * exp(-outer(a, b, "-")^2 / (2 * ls^2))
  L <- t(chol(k(x, x) + diag(1e-12, length(x))))
  f <- drop(L %*% z)
  Li <- solve(L)
  v <- Li %*% t(k(xnew / dmax, x))
  list(f = f, x = x * dmax, sd = sd, ls = ls * dmax,
       mean = drop(t(v) %*% (Li %*% f)),
       var = diag(k(xnew / dmax, xnew / dmax)) - colSums(v^2) + 1e-8)
}
frm_at <- function(s) {
  g <- fr
  est <- g$estimates
  for (j in seq_along(bks)) {
    bc <- brms_cond(s, j, 0)
    bk <- bks[[j]]
    est$theta[bk$theta_idx] <- c(log(bc$sd), log(bc$ls))
    pos <- gis[[j]]$positions[, 1]
    est$b[bk$b_idx] <- bc$f[vapply(pos, function(p) {
      which.min(abs(bc$x - p))
    }, 1L)]
  }
  est$beta[] <- 0
  g$estimates <- est
  g$cache <- new.env(parent = emptyenv())
  g
}
ids <- seq_len(nrow(dm))
res <- NULL
for (s in ids) {
  lb <- frm_lp_basis(frm_at(s), newdata = nd, re_formula = NA)
  for (r in seq_len(nrow(nd))) {
    bc <- brms_cond(s, match(as.character(nd$f[r]), lev), nd$x[r])
    res <- rbind(res, data.frame(row = r, dmean = (lb$eta[r] - bc$mean) /
                                   sqrt(bc$var),
                                 rsd = sqrt(lb$extra_var[r] / bc$var)))
  }
}
cat("== per draw (", length(ids), "draws): frmtmb's conditional law of",
    "the field at a new position against brms's\n")
for (r in seq_len(nrow(nd))) {
  q <- res[res$row == r, ]
  qa <- stats::quantile(abs(q$dmean), c(0.5, 0.9))
  qr <- stats::quantile(q$rsd, c(0.05, 0.5, 0.95))
  cat(sprintf(paste0("row %d x %4.2f f %s | |mean gap| / brms sd median ",
                     "%.1e q90 %.1e | sd ratio q05 %.4f median %.4f q95 ",
                     "%.4f | agree to 1%% on %.3f\n"),
              r, nd$x[r], nd$f[r], qa[1], qa[2], qr[1], qr[2], qr[3],
              mean(abs(q$dmean) < 0.01 & abs(q$rsd - 1) < 0.01)))
}
set.seed(11)
ep <- posterior_epred(bfit, newdata = nd) - dm[, "b_Intercept"]
set.seed(12)
on <- t(vapply(ids, function(s) {
  g <- frm_at(s)
  g[["krige_draw"]] <- TRUE
  as.numeric(frm_linpred(g, newdata = nd, re_formula = NA))
}, numeric(nrow(nd))))
off <- t(vapply(ids, function(s) {
  as.numeric(frm_linpred(frm_at(s), newdata = nd, re_formula = NA))
}, numeric(nrow(nd))))
cat("== sampled at the same", length(ids), "draws: sd across draws of the",
    "field (brms | frmtmb draw on | ratio | frmtmb draw off | ratio)\n")
for (r in seq_len(nrow(nd))) {
  cat(sprintf("row %d x %4.2f f %s | %.5f | %.5f | %.4f | %.5f | %.4f\n",
              r, nd$x[r], nd$f[r], stats::sd(ep[, r]), stats::sd(on[, r]),
              stats::sd(on[, r]) / stats::sd(ep[, r]), stats::sd(off[, r]),
              stats::sd(off[, r]) / stats::sd(ep[, r])))
}
# the per-draw gaps in units of the field's POSTERIOR spread at the row,
# which is what a prediction interval is made of. A conditional sd near
# brms's own 1e-8 nugget makes a ratio of two conditional sds large and
# meaningless, so the conditional sds are reported as numbers too.
res2 <- NULL
for (s in ids) {
  lb <- frm_lp_basis(frm_at(s), newdata = nd, re_formula = NA)
  for (r in seq_len(nrow(nd))) {
    bc <- brms_cond(s, match(as.character(nd$f[r]), lev), nd$x[r])
    res2 <- rbind(res2, data.frame(row = r, gap = lb$eta[r] - bc$mean,
                                   bsd = sqrt(bc$var),
                                   fsd = sqrt(lb$extra_var[r])))
  }
}
cat("== per draw, in units of brms's posterior sd of the field at the row\n")
for (r in seq_len(nrow(nd))) {
  q <- res2[res2$row == r, ]
  ps <- stats::sd(ep[, r])
  cat(sprintf(paste0("row %d x %4.2f f %s | |mean gap| median %.2e q90 ",
                     "%.2e | conditional sd brms median %.2e, frmtmb ",
                     "median %.2e\n"),
              r, nd$x[r], nd$f[r], stats::median(abs(q$gap)) / ps,
              stats::quantile(abs(q$gap), 0.9) / ps,
              stats::median(q$bsd) / ps, stats::median(q$fsd) / ps))
}
