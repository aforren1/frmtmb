# frmtmb's prediction of an exact gp(by = ) at NEW positions against
# brms's, on brms's own stored fit6 (volume ~ Trt + gp(Age, by = Trt)),
# draw by draw: each brms draw is mapped into frmtmb's parameter space
# (theta = log sdgp, log(lscale * dmax); b = the field at the fitted
# positions, L zgp with brms's own kernel), and frmtmb's conditional law
# of the field at the new positions is compared with brms's at the same
# draw. Then both sample it: brms's posterior_epred() against frmtmb's
# per-draw predictor with the kriging draw on (krige_draw), the route
# frmtmb.sample's draws methods take.
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(brms); library(frmtmb)})
cat("lib:", find.package("frmtmb"), "| brms", format(packageVersion("brms")),
    "\n")
fit6 <- brms:::rename_pars(brms:::brmsfit_example6)
dat <- fit6$data
fr <- frm(bf(volume ~ Trt + gp(Age, by = Trt, gr = TRUE)) + gaussian() +
            bf(count ~ Trt + Age) + poisson() + set_rescor(FALSE),
          data = dat)
dm <- as.matrix(fit6)
nd <- data.frame(Age = c(-1.9, 0.33, 1.1, 3.1, 4),
                 Trt = factor(c("0", "0", "1", "1", "0"),
                              levels = levels(dat$Trt)),
                 count = 0)
gis <- fr$frame$linpreds[["volume.mu"]]$gps
bks <- lapply(gis, function(g) fr$frame$re_blocks[[g$block_id]])
sdat <- standata(fit6, internal = TRUE)
lev <- levels(dat$Trt)
ids <- seq_len(nrow(dm))
# brms's conditional law at one draw, its own formulas
# (brms:::.predictor_gp_new: nug 1e-12 on the old points, 1e-8 on new)
brms_cond <- function(s, j, xnew) {
  sd <- dm[s, paste0("sdgp_volume_gpAgeTrt", lev[j])]
  ls <- dm[s, paste0("lscale_volume_gpAgeTrt", lev[j])]
  dmax <- sdat[[paste0("dmax_volume_1_", j)]]
  x <- as.numeric(sdat[[paste0("Xgp_volume_1_", j)]])
  z <- dm[s, grep(paste0("^zgp_volume_gpAgeTrt", lev[j], "\\["),
                  colnames(dm))]
  k <- function(a, b) sd^2 * exp(-outer(a, b, "-")^2 / (2 * ls^2))
  K <- k(x, x) + diag(1e-12, length(x))
  L <- t(chol(K))
  f <- drop(L %*% z)
  xs <- xnew / dmax
  # brms's own arithmetic, L^-1 and v'v, which stays positive where
  # k** - k' K^-1 k through solve() does not
  Li <- solve(L)
  v <- Li %*% t(k(xs, x))
  list(f = f, x = x * dmax, sd = sd, ls = ls * dmax,
       mean = drop(t(v) %*% (Li %*% f)),
       var = diag(k(xs, xs)) - colSums(v^2) + 1e-8)
}
# frmtmb's conditional law at the same draw, through frm_lp_basis on a
# fit carrying that draw's parameters
frm_at <- function(s) {
  g <- fr
  est <- g$estimates
  for (j in seq_along(bks)) {
    bc <- brms_cond(s, j, 0)
    bk <- bks[[j]]
    est$theta[bk$theta_idx] <- c(log(bc$sd), log(bc$ls))
    pos <- gis[[j]]$positions[, 1]
    est$b[bk$b_idx] <- bc$f[match(pos, round(bc$x, 12))]
    if (anyNA(est$b[bk$b_idx])) {
      est$b[bk$b_idx] <- bc$f[sapply(pos, function(p) which.min(abs(bc$x - p)))]
    }
  }
  est$beta[] <- 0
  g$estimates <- est
  g$cache <- new.env(parent = emptyenv())
  g
}
res <- NULL
for (s in ids) {
  g <- frm_at(s)
  lb <- frm_lp_basis(g, newdata = nd, resp = "volume", re_formula = NA)
  for (r in seq_len(nrow(nd))) {
    j <- match(as.character(nd$Trt[r]), lev)
    bc <- brms_cond(s, j, nd$Age[r])
    res <- rbind(res, data.frame(
      draw = s, row = r, brms_mean = bc$mean, frm_mean = lb$eta[r],
      brms_sd = sqrt(bc$var), frm_sd = sqrt(lb$extra_var[r])))
  }
}
res$dmean <- (res$frm_mean - res$brms_mean) / res$brms_sd
res$rsd <- res$frm_sd / res$brms_sd
cat("== per draw, every draw x 5 new rows: frmtmb's conditional law of the",
    "field against brms's\n")
for (r in seq_len(nrow(nd))) {
  q <- res[res$row == r, ]
  qa <- stats::quantile(abs(q$dmean), c(0.5, 0.9, 1))
  qr <- stats::quantile(q$rsd, c(0.05, 0.5, 0.95))
  cat(sprintf(paste0("row %d Age %5.2f Trt %s | |mean gap| / brms sd: ",
                     "median %.1e q90 %.1e max %.1e | sd ratio q05 %.4f ",
                     "median %.6f q95 %.4f\n"),
              r, nd$Age[r], nd$Trt[r], qa[1], qa[2], qa[3], qr[1], qr[2],
              qr[3]))
}
# sampled: brms's posterior_epred against frmtmb's kriging draws at the
# SAME parameter draws, the field part only (beta held at 0 on our side,
# so brms's fixed part is subtracted from its draws)
set.seed(11)
ep <- posterior_epred(fit6, newdata = nd, resp = "volume", draw_ids = ids)
fx <- dm[ids, "b_volume_Intercept"] +
  outer(dm[ids, "b_volume_Trt1"], as.numeric(nd$Trt == "1"))
ep_field <- ep - fx
set.seed(12)
fe <- t(vapply(ids, function(s) {
  g <- frm_at(s)
  g[["krige_draw"]] <- TRUE
  as.numeric(frm_linpred(g, newdata = nd, resp = "volume", re_formula = NA))
}, numeric(nrow(nd))))
cat("== sampled at the same draws: sd across draws of the field at",
    "each new row\n")
for (r in seq_len(nrow(nd))) {
  cat(sprintf(paste0("row %d Age %5.2f | brms %.5f | frmtmb %.5f | ratio ",
                     "%.4f | means %+.5f %+.5f\n"),
              r, nd$Age[r], stats::sd(ep_field[, r]), stats::sd(fe[, r]),
              stats::sd(fe[, r]) / stats::sd(ep_field[, r]),
              mean(ep_field[, r]), mean(fe[, r])))
}
# the same draws with the kriging draw OFF: the conditional mean only,
# which is what frmtmb.sample returned before
fm <- t(vapply(ids, function(s) {
  as.numeric(frm_linpred(frm_at(s), newdata = nd, resp = "volume",
                         re_formula = NA))
}, numeric(nrow(nd))))
cat("== kriging draw off (the 0.67.0 route): sd ratio to brms per row:",
    sprintf("%.4f", apply(fm, 2, stats::sd) / apply(ep_field, 2, stats::sd)),
    "\n")
cat("== share of draws whose two conditional laws agree to 1 percent in",
    "mean (in brms sd) and in sd, per row:",
    vapply(seq_len(nrow(nd)), function(r) {
      q <- res[res$row == r, ]
      sprintf("%.3f", mean(abs(q$dmean) < 0.01 & abs(q$rsd - 1) < 0.01))
    }, ""), "\n")
cat("== draws:", length(ids), "\n")
