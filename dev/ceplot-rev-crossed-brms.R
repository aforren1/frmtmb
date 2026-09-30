# Reviewer check (lane ceplot): the renamed-level draws against brms
# 2.23.0 at the same draws, on shapes the worker did not try:
#   B  two unseen g:h combinations and one observed in one call
#   F  mm(g1, g2) beside (1 | h) + (1 | g1:h), g1:h unseen
#   G  (1 | g:h) and (1 | g:k) both unseen in the same rows
# brms is initialized at frmtmb's hand-built draws, algorithm =
# "fixed_param", one chain per draw (dev/ceplot-crossed-brms.R's
# construction), and the level order of each brms term is read from a
# probe fit rather than assumed. Seeds: shapes as in
# dev/ceplot-rev-shapes.R, draws 1, curves 1 to 3.
#   Rscript dev/ceplot-rev-crossed-brms.R > dev/ceplot-rev-log/crossed-brms.txt
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(brms)
  library(frmtmb)
  library(frmtmb.sample)
})
options(mc.cores = 1)
source("C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-rev-shapes.R")
cat("frmtmb from", find.package("frmtmb"), "\n")
say <- function(...) cat(sprintf(...), "\n", sep = "")
f4 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
N <- 200
S <- shape_fits()
bform <- list(
  B = y ~ x + (1 | g) + (1 | h) + (1 | g:h),
  F = y ~ x + (1 | mm(g1, g2)) + (1 | h) + (1 | g1:h),
  G = y ~ x + (1 | g) + (1 | h) + (1 | k) + (1 | g:h) + (1 | g:k))
qb <- function(expr) suppressMessages(suppressWarnings(expr))
for (nm in names(bform)) {
  fit <- S[[nm]]$fit
  dat <- fit$frame$data_frame
  ds <- hand(fit, N)
  M <- ds$draws
  sdat <- standata(bform[[nm]], data = dat)
  K <- sum(grepl("^N_[0-9]+$", names(sdat)))
  # probe: sd_k = k, z_k[j] = j / 1000, so r = k * j / 1000 names the
  # group of k and the level of j
  pinit <- list(c(list(b = array(0, 1), Intercept = 0, sigma = 1),
                  setNames(lapply(seq_len(K), function(k) array(k, 1)),
                           paste0("sd_", seq_len(K))),
                  setNames(lapply(seq_len(K), function(k) {
                    matrix(seq_len(sdat[[paste0("N_", k)]]) / 1000, 1)
                  }), paste0("z_", seq_len(K)))))
  pb <- qb(brm(bform[[nm]], data = dat, algorithm = "fixed_param",
               chains = 1, iter = 1, warmup = 0, init = pinit, refresh = 0,
               seed = 1, silent = 2))
  pm <- as_draws_matrix(pb)
  rn <- grep("^r_", colnames(pm), value = TRUE)
  rv <- as.numeric(pm[1, rn])
  kk <- round(rv * 1000) %/% vapply(round(rv * 1000), function(v) 1, 1)
  # r = k * j / 1000: recover k from the sd names instead, then j
  grp_of <- sub("^r_(.*)\\[.*$", "\\1", rn)
  sdn <- grep("^sd_", colnames(pm), value = TRUE)
  k_of_grp <- setNames(as.numeric(pm[1, sdn]), sub("^sd_(.*)__Intercept$",
                                                   "\\1", sdn))
  j_of <- round(rv * 1000 / k_of_grp[grp_of])
  # frmtmb: each block's sd column and its group's brms name
  th <- exp(M[, grep("^theta_", colnames(M))])
  gname <- vapply(fit$frame$re_blocks, function(bk) {
    if (!is.null(bk$components[[1]]$mm)) {
      paste0("mm", paste(bk$components[[1]]$mm$gvars, collapse = ""))
    } else bk$group_name
  }, "")
  tidx <- vapply(fit$frame$re_blocks, function(bk) as.integer(bk$theta_idx[1]), 1L)
  inits <- lapply(seq_len(N), function(i) {
    out <- list(b = array(M[i, "b_x"], 1),
                Intercept = M[i, "b_Intercept"] + mean(dat$x) * M[i, "b_x"],
                sigma = M[i, "sigma"])
    for (g in names(k_of_grp)) {
      k <- k_of_grp[[g]]
      s <- th[i, tidx[match(g, gname)]]
      labs <- rn[grp_of == g][order(j_of[grp_of == g])]
      out[[paste0("sd_", k)]] <- array(s, 1)
      out[[paste0("z_", k)]] <- matrix(M[i, labs] / s, 1)
    }
    out
  })
  b <- qb(brm(bform[[nm]], data = dat, algorithm = "fixed_param",
              chains = N, iter = 1, warmup = 0, init = inits, refresh = 0,
              seed = 1, silent = 2))
  bm <- as_draws_matrix(b)
  say("== %s: brms draws %d; max |r diff| %.3g", nm, ndraws(b),
      max(abs(bm[, rn] - M[, rn])))
  cond <- S[[nm]]$cond
  for (s in 1:3) {
    set.seed(s)
    a <- qb(conditional_effects(b, "x", resolution = 3, re_formula = NULL,
                                conditions = cond,
                                sample_new_levels = "gaussian"))$x
    d <- qb(conditional_effects(ds, "x", resolution = 3, re_formula = NULL,
                                conditions = cond, seed = s))$x
    for (cv in levels(factor(a$cond__))) {
      aa <- a[a$cond__ == cv, ]
      dd <- d[d$cond__ == cv, ]
      say("  seed %d cond %s: brms est %s lo %s | frmtmb est %s lo %s | max diff est %.3g lo %.3g up %.3g",
          s, cv, f4(aa$estimate__), f4(aa$lower__), f4(dd$estimate__),
          f4(dd$lower__), max(abs(aa$estimate__ - dd$estimate__)),
          max(abs(aa$lower__ - dd$lower__)),
          max(abs(aa$upper__ - dd$upper__)))
    }
  }
}
say("done")
