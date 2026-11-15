# Review of lane parity2, item 1497: sample_new_levels = "uncertainty".
# Part A: layout of the picked level on hand-built draws (correlated
#   slopes, gr(by =), two grouping factors), per draw.
# Part B: the boot band on a fit and the band on NUTS draws, against
#   brms (dev/reviews/parity2-review-snl-brms.R), on data with strong
#   shrinkage: 30 groups of 2 rows, sd_g 1, sigma 2.5. Seed 31.
#   Rscript dev/reviews/parity2-review-snl.R
.libPaths(c("C:/Users/adf44/source/r/wt-parity2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-parity2/dev/stan-cache")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    find.package("frmtmb.sample"), "\n")
part <- commandArgs(TRUE)[1]
if (is.na(part)) part <- "AB"

hand_draws <- function(fit, n = 40L, seed = 1L, sd = 0.05) {
  tpl <- fit$frame[["par_template"]]
  fx <- fit$frame[["betad_fixed_idx"]]
  est <- unlist(lapply(names(tpl), function(cp) {
    v <- fit$estimates[[cp]]
    if (cp == "betad" && length(fx)) v <- v[-fx]
    v
  }))
  lab <- frmtmb::brms_par_labels(fit)
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, sd),
              n, dimnames = list(NULL, lab))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  structure(list(stanfit = NULL, draws = M, fit = fit, data = fit$data),
            class = "frmtmb_draws")
}
spag <- function(ds, cond, xs = c(-1, 0, 1), ...) {
  ce <- conditional_effects(ds, "x", int_conditions = list(x = xs),
                            conditions = cond, spaghetti = TRUE, seed = 4,
                            ...)
  s <- attr(ce$x, "spaghetti")
  matrix(s$estimate__, ncol = length(xs), byrow = TRUE)
}

if (grepl("A", part)) {
  cat("\n#### A1 correlated slopes (1 + x | g)\n")
  set.seed(5)
  G <- 8
  d1 <- data.frame(x = rnorm(160), g = factor(rep(1:G, each = 20)))
  u <- MASS::mvrnorm(G, c(0, 0), matrix(c(1, 0.6, 0.6, 0.5), 2))
  d1$y <- 1 + 0.5 * d1$x + u[d1$g, 1] + u[d1$g, 2] * d1$x + rnorm(160, sd = 0.5)
  f1 <- frm(bf(y ~ x + (1 + x | g)), family = gaussian(), data = d1)
  bk <- f1$frame$re_blocks[[1]]
  bb <- f1$estimates$b[bk$b_idx]
  re <- ranef(f1)$g
  est <- if (length(dim(re)) == 3L) re[, "Estimate", ] else
    as.matrix(re)[, 1:2]
  cat("b level-major equals ranef:",
      isTRUE(all.equal(unname(matrix(bb, ncol = 2, byrow = TRUE)),
                       unname(est))), " ranef cols:",
      paste(colnames(as.matrix(re)), collapse = " "), "\n")
  ds <- hand_draws(f1)
  xs <- c(-1, 0, 1)
  pop <- spag(ds, list(), re_formula = NA)
  nu <- spag(ds, list(), re_formula = NULL)
  sh <- nu - pop
  ri <- ds$draws[, grep("^r_g\\[.*,Intercept\\]", colnames(ds$draws))]
  rs <- ds$draws[, grep("^r_g\\[.*,x\\]", colnames(ds$draws))]
  ok <- vapply(seq_len(nrow(sh)), function(i) {
    pred <- outer(ri[i, ], rep(1, 3)) + outer(rs[i, ], xs)
    any(apply(abs(sweep(pred, 2, sh[i, ])), 1, max) < 1e-10)
  }, NA)
  cat("draws whose shift is one level's (Intercept, x) pair:", sum(ok),
      "of", length(ok), "\n")
  mixed <- vapply(seq_len(nrow(sh)), function(i) {
    any(abs(ri[i, ] - sh[i, 2]) < 1e-10)
  }, NA)
  cat("draws whose shift at x = 0 is some level's intercept:", sum(mixed),
      "\n")
  a1 <- spag(ds, list(), re_formula = NULL)
  cat("same seed twice identical:", identical(a1, nu), "\n")

  cat("\n#### A2 gr(g, by = f): the pick is within the row's f level\n")
  set.seed(6)
  d2 <- data.frame(x = rnorm(240), g = factor(rep(1:12, each = 20)))
  d2$f <- factor(ifelse(as.integer(d2$g) <= 6, "a", "b"))
  ug <- ifelse(as.integer(levels(d2$g)) <= 6, rnorm(12, 0, 0.3),
               rnorm(12, 3, 2))
  d2$y <- 1 + 0.5 * d2$x + ug[d2$g] + rnorm(240, sd = 0.5)
  f2 <- frm(bf(y ~ x + f + (1 | gr(g, by = f))), family = gaussian(),
            data = d2)
  cat("blocks:", length(f2$frame$re_blocks), " b_idx lengths:",
      vapply(f2$frame$re_blocks, function(b) length(b$b_idx), 1L), "\n")
  ds2 <- hand_draws(f2)
  rg <- ds2$draws[, grep("^r_g\\[", colnames(ds2$draws))]
  cat("r_g columns:", paste(colnames(rg), collapse = " "), "\n")
  la <- paste0("r_g[", 1:6, ",Intercept]")
  lb <- paste0("r_g[", 7:12, ",Intercept]")
  for (fl in c("a", "b")) {
    pop <- spag(ds2, list(f = fl), re_formula = NA)
    nu <- spag(ds2, list(f = fl), re_formula = NULL)
    sh <- (nu - pop)[, 2]
    ina <- vapply(seq_along(sh), function(i) any(abs(rg[i, la] - sh[i]) < 1e-10), NA)
    inb <- vapply(seq_along(sh), function(i) any(abs(rg[i, lb] - sh[i]) < 1e-10), NA)
    cat("f =", fl, ": shift among f=a levels", sum(ina), " among f=b levels",
        sum(inb), " of", length(sh), "\n")
  }

  cat("\n#### A3 two grouping factors, both new\n")
  set.seed(7)
  d3 <- data.frame(x = rnorm(200), g = factor(sample(1:6, 200, TRUE)),
                   h = factor(sample(1:5, 200, TRUE)))
  d3$y <- 1 + 0.5 * d3$x + rnorm(6)[d3$g] + rnorm(5, 0, 0.7)[d3$h] +
    rnorm(200, sd = 0.5)
  f3 <- frm(bf(y ~ x + (1 | g) + (1 | h)), family = gaussian(), data = d3)
  ds3 <- hand_draws(f3)
  rg <- ds3$draws[, grep("^r_g\\[", colnames(ds3$draws))]
  rh <- ds3$draws[, grep("^r_h\\[", colnames(ds3$draws))]
  pop <- spag(ds3, list(), re_formula = NA)
  nu <- spag(ds3, list(), re_formula = NULL)
  sh <- (nu - pop)[, 2]
  ok <- vapply(seq_along(sh), function(i) {
    any(abs(outer(rg[i, ], rh[i, ], "+") - sh[i]) < 1e-10)
  }, NA)
  cat("shift = r_g[j] + r_h[k] on", sum(ok), "of", length(ok), "draws\n")
  nug <- spag(ds3, list(g = "1"), re_formula = NULL)
  sh1 <- (nug - pop)[, 2]
  ok1 <- vapply(seq_along(sh1), function(i) {
    any(abs(rg[i, "r_g[1,Intercept]"] + rh[i, ] - sh1[i]) < 1e-10)
  }, NA)
  cat("g = 1 set, h new: shift = r_g[1] + r_h[k] on", sum(ok1), "of",
      length(ok1), "draws\n")
}

if (grepl("B", part)) {
  cat("\n#### B shrinkage: 30 groups of 2 rows, sd_g 1, sigma 2.5\n")
  set.seed(31)
  G <- 30
  ds_ <- data.frame(x = rnorm(2 * G), g = factor(rep(1:G, each = 2)))
  ds_$y <- 1 + 0.5 * ds_$x + rnorm(G)[ds_$g] + rnorm(2 * G, sd = 2.5)
  saveRDS(ds_, "C:/Users/adf44/source/r/frmtmb-wt-parity2/dev/reviews/parity2-review-out/snl-data.rds")
  fb <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = ds_)
  print(VarCorr(fb))
  bk <- fb$frame$re_blocks[[1]]
  modes <- fb$estimates$b[bk$b_idx]
  cat("sd of the 30 modes:", sd(modes), "\n")
  ic <- list(x = c(-1, 0, 1))
  wid <- function(ce) ce$x$upper__ - ce$x$lower__
  for (snl in c("uncertainty", "gaussian")) {
    t0 <- proc.time()[["elapsed"]]
    ce <- suppressMessages(conditional_effects(
      fb, "x", int_conditions = ic, re_formula = NULL, band = "boot",
      boot = 200, seed = 2, sample_new_levels = snl))
    cat("fit boot 200 seed 2", snl, "width:", format(wid(ce), digits = 4),
        " (", round(proc.time()[["elapsed"]] - t0, 1), "s)\n")
  }
  cw <- conditional_effects(fb, "x", int_conditions = ic, re_formula = NULL)
  cat("fit wald width:", format(wid(cw), digits = 4), "\n")
  set.seed(1)
  smp <- suppressMessages(suppressWarnings(
    frm_sample(fb, chains = 2, iter = 2000, seed = 1, refresh = 0)))
  for (snl in c("uncertainty", "gaussian")) {
    ce <- conditional_effects(smp, "x", int_conditions = ic,
                              re_formula = NULL, seed = 2,
                              sample_new_levels = snl)
    cat("NUTS draws", snl, "width:", format(wid(ce), digits = 4), "\n")
  }
  r <- as.matrix(smp$draws)[, grep("^r_g\\[", colnames(smp$draws))]
  cat("sd over draws and levels of r_g draws:", sd(as.vector(r)), "\n")
  cat("sd_g posterior median:",
      median(smp$draws[, grep("^sd_g", colnames(smp$draws))[1]]), "\n")
}
