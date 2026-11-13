# Review of lane names: draws names vs values, every accessor's order.
# Usage: Rscript dev/reviews/names-review-draws.R [case ...]
# Output: stdout (log), dev/reviews/names-review-draws.rds
a <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-names-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
options(mc.cores = 1)
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
set.seed(20261111)
n <- 72
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(rep(sprintf("g%d", 1:6), each = 12)),
                h = factor(rep(1:4, length.out = n)),
                t = rep(1:12, 6))
ug <- rnorm(6, 0, 0.6)[d$g]
d$y <- 1 + 0.5 * d$x + ug + rnorm(n, 0, exp(0.2 * d$z))
d$y2 <- 0.5 * d$y + rnorm(n)
d$ys <- sin(d$x) + rnorm(n, 0, 0.3)
d$k <- sample(1:4, n, TRUE)
d$ymo <- 1 + 0.4 * d$k + ug + rnorm(n)
d$xg <- rnorm(6)[d$g]
d$wg <- rnorm(6)[d$g]
d$sg <- 0.25
d$sz <- 0.3
d$yme <- 1 + 0.5 * d$xg + 0.3 * d$z - 0.2 * d$wg + rnorm(n)
d$cnt <- rpois(n, exp(0.5 + 0.3 * d$x))
d$yc <- factor(sample(c("a", "b", "c"), n, TRUE))
d$ymix <- ifelse(runif(n) < 0.4, rnorm(n, -1), rnorm(n, 2))
cases <- list(
  mv = list(bf(y ~ x + (1 | g)) + bf(y2 ~ z + (1 | g)) + set_rescor(FALSE),
            gaussian()),
  skew = list(bf(y ~ x + (1 | h), alpha ~ (1 | g), sigma ~ z + (1 | g)),
              skew_normal()),
  gp = list(bf(ys ~ gp(x)), gaussian()),
  smooth = list(bf(ys ~ s(x) + (1 | g)), gaussian()),
  mo = list(bf(ymo ~ mo(k) + (1 | g)), gaussian()),
  me = list(bf(yme ~ me(xg, sg, gr = g) + me(z, sz) + me(wg, sg, gr = g)),
            gaussian()),
  ar = list(bf(cnt ~ x + ar(t, g)), poisson()),
  cat = list(bf(yc ~ x + (1 | g)), categorical()),
  mix = list(bf(ymix ~ 1, mu1 ~ x, mu2 ~ (1 | g)),
             mixture(gaussian(), gaussian())),
  nl = list(bf(y ~ a + b * x, a ~ 1 + z + (1 | g), b ~ 1, nl = TRUE),
            gaussian())
)
if (length(a)) cases <- cases[a]
out <- list()
chk <- function(case, label, ok) {
  cat(sprintf("%-8s %-58s %s\n", case, label, if (isTRUE(ok)) "ok" else
    paste("FAIL", paste(format(ok), collapse = " "))))
  out[[case]][[label]] <<- ok
}
# brms's r_ name to (group, level, coef as ranef() names it)
r_parse <- function(nm) {
  m <- regmatches(nm, regexec("^r_([^[]+)\\[([^,]+),(.+)\\]$", nm))[[1]]
  grp <- m[2]
  pre <- ""
  if (grepl("__", grp)) {
    pre <- sub("^.*?__", "", grp)
    grp <- sub("__.*$", "", grp)
  }
  list(group = grp, level = m[3],
       coef = if (nzchar(pre)) paste0(pre, "_", m[4]) else m[4])
}
for (cs in names(cases)) {
  t0 <- proc.time()[["elapsed"]]
  ds <- tryCatch(suppressWarnings(suppressMessages(frm_sample(
    cases[[cs]][[1]], family = cases[[cs]][[2]], data = d, chains = 2,
    iter = 100, warmup = 50, refresh = 0, seed = 5, save_warmup = TRUE))),
    error = function(e) e)
  if (inherits(ds, "error")) {
    cat(cs, "SAMPLE ERROR:", conditionMessage(ds), "\n")
    next
  }
  cat(sprintf("\n== %s (%.0f s)\n", cs, proc.time()[["elapsed"]] - t0))
  v <- variables(ds)
  cat("variables:", head(v, 40), if (length(v) > 40) "...", "\n")
  cat("x$draws  :", head(colnames(ds$draws), 40), "\n")
  chk(cs, "setequal(colnames(x$draws), variables)",
      setequal(colnames(ds$draws), v) && !anyDuplicated(v))
  adf <- posterior::as_draws_df(ds)
  chk(cs, "as_draws_df names == variables",
      identical(setdiff(names(adf), c(".chain", ".iteration", ".draw")), v))
  m <- as.matrix(ds)
  chk(cs, "as.matrix colnames == variables", identical(colnames(m), v))
  bad <- v[!vapply(v, function(nm) identical(unname(m[, nm]),
                                            unname(ds$draws[, nm])), NA)]
  chk(cs, "as.matrix values by name == x$draws", if (!length(bad)) TRUE
      else head(bad))
  arr <- as.array(ds)
  chk(cs, "as.array dimnames[[3]] == variables",
      identical(dimnames(arr)[[3]], v))
  am <- unclass(posterior::as_draws_matrix(posterior::as_draws_array(ds)))
  bad <- v[!vapply(v, function(nm) identical(unname(am[, nm]),
                                            unname(ds$draws[, nm])), NA)]
  chk(cs, "as_draws_array values by name", if (!length(bad)) TRUE
      else head(bad))
  ad <- as.data.frame(ds)
  chk(cs, "as.data.frame names == variables", identical(names(ad), v))
  sm <- posterior::summarise_draws(adf)
  chk(cs, "summarise_draws order", identical(sm$variable, v))
  s <- tryCatch(summary(ds), error = function(e) e)
  if (!inherits(s, "error")) {
    rn <- if (is.matrix(s)) rownames(s) else NULL
    if (!is.null(rn)) {
      chk(cs, "summary rows in variables order",
          identical(rn, v[v %in% rn]))
    }
  }
  mc <- tryCatch(coda::as.mcmc(ds, combine_chains = TRUE),
                 error = function(e) e)
  if (!inherits(mc, "error")) {
    chk(cs, "as.mcmc(combine) colnames == variables",
        identical(colnames(mc), v))
    bad <- v[!vapply(v, function(nm) identical(unname(as.vector(mc[, nm])),
                                              unname(ds$draws[, nm])), NA)]
    chk(cs, "as.mcmc(combine) values by name", if (!length(bad)) TRUE
        else head(bad))
  } else cat("as.mcmc error:", conditionMessage(mc), "\n")
  # regex selection keeps variables() order (brms extract_pars())
  pm <- suppressWarnings(as.matrix(ds, pars = "^(b|sd|r)_"))
  chk(cs, "as.matrix(pars = regex) in variables order",
      identical(colnames(pm), v[grepl("^(b|sd|r)_", v)]))
  if (!identical(colnames(pm), v[grepl("^(b|sd|r)_", v)])) {
    cat("  pars= order:", head(colnames(pm), 8), "\n")
  }
  aw <- tryCatch(posterior::as_draws_array(ds, inc_warmup = TRUE),
                 error = function(e) e)
  if (inherits(aw, "error")) {
    chk(cs, "inc_warmup", conditionMessage(aw))
  } else {
    chk(cs, "inc_warmup variables order",
        identical(posterior::variables(aw), v))
    # post-warmup tail equals the stored draws, by name
    nw <- posterior::niterations(aw) - dim(arr)[1]
    nm1 <- v[1]
    ch1 <- unname(as.vector(aw[(nw + 1):posterior::niterations(aw), 1, nm1]))
    chk(cs, paste("inc_warmup tail by name", nm1),
        isTRUE(all.equal(ch1, unname(arr[, 1, nm1]))))
  }
  # ranef(summary = FALSE) of the draws against the r_ columns by name
  rn_cols <- v[startsWith(v, "r_")]
  if (length(rn_cols)) {
    re <- tryCatch(ranef(ds, summary = FALSE), error = function(e) e)
    if (inherits(re, "error")) {
      chk(cs, "draws ranef", conditionMessage(re))
    } else {
      bad <- character(0)
      for (nm in rn_cols) {
        p <- r_parse(nm)
        A <- re[[p$group]]
        val <- tryCatch(A[, p$level, p$coef], error = function(e) NULL)
        if (is.null(val) || !isTRUE(all.equal(unname(val),
                                              unname(ds$draws[, nm])))) {
          bad <- c(bad, nm)
        }
      }
      chk(cs, "draws ranef(summary=FALSE) == r_ columns by name",
          if (!length(bad)) TRUE else head(bad))
      cat("draws ranef names:", names(re), "; coefs:",
          vapply(re, function(A) paste(dimnames(A)[[3]], collapse = ","), ""),
          "\n")
    }
    vc <- tryCatch(VarCorr(ds), error = function(e) e)
    sdc <- v[startsWith(v, "sd_")]
    if (!inherits(vc, "error")) {
      bad <- character(0)
      for (nm in sdc) {
        gp <- sub("^sd_(.*)__([^_].*)$", "\\1", nm)
        cf <- sub("^sd_(.*)__([^_].*)$", "\\2", nm)
        if (grepl("__", gp)) {
          cf <- paste0(sub("^.*?__", "", gp), "_", cf)
          gp <- sub("__.*$", "", gp)
        }
        est <- tryCatch(vc[[gp]]$sd[cf, "Estimate"], error = function(e) NA)
        if (!isTRUE(all.equal(unname(est), mean(ds$draws[, nm])))) {
          bad <- c(bad, nm)
        }
      }
      chk(cs, "VarCorr sd Estimate == mean(sd_ column)",
          if (!length(bad)) TRUE else head(bad))
    }
    co <- tryCatch(coef(ds), error = function(e) e)
    if (!inherits(co, "error")) {
      cat("draws coef names:", names(co), "; coefs:",
          vapply(co, function(A) paste(dimnames(A)[[3]], collapse = ","), ""),
          "\n")
    }
  }
  bx <- v[startsWith(v, "b_") & !grepl("Intercept", v)][1]
  if (!is.na(bx)) {
    hy <- tryCatch(hypothesis(ds, paste(bx, "> 0"), class = NULL),
                   error = function(e) e)
    if (!inherits(hy, "error")) {
      chk(cs, paste("hypothesis", bx, "Estimate"),
          isTRUE(all.equal(hy$hypothesis$Estimate, mean(ds$draws[, bx]))))
    } else chk(cs, "hypothesis", conditionMessage(hy))
  }
  if (requireNamespace("bayesplot", quietly = TRUE)) {
    p <- tryCatch(bayesplot::mcmc_intervals_data(as.array(ds),
                                                  pars = head(v, 3)),
                  error = function(e) e)
    if (!inherits(p, "error")) {
      chk(cs, "bayesplot intervals m == mean? (median)",
          isTRUE(all.equal(p$m, unname(apply(m[, head(v, 3), drop = FALSE],
                                             2, stats::median)))))
    } else chk(cs, "bayesplot", conditionMessage(p))
  }
  lo <- tryCatch(suppressWarnings(loo::loo(ds)), error = function(e) e)
  chk(cs, "loo() runs", if (inherits(lo, "error")) conditionMessage(lo) else
    TRUE)
  ll <- tryCatch(log_lik(ds), error = function(e) e)
  if (!inherits(ll, "error")) cat("log_lik dimnames:",
                                  format(dimnames(ll)), "\n")
  # ML fit's variables() vs draws' variables(): same relative order
  vf <- variables(ds$fit)
  com <- intersect(vf, v)
  chk(cs, "ML variables() order == draws order on common names",
      identical(com, v[v %in% com]))
  cat("ML-only:", head(setdiff(vf, v), 10), " draws-only (non-r_):",
      head(setdiff(v[!grepl("^(r|s|zgp|Xme|zerr)_|^(zerr|lp__|lprior)",
                                v)], vf), 10), "\n")
  rm(ds)
  gc()
}
saveRDS(out, "dev/reviews/names-review-draws.rds")
cat("\nDONE\n")
