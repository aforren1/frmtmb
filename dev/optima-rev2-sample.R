# Reviewer of lane optima, re-check (a): frm_sample() on the lane's
# punch-round build against brms 2.23.0, the review's weak case and
# brms_monotonic fit1, 4 chains of 2000 (1000 warmup) on each side,
# sampler seed 7 (the lane used 1).
#   Rscript dev/optima-rev2-sample.R
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb", find.package("frmtmb"), "| frmtmb.sample",
    find.package("frmtmb.sample"), "\n")
out <- "C:/Users/adf44/source/r/frmtmb-wt-optima/dev/optima-rev2-out"
set.seed(1)
n <- 100
x <- sample(0:3, n, TRUE)
weak <- data.frame(x = x, y = 0.05 * x + rnorm(n))
set.seed(1234)
lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
strong <- data.frame(income, ls)
cases <- list(weak = list(f = y ~ mo(x), d = weak, nm = "simo_mox1"),
              strong = list(f = ls ~ mo(income), d = strong,
                            nm = "simo_moincome1"))
smry <- function(arr, cols) {
  t(vapply(cols, function(cn) {
    v <- arr[, , cn]
    c(mean = mean(v), sd = stats::sd(c(v)),
      mcse_mean = posterior::mcse_mean(v), mcse_sd = posterior::mcse_sd(v),
      rhat = posterior::rhat(v), ess = posterior::ess_bulk(v))
  }, numeric(6)))
}
res <- list()
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fit <- frm(bf(cs$f), data = cs$d, family = gaussian())
  msg <- character()
  s <- withCallingHandlers(
    frm_sample(fit, chains = 4, iter = 2000, warmup = 1000, seed = 7,
               cores = 1, refresh = 0),
    message = function(m) {
      msg <<- c(msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      msg <<- c(msg, paste("WARN", conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
  cat("\n==", nm, "frmtmb messages:\n")
  for (m in msg) cat("  ", substr(gsub("\n", " ", m), 1, 220), "\n")
  a <- posterior::as_draws_array(s)
  cat("variables:", paste(posterior::variables(a), collapse = " "), "\n")
  cols <- grep("^simo_", posterior::variables(a), value = TRUE)
  sp <- rstan::get_sampler_params(s$stanfit, inc_warmup = FALSE)
  divf <- sum(sapply(sp, function(z) sum(z[, "divergent__"])))
  sf <- smry(a, cols)
  b <- brms::brm(cs$f, data = cs$d, chains = 4, iter = 2000, seed = 7,
                 cores = 1, refresh = 0)
  ab <- posterior::as_draws_array(b)
  colsb <- grep("^simo_", posterior::variables(ab), value = TRUE)
  spb <- rstan::get_sampler_params(b$fit, inc_warmup = FALSE)
  divb <- sum(sapply(spb, function(z) sum(z[, "divergent__"])))
  sb <- smry(ab, colsb)
  cat("brms variables:", paste(posterior::variables(ab), collapse = " "),
      "\n")
  cat("names identical:", identical(cols, colsb), "\n")
  zm <- (sf[, "mean"] - sb[, "mean"]) /
    sqrt(sf[, "mcse_mean"]^2 + sb[, "mcse_mean"]^2)
  zs <- (sf[, "sd"] - sb[, "sd"]) /
    sqrt(sf[, "mcse_sd"]^2 + sb[, "mcse_sd"]^2)
  tab <- data.frame(weight = cols, mean_f = sf[, "mean"],
                    mean_b = sb[, "mean"], z_mean = zm, sd_f = sf[, "sd"],
                    sd_b = sb[, "sd"], z_sd = zs, rhat_f = sf[, "rhat"],
                    rhat_b = sb[, "rhat"], ess_f = sf[, "ess"],
                    ess_b = sb[, "ess"])
  cat("divergences frmtmb", divf, "brms", divb, "\n")
  print(format(tab, digits = 4), row.names = FALSE)
  # every non-simplex parameter too
  com <- intersect(posterior::variables(a), posterior::variables(ab))
  com <- setdiff(com, c(cols, "lp__", "lprior"))
  for (cn in com) {
    zf <- (mean(a[, , cn]) - mean(ab[, , cn])) /
      sqrt(posterior::mcse_mean(a[, , cn])^2 +
             posterior::mcse_mean(ab[, , cn])^2)
    cat(sprintf("  %-16s mean f %.4f b %.4f z %.2f | sd f %.4f b %.4f\n",
                cn, mean(a[, , cn]), mean(ab[, , cn]), zf,
                stats::sd(c(a[, , cn])), stats::sd(c(ab[, , cn]))))
  }
  res[[nm]] <- list(fit = fit, s = s, b = b)
  cat("\n-- brms summary --\n")
  print(summary(b))
  cat("\n-- frmtmb summary (ML fit) --\n")
  print(summary(fit))
  cat("\n-- frmtmb summary (draws) --\n")
  print(tryCatch(summary(s), error = function(e) conditionMessage(e)))
}
saveRDS(res, file.path(out, "sample-res.rds"))
