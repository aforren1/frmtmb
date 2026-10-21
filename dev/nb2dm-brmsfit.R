# Fit agreement with brms 2.23.0 on the same data: negbinomial2 and
# dirichlet_multinomial, one fixed-effects and one random-intercept model
# each. Three arms per model: brms (Stan, 4 chains x 2000, seed 101),
# frmtmb's ML fit, and frmtmb.sample's frm_sample() of that fit under
# brms's default priors (4 chains x 2000, seed 101). Each estimate is
# reported in units of brms's posterior sd; the two posteriors are also
# compared against the MCSE of the difference.
#   Rscript dev/nb2dm-brmsfit.R <case> > dev/nb2dm-out/brmsfit-<case>.txt
# case: nb2, nb2re, dm, dmre
.libPaths(c("C:/Users/adf44/source/r/wt-nb2dm-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache",
                                               winslash = "/"))
}
case <- commandArgs(trailingOnly = TRUE)[1]
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", dirname(system.file(package = "frmtmb")),
    format(packageVersion("frmtmb")), "sample",
    format(packageVersion("frmtmb.sample")), "brms",
    format(packageVersion("brms")), "case", case, "\n")

set.seed(2026)
n <- 400
d <- data.frame(x = rnorm(n), g = gl(20, 20))
u <- rnorm(20, 0, 0.5)[d$g]
if (startsWith(case, "nb2")) {
  eta <- 0.8 + 0.5 * d$x + if (case == "nb2re") u else 0
  d$y <- rnbinom(n, mu = exp(eta), size = 1 / 0.4)
  f_frm <- if (case == "nb2re") y ~ x + (1 | g) else y ~ x
  fam_frm <- frmtmb::brmsfamily("negbinomial2")
  fam_brm <- brms::brmsfamily("negbinomial2")
  map <- c(b_Intercept = "Intercept", b_x = "x", sigma = "sigma",
           sd_g__Intercept = "sd_g__Intercept")
} else {
  K <- 3
  tr <- sample(5:40, n, TRUE)
  e2 <- 0.4 + 0.6 * d$x + if (case == "dmre") u else 0
  e3 <- -0.3 - 0.5 * d$x + if (case == "dmre") rnorm(20, 0, 0.4)[d$g] else 0
  P <- exp(cbind(0, e2, e3))
  P <- P / rowSums(P)
  phi <- 6
  Y <- t(vapply(seq_len(n), function(i) {
    p <- rgamma(K, P[i, ] * phi)
    as.vector(rmultinom(1, tr[i], p / sum(p)))
  }, numeric(K)))
  colnames(Y) <- c("a", "b", "c")
  d$Y <- Y
  d$tr <- tr
  f_frm <- if (case == "dmre") Y | trials(tr) ~ x + (1 | g) else
    Y | trials(tr) ~ x
  fam_frm <- frmtmb::dirichlet_multinomial()
  fam_brm <- brms::dirichlet_multinomial()
  map <- c(b_mub_Intercept = "mub_Intercept", b_mub_x = "mub_x",
           b_muc_Intercept = "muc_Intercept", b_muc_x = "muc_x", phi = "phi",
           sd_g__mub_Intercept = "sd_g__mub_Intercept",
           sd_g__muc_Intercept = "sd_g__muc_Intercept")
}

t0 <- proc.time()[["elapsed"]]
fit <- frm(f_frm, family = fam_frm, data = d)
cat("frmtmb ML fit:", round(proc.time()[["elapsed"]] - t0, 2), "s, logLik",
    format(logLik(fit), digits = 12), "convergence",
    fit$opt$convergence %||% NA, "\n")
t0 <- proc.time()[["elapsed"]]
dr <- frm_sample(fit, chains = 4, iter = 2000, seed = 101, refresh = 0)
cat("frm_sample:", round(proc.time()[["elapsed"]] - t0, 1), "s\n")
t0 <- proc.time()[["elapsed"]]
rds <- file.path("dev/nb2dm-out", paste0("brmsfit-", case, ".rds"))
bfit <- if (file.exists(rds)) readRDS(rds) else {
  b <- brms::brm(f_frm, family = fam_brm, data = d, chains = 4,
                 iter = 2000, seed = 101, refresh = 0, silent = 2,
                 init = 0)
  saveRDS(b, rds)
  b
}
cat("brms:", round(proc.time()[["elapsed"]] - t0, 1), "s\n")

print(summary(bfit))
bd <- posterior::as_draws_df(bfit)
fd <- posterior::as_draws_df(dr)
cat("frm_sample variables:", paste(head(posterior::variables(fd), 20),
                                   collapse = " "), "\n")
# the ML point, on brms's names
ml <- c()
fx <- fixef(fit)
for (b in names(map)) {
  nm <- map[[b]]
  v <- if (nm %in% rownames(fx)) fx[nm, "Estimate"] else NA_real_
  if (is.na(v) && nm %in% c("sigma", "phi")) {
    lp <- Filter(function(l) identical(l$dpar, nm), fit$frame$linpreds)[[1]]
    v <- lp$link$linkinv(fit$estimates[[lp$par]][lp$idx])
  }
  if (is.na(v) && startsWith(nm, "sd_")) {
    vc <- tryCatch(VarCorr(fit), error = function(e) NULL)
    v <- tryCatch({
      s <- vc$g$sd
      key <- sub("^sd_g__", "", nm)
      if (is.matrix(s)) s[key, "Estimate"] else s[[1]]
    }, error = function(e) NA_real_)
  }
  ml[b] <- if (length(v) == 1L) v else NA_real_
}
rows <- lapply(names(map), function(b) {
  if (!b %in% names(bd)) return(NULL)
  bx <- bd[[b]]
  fv <- if (b %in% names(fd)) fd[[b]] else rep(NA_real_, 2)
  mcse_b <- posterior::mcse_mean(bx)
  mcse_f <- if (all(is.na(fv))) NA else posterior::mcse_mean(fv)
  data.frame(par = b, brms_mean = mean(bx), brms_sd = sd(bx),
             ml = ml[[b]], ml_z = (ml[[b]] - mean(bx)) / sd(bx),
             sample_mean = mean(fv), sample_sd = sd(fv),
             diff_z_sd = (mean(fv) - mean(bx)) / sd(bx),
             diff_z_mcse = (mean(fv) - mean(bx)) / sqrt(mcse_b^2 + mcse_f^2),
             sd_ratio = sd(fv) / sd(bx))
})
tab <- do.call(rbind, rows)
print(tab, digits = 4, row.names = FALSE)
cat("max |ml - brms mean| / brms sd:", format(max(abs(tab$ml_z), na.rm = TRUE),
                                             digits = 3), "\n")
cat("max |sample - brms mean| / brms sd:",
    format(max(abs(tab$diff_z_sd), na.rm = TRUE), digits = 3), "\n")
cat("max |sample - brms mean| / mcse:",
    format(max(abs(tab$diff_z_mcse), na.rm = TRUE), digits = 3), "\n")
cat("brms rhat max:", format(max(brms::rhat(bfit), na.rm = TRUE), digits = 4),
    " divergences:", sum(brms::nuts_params(bfit)$Value[
      brms::nuts_params(bfit)$Parameter == "divergent__"]), "\n")
