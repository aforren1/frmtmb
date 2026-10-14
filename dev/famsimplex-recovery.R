# ML recovery of the simplex families on simulated data across seeds.
# Dirichlet data by gamma normalization, logistic-normal by an MVN draw on
# the additive log-ratio scale (dev/famsimplex-sim.R). Every warning a
# fit raises is recorded, which is the SE check's false-alarm count.
# Usage: Rscript dev/famsimplex-recovery.R <config> <first seed> <nseeds>
# Output: dev/famsimplex-recovery-<config>.rds and .log
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-famsimplex"
source(file.path(wt, "dev/famsimplex-sim.R"))
ns <- asNamespace("frmtmb")
args <- commandArgs(TRUE)
cfg <- args[1]
s0 <- as.integer(args[2])
ns_ <- as.integer(args[3])
truth <- list(
  dir = c(b_mub_Intercept = 0.4, b_mub_x = 0.6, b_muc_Intercept = -0.3,
          b_muc_x = -0.5, b_phi_Intercept = 2, b_phi_z = 0.5),
  dirre = c(b_mub_Intercept = 0.4, b_mub_x = 0.6, b_muc_Intercept = -0.3,
            b_muc_x = -0.5, phi = exp(2), sd_g__mub_Intercept = 0.5,
            sd_g__muc_Intercept = 0.3,
            cor_g__mub_Intercept__muc_Intercept = 0.4),
  ln = c(b_mub_Intercept = 0.4, b_mub_x = 0.6, b_muc_Intercept = -0.3,
         b_muc_x = -0.5, b_sigmab_Intercept = -0.5, b_sigmab_z = 0.4,
         sigmac = exp(-0.2), lncor__b__c = 0.3),
  lnre = c(b_mub_Intercept = 0.4, b_mub_x = 0.6, b_muc_Intercept = -0.3,
           b_muc_x = -0.5, sigmab = exp(-0.5), sigmac = exp(-0.2),
           lncor__b__c = 0.3, sd_g__mub_Intercept = 0.5,
           sd_g__muc_Intercept = 0.3,
           cor_g__mub_Intercept__muc_Intercept = 0.4),
  dir2 = c(b_mua_Intercept = 1, b_mua_x = 0.3, b_mub_Intercept = 1.5,
           b_mub_x = -0.4, b_muc_Intercept = 0.8, b_muc_x = 0.2))[[cfg]]
one <- function(seed) {
  d <- switch(cfg,
    dir = fs_data_dirichlet(seed, n = 300),
    dirre = fs_data_dirichlet(seed, n = 400, ng = 20, phi_formula = FALSE),
    ln = fs_data_logistic_normal(seed, n = 300),
    lnre = fs_data_logistic_normal(seed, n = 400, ng = 20,
                                   sigma_formula = FALSE),
    dir2 = fs_data_dirichlet2(seed, n = 300))
  form <- switch(cfg, dir = bf(Y ~ x, phi ~ z),
                 dirre = bf(Y ~ x + (1 | ID | g)),
                 ln = bf(Y ~ x, sigmab ~ z),
                 lnre = bf(Y ~ x + (1 | ID | g)), dir2 = bf(Y ~ x))
  fam <- switch(cfg, dir = , dirre = dirichlet(), ln = , lnre =
                  logistic_normal(), dir2 = brmsfamily("dirichlet2"))
  warns <- character(0)
  t0 <- proc.time()[["elapsed"]]
  fit <- withCallingHandlers(
    tryCatch(frm(form, family = fam, data = d),
             error = function(e) structure(conditionMessage(e),
                                           class = "fit_error")),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(fit, "fit_error")) {
    return(list(seed = seed, error = unclass(fit), warns = warns))
  }
  vo <- ns$hyp_vals_only(fit)
  est <- unlist(ns$hyp_env_vals(fit, vo$vals, vo$comp))[names(truth)]
  sm <- summary(fit)
  se <- rep(NA_real_, length(truth))
  names(se) <- names(truth)
  fx <- sm$fixed
  for (v in names(truth)) {
    r <- sub("^b_", "", v)
    if (startsWith(v, "b_") && r %in% rownames(fx)) {
      se[v] <- fx[r, "Est.Error"]
    } else if (v %in% rownames(sm$spec_pars)) {
      se[v] <- sm$spec_pars[v, "Est.Error"]
    } else if (startsWith(v, "lncor__")) {
      lab <- paste0("lncor(", sub("__", ",", sub("^lncor__", "", v)), ")")
      se[v] <- sm$spec_pars[lab, "Est.Error"]
    }
  }
  list(seed = seed, est = est, se = se, conv = fit$opt$convergence,
       warns = warns, time = el, ll = as.numeric(logLik(fit)))
}
res <- lapply(seq(s0, length.out = ns_), one)
saveRDS(list(cfg = cfg, truth = truth, res = res),
        file.path(wt, "dev", paste0("famsimplex-recovery-", cfg, ".rds")))
ok <- Filter(function(r) is.null(r$error), res)
cat("config", cfg, ": seeds", s0, "to", s0 + ns_ - 1, "; fitted", length(ok),
    "of", length(res), "; errors", length(res) - length(ok), "\n")
cat("nonzero convergence codes:", sum(vapply(ok, `[[`, 0, "conv") != 0), "\n")
nw <- vapply(res, function(r) length(r$warns), 0L)
cat("fits with a warning:", sum(nw > 0), "\n")
for (r in res[nw > 0]) cat("  seed", r$seed, ":", r$warns, sep = " ", "\n")
E <- do.call(rbind, lapply(ok, `[[`, "est"))
S <- do.call(rbind, lapply(ok, `[[`, "se"))
Z <- sweep(E, 2, truth) / S
tab <- data.frame(truth = truth, mean_est = colMeans(E),
                  bias = colMeans(E) - truth,
                  mc_se = apply(E, 2, stats::sd) / sqrt(nrow(E)),
                  sd_est = apply(E, 2, stats::sd),
                  mean_se = colMeans(S),
                  cover95 = colMeans(abs(Z) < stats::qnorm(0.975)))
print(tab, digits = 4)
cat("median seconds per fit:", stats::median(vapply(ok, `[[`, 0, "time")),
    "\n")
