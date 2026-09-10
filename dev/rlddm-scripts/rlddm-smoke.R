# The CONTROL: a model without ndt_group() must be the model 0.3.0
# fitted, bit for bit.
#
#   Rscript dev/rlddm-scripts/rlddm-smoke.R <lib> <out.rds>
#
# Seed 4242 throughout. Run once against the round's shared reference
# build of the base commit and once against this lane's library, then
# compare the two RDS files with identical().

args <- commandArgs(trailingOnly = TRUE)
lib <- args[[1L]]
out <- args[[2L]]

.libPaths(c(lib,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.eam)
  library(frmtmb.learn)
})
cat("lib  :", lib, "\n")
cat("learn:", format(packageVersion("frmtmb.learn")), "at",
    dirname(system.file(package = "frmtmb.learn")), "\n")

seed <- 4242L
d <- frm_task_design("bandit2arm", n_subject = 6, n_trial = 60,
                     seed = seed)
s <- frm_task_simulate(rlddm(subject = id, trial = trial), d,
                       pars = list(alpha = 0.4, drift = 3, bs = 1.6,
                                   ndt = 0.2, bias = 0.5),
                       seed = seed)[[1L]]
got <- list()
got$n <- nrow(s)
got$min_rt <- min(s$rt)
got$sim_rt <- s$rt
got$sim_choice <- s$choice

# ---- 1. the plain model, which must not move
f0 <- bf(rt | dec(choice) + reward(pay1, pay2) ~ 1 + (1 | id),
         drift ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5)
fit0 <- frm(f0, family = rlddm(subject = id, trial = trial), data = s,
            se = TRUE)
got$ll_plain <- as.numeric(stats::logLik(fit0))
got$coef_plain <- unlist(fixef(fit0))
got$sd_plain <- as.numeric(fit0$sdr$sd)
got$vcov_plain <- as.numeric(stats::vcov(fit0))
got$ci_plain <- as.matrix(suppressWarnings(stats::confint(fit0)))
got$varcorr_plain <- as.numeric(VarCorr(fit0)[[1L]])
got$ndt_response <- as.numeric(suppressWarnings(stats::predict(
  fit0, dpar = "ndt", type = "response")))
got$eta_plain <- as.numeric(suppressWarnings(stats::predict(
  fit0, dpar = "drift", type = "link")))
got$trace <- as.matrix(frm_value_trace(fit0)[,
  c("q1", "q2", "drift_t", "pe", "dens")])
got$conv_plain <- fit0$opt$convergence
got$par_plain <- as.numeric(fit0$opt$par)

# ---- 2. a fixed-parameter probe of the likelihood, which is the
# density rather than an optimizer path
probe <- function(ndt) {
  fp <- bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift = 3,
           bs = 1.6, ndt = ndt, bias = 0.5)
  o <- tryCatch(frm(fp, family = rlddm(subject = id, trial = trial),
                    data = s, dry_run = "objective"),
                error = function(e) conditionMessage(e))
  if (is.character(o)) return(o)
  as.numeric(o$obj$fn(o$obj$par))
}
got$pin_02 <- probe(0.2)
got$pin_015 <- probe(0.15)

# ---- 3. a max_ndt model, where the bound is stated up front
f_mx <- bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift ~ 1,
           bs ~ 1, ndt ~ 1, bias = 0.5)
fit_mx <- frm(f_mx, family = rlddm(subject = id, trial = trial,
                                   max_ndt = min(s$rt)), data = s)
got$ll_maxndt <- as.numeric(stats::logLik(fit_mx))
got$coef_maxndt <- unlist(fixef(fit_mx))

# ---- 4. the refusals
got$err_maxndt_high <- tryCatch({
  frm(f_mx, family = rlddm(subject = id, trial = trial,
                           max_ndt = max(s$rt) + 1), data = s)
  "NO ERROR"
}, error = function(e) conditionMessage(e))

# ---- 5. does a bare wiener() score outside frm()? The round-two
# review's fifth-entry question, answered rather than argued.
fam_w <- wiener()
got$wiener_lpdf_bare <- tryCatch(
  as.numeric(fam_w$lpdf(c(0.4, 0.9),
                        list(mu = 1.2, bs = 1.5, ndt = 0.15, bias = 0.3),
                        list(dec = c(1, 0)))),
  error = function(e) paste("ERR:", conditionMessage(e)))
got$wiener_link_bare <- tryCatch(fam_w$links$ndt$linkinv(0),
                                 error = function(e) "ERR")

saveRDS(got, out)
cat("\nlogLik plain      :", sprintf("%.12f", got$ll_plain), "\n")
cat("coef plain        :",
    paste(names(got$coef_plain), sprintf("%.12f", got$coef_plain),
          collapse = "; "), "\n")
cat("ndt response[1:3] :", sprintf("%.12f", got$ndt_response), "\n")
cat("pinned ndt=0.2    :", if (is.character(got$pin_02)) got$pin_02 else
  sprintf("%.12f", got$pin_02), "\n")
cat("pinned ndt=0.15   :", if (is.character(got$pin_015)) got$pin_015 else
  sprintf("%.12f", got$pin_015), "\n")
cat("logLik max_ndt    :", sprintf("%.12f", got$ll_maxndt), "\n")
cat("max_ndt too high  :", substr(got$err_maxndt_high, 1, 90), "\n")
cat("wiener()$lpdf bare:",
    if (is.character(got$wiener_lpdf_bare)) got$wiener_lpdf_bare else
      sprintf("%.12f", got$wiener_lpdf_bare), "\n")
cat("wiener() linkinv  :", format(got$wiener_link_bare), "\n")
cat("wrote", out, "\n")
