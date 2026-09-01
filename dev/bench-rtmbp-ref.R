## Reference pass (RTMB only):
##  - frm() fit of the InstEval model, for the reference logLik
##  - hand-rolled RTMB fit of shape A and shape B, to (i) check the
##    hand-rolled optimum against frm() and (ii) freeze an optimum
##    parameter vector that every later configuration reuses, so all
##    per-call timings and identity checks happen at identical inputs.
.libPaths(c("C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/libp",
            .libPaths()))
suppressPackageStartupMessages({
  library(RTMB); library(Matrix); library(lme4)
})
source("dev/bench-rtmbp-data.R")

res <- list()

## --- frm() reference on InstEval --------------------------------------
suppressPackageStartupMessages(library(frmtmb))
t0 <- proc.time()[["elapsed"]]
fit <- suppressWarnings(frm(bf(y ~ service + (1 | s) + (1 | d)) + gaussian(),
                            data = lme4::InstEval, se = FALSE))
res$frm_wall <- proc.time()[["elapsed"]] - t0
res$frm_logLik <- as.numeric(stats::logLik(fit))
res$frm_coef <- stats::coef(fit)
res$frm_sigma <- tryCatch(sigma(fit), error = function(e) NA_real_)
res$frm_vc <- tryCatch(VarCorr(fit), error = function(e) NULL)
cat(sprintf("frm() logLik = %.6f  (wall %.1f s)\n", res$frm_logLik, res$frm_wall))
print(res$frm_coef)

## --- hand-rolled fits --------------------------------------------------
for (shape in c("A", "B")) {
  dat <- build_shape(shape)
  f <- make_f(shape, dat)
  obj <- MakeADFun(f, dat$parameters, random = dat$random, silent = TRUE)
  t0 <- proc.time()[["elapsed"]]
  op <- nlminb(obj$par, obj$fn, obj$gr)
  w <- proc.time()[["elapsed"]] - t0
  cat(sprintf("shape %s: nll = %.6f  logLik = %.6f  wall %.1f s  conv %d\n",
              shape, op$objective, -op$objective, w, op$convergence))
  print(round(op$par, 6))
  res[[paste0("opt_", shape)]] <- op$par
  res[[paste0("nll_", shape)]] <- op$objective
  res[[paste0("wall_", shape)]] <- w
  res[[paste0("iter_", shape)]] <- op$iterations
  res[[paste0("eval_", shape)]] <- op$evaluations
}

cat(sprintf("\nlogLik check shape A: frm %.6f  hand-rolled %.6f  diff %.3e\n",
            res$frm_logLik, -res$nll_A, abs(res$frm_logLik - (-res$nll_A))))

saveRDS(res, "dev/rtmbp-ref.rds")
cat("wrote dev/rtmbp-ref.rds\n")
