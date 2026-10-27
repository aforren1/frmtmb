# Reviewer, fix round 1, item 2: what frmtmb REPORTS as `sigma` under a
# cov = TRUE ar(), against brms's `sigma` on the same data (the AR(1)
# 0.9 data of dev/acprior-agree.R). brms's reference is its own draws
# from the lane's run (dev/acprior-out/agree-brms.rds, seed 3,
# 4 x 4000, brms defaults). frmtmb: the ML fit and frm_sample (4 x 2000,
# seed 3), every reporter that names sigma.
# Usage: Rscript dev/review-acprior-r1-report.R
# Log: dev/review-acprior-log/r1-report.txt
wt <- "C:/Users/adf44/source/r/frmtmb-wt-acprior"
.libPaths(c("C:/Users/adf44/source/r/wt-acprior-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
Sys.setenv(FRMTMB_STAN_CACHE = file.path(wt, "dev", "stan-cache"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
src <- readLines(file.path(wt, "dev", "acprior-agree.R"))
i0 <- grep("^make_data <- function", src)
i1 <- i0 + which(src[(i0 + 1):length(src)] == "}")[1]
eval(parse(text = src[i0:i1]))
da <- make_data()$ar
R <- readRDS(file.path(wt, "dev", "acprior-out", "agree-brms.rds"))$ar
cat(sprintf("brms sigma (innovation sd): mean %.4f sd %.4f; ar mean %.4f\n",
            mean(R$sigma), sd(R$sigma), mean(R$par)))
cat(sprintf("brms draws carried to the marginal sd sigma / sqrt(1 - ar^2): mean %.4f\n",
            mean(R$sigma / sqrt(1 - R$par^2))))
fit <- frm(y ~ x + ar(time = t, gr = g, cov = TRUE), data = da)
show <- function(lab, expr) {
  cat("\n---", lab, "\n")
  out <- tryCatch(utils::capture.output(print(expr)),
                  error = function(e) paste("ERROR:", conditionMessage(e)))
  hit <- grep("sigma|Residual|ar", out, value = TRUE)
  cat(if (length(hit)) hit else out[seq_len(min(15, length(out)))], sep = "\n")
}
show("ML summary(fit)", summary(fit))
show("ML fixef(fit)", fixef(fit))
show("ML VarCorr(fit)", VarCorr(fit))
show("ML sigma(fit)", stats::sigma(fit))
th <- fit$estimates$thetaac
phi <- unlist(frmtmb:::autocor_natural(th, fit$frame$autocor[[1]]))
cat(sprintf("\nML ar %.4f; sqrt(1 - ar^2) %.4f\n", phi, sqrt(1 - phi^2)))
ds <- suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 4, iter = 2000, seed = 3, refresh = 0)))
show("frm_sample variables(ds)", grep("sigma", posterior::variables(ds),
                                       value = TRUE))
show("frm_sample posterior_summary(ds)", posterior_summary(ds))
show("frm_sample summary(ds)", summary(ds))
mm <- as.matrix(ds)
cat(sprintf("\nframe draws 'sigma' mean %.4f; carried to innovation sd %.4f\n",
            mean(mm[, "sigma"]), mean(mm[, "sigma"] * sqrt(1 - mm[, "ar[1]"]^2))))
