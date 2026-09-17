# Lane wt-conditions: candidate representative refusals, one per
# package, and their class vectors, on a named library arm.
#   Rscript dev/conditions-probe.R lane|base
arm <- commandArgs(trailingOnly = TRUE)[1L]
libs <- c(if (arm == "lane") "C:/Users/adf44/source/r/conditions-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/source/r/pinlib",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(libs)
suppressMessages({
  library(frmtmb); library(frmtmb.coupling); library(frmtmb.eam)
  library(frmtmb.latent); library(frmtmb.learn); library(frmtmb.ode)
  library(frmtmb.sample); library(frmtmb.spline)
})
show <- function(label, expr) {
  e <- tryCatch({ expr; NULL }, error = function(e) e)
  if (is.null(e)) { cat(label, ": NO ERROR\n"); return(invisible()) }
  cat(sprintf("%-10s %s\n           call: %s\n           msg: %s\n",
              label, paste(class(e), collapse = "/"),
              paste(deparse(conditionCall(e)), collapse = " "),
              substr(conditionMessage(e), 1, 90)))
}
d <- data.frame(y = c(1, 2, 3))
show("core", frm(y ~ 1, data = d, family = "not_a_family"))
show("coupling", frm_cross_spectrum(rnorm(20), rnorm(10)))
show("eam", wiener()$links[["ndt"]]$linkinv(0))
show("latent", hmm_starts(1))
show("learn", bandit4arm2_kalman_filter(subject = 1, sigma_o = -1))
tt <- c(0, 1, 2, 4)
show("ode", frm_lincmt(parms = list(ka = 1.1, ke = 0.2, V = 10), times = tt,
                       ncmt = 1, depot = TRUE,
                       events = data.frame(time = 2, state = "depot",
                                           value = 100,
                                           method = "replace")))
show("sample", frm_sample(y ~ 1, data = d, priors = NULL))
show("spline", frm_curve(data.frame(), newdata = data.frame()))
w <- tryCatch(frmtmb:::frm_warning("w"), warning = function(w) w)
cat("warning class:", class(w), "\n")
