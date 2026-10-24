# A third-party family that declares a simplex (simplex_pars =), on
# frmtmb.sample's sampler and through the prior tables. Lane lib only.
.libPaths(c("C:/Users/adf44/source/r/wt-splitsurv-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb", format(packageVersion("frmtmb")), dirname(find.package("frmtmb")),
    "; sample", dirname(find.package("frmtmb.sample")), "\n")
tsimplex_family <- function(centers = c(-3, 0, 3)) {
  frmtmb_family(
    "tsimplex", dpars = "mu", links = list(mu = "identity"),
    lpdf = function(y, dpars, aterms, extra) {
      w <- mo_simplex(extra$tw_raw, aterms[["tw_chart"]])
      cen <- aterms[["tw_centers"]]
      f <- 0
      for (j in seq_along(cen)) {
        f <- f + w[j] * exp(-0.5 * (y - dpars[["mu"]] - cen[j])^2)
      }
      log(f) - 0.5 * log(2 * pi)
    },
    init_dpars = list(mu = function(y, aterms) 0),
    extra_pars = function(y, aterms) {
      list(tw_raw = rep(0, length(aterms[["tw_centers"]]) - 1L))
    },
    aterm_data = function(y, aterms) list(tw_centers = centers),
    simplex_pars = c(tw = "tw_raw"))
}
set.seed(5)
k <- sample(1:2, 300, replace = TRUE, prob = c(0.6, 0.4))
dd <- data.frame(y = 0.5 + c(-3, 0)[k] + rnorm(300))
fit <- frm(y ~ 1, family = tsimplex_family(), data = dd)
cat("ML weights:", format(mo_simplex(fit$estimates$tw_raw), digits = 4), "\n")
cat("== default_prior(route sample)\n")
print(tryCatch(frmtmb.sample::frm_sample_prior(fit), error = function(e)
  conditionMessage(e)))
msgs <- character()
ds <- withCallingHandlers(
  frm_sample(fit, chains = 2, iter = 600, seed = 1, refresh = 0),
  message = function(m) { msgs <<- c(msgs, conditionMessage(m))
                          invokeRestart("muffleMessage") },
  warning = function(w) { msgs <<- c(msgs, paste("W:", conditionMessage(w)))
                          invokeRestart("muffleWarning") })
cat("== messages from frm_sample (default prior)\n"); cat(msgs, sep = "\n")
v <- posterior::variables(ds)
cat("variables:", paste(v, collapse = " "), "\n")
W <- as.matrix(ds, variable = paste0("tw[", 1:3, "]"))
cat("draw means:", format(colMeans(W), digits = 4), "\n")
msgs <- character()
ds2 <- withCallingHandlers(
  frm_sample(fit, chains = 2, iter = 600, seed = 1, refresh = 0,
             prior = set_prior("dirichlet(20)", class = "tw")),
  message = function(m) { msgs <<- c(msgs, conditionMessage(m))
                          invokeRestart("muffleMessage") },
  warning = function(w) { msgs <<- c(msgs, paste("W:", conditionMessage(w)))
                          invokeRestart("muffleWarning") })
cat("== messages from frm_sample (dirichlet(20) on tw)\n"); cat(msgs, sep = "\n")
W2 <- as.matrix(ds2, variable = paste0("tw[", 1:3, "]"))
cat("draw means under dirichlet(20):", format(colMeans(W2), digits = 4), "\n")
# update() of a fit with a declared simplex
up <- tryCatch(update(fit, prior = set_prior("dirichlet(2)", class = "tw")),
               error = function(e) conditionMessage(e))
cat("update with tw prior:", if (is.character(up)) up else
  format(mo_simplex(up$estimates$tw_raw), digits = 4), "\n")
# summary() and its simplex table
s <- tryCatch(capture.output(summary(fit)), error = function(e) conditionMessage(e))
cat("== summary(fit)\n"); cat(s, sep = "\n")
