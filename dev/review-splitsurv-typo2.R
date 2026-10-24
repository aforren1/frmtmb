# simplex_pars = naming an extra parameter the family does not have,
# and a mv model with a declared simplex. Lane lib.
.libPaths(c("C:/Users/adf44/source/r/wt-splitsurv-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
fam <- function(sp) frmtmb_family(
  "tsimplex", dpars = "mu", links = list(mu = "identity"),
  lpdf = function(y, dpars, aterms, extra) {
    w <- mo_simplex(extra$tw_raw, aterms[["tw_chart"]])
    cen <- c(-3, 0, 3)
    f <- 0
    for (j in 1:3) f <- f + w[j] * exp(-0.5 * (y - dpars[["mu"]] - cen[j])^2)
    log(f) - 0.5 * log(2 * pi)
  },
  init_dpars = list(mu = function(y, aterms) 0),
  extra_pars = function(y, aterms) list(tw_raw = rep(0, 2)),
  simplex_pars = sp)
set.seed(5)
k <- sample(1:2, 300, replace = TRUE, prob = c(0.6, 0.4))
dd <- data.frame(y = 0.5 + c(-3, 0)[k] + rnorm(300))
show <- function(lab, expr) {
  w <- character()
  r <- tryCatch(withCallingHandlers({
    v <- expr
    if (inherits(v, "frmtmb_fit")) paste("FIT logLik", format(as.numeric(logLik(v)), digits = 10))
    else if (is.data.frame(v)) paste("TABLE classes:", paste(unique(v$class), collapse = ","))
    else "ok"
  }, warning = function(z) { w <<- c(w, conditionMessage(z)); invokeRestart("muffleWarning") }),
  error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-28s %s%s\n", lab, substr(r, 1, 300),
              if (length(w)) paste0("  [warnings: ", paste(substr(w, 1, 120), collapse = " | "), "]") else ""))
}
show("fit, correct declaration", frm(y ~ 1, family = fam(c(tw = "tw_raw")), data = dd))
show("fit, typo in declaration", frm(y ~ 1, family = fam(c(tw = "tw_rwa")), data = dd))
show("get_prior, typo", get_prior(y ~ 1, data = dd, family = fam(c(tw = "tw_rwa"))))
show("prior on typo class", frm(y ~ 1, family = fam(c(tw = "tw_rwa")), data = dd,
                               prior = set_prior("dirichlet(2)", class = "tw")))
show("sbhaz name, 3rd-party", frm(y ~ 1, family = fam(c(sbhaz = "tw_raw")), data = dd,
                                  prior = set_prior("dirichlet(2)", class = "sbhaz")))
fam0 <- frmtmb_family("tzero", dpars = "mu", links = list(mu = "identity"),
  lpdf = function(y, dpars, aterms, extra) -0.5 * (y - dpars[["mu"]])^2,
  extra_pars = function(y, aterms) list(tw_raw = numeric(0)),
  simplex_pars = c(tw = "tw_raw"))
show("no coordinate", frm(y ~ 1, family = fam0, data = dd))
show("get_prior, no coordinate", get_prior(y ~ 1, data = dd, family = fam0))
show("get_prior, correct", get_prior(y ~ 1, data = dd, family = fam(c(tw = "tw_raw"))))
