# Reviewer: frm_sample() on each simplex family; variables and priors.
# Output: dev/famsimplex-rev-sample.log
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/famsimplex-sim.R")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev"
smp <- function(lab, fit, ...) {
  cat("\n==========", lab, "\n")
  msg <- character()
  w <- character()
  ds <- withCallingHandlers(
    tryCatch(frm_sample(fit, chains = 1, iter = 400, refresh = 0, seed = 3,
                        ...),
             error = function(e) {
               cat("ERROR:", conditionMessage(e), "\n")
               NULL
             }),
    message = function(m) {
      msg <<- c(msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  cat("MESSAGES:\n", paste(msg, collapse = "\n"), "\n")
  if (length(w)) cat("WARNINGS:", unique(w), sep = "\n  ")
  if (!is.null(ds)) {
    cat("variables:", variables(ds), "\n")
    s <- posterior::summarise_draws(posterior::as_draws_df(ds), "median",
                                    "sd")
    print(as.data.frame(s), digits = 4)
  }
  invisible(ds)
}
b_ln <- readRDS(file.path(wt, "famsimplex-brms-ln.rds"))
cat("brms ln variables:", b_ln$variables, "\n")
fl <- frm(bf(Y ~ x, sigmab ~ z), family = logistic_normal(), data = b_ln$data)
dl <- smp("ln sigmab ~ z", fl)
dl2 <- smp("ln with lkj(4) on lncor", fl,
           prior = set_prior("lkj(4)", class = "lncor"))
b_d2 <- readRDS(file.path(wt, "famsimplex-brms-dir2.rds"))
cat("brms dir2 variables:", b_d2$variables, "\n")
f2 <- frm(bf(Y ~ x), family = brmsfamily("dirichlet2"), data = b_d2$data)
d2 <- smp("dirichlet2", f2)
if (!is.null(d2)) {
  pp <- posterior_predict(d2, ndraws = 5)
  cat("dir2 pp dims:", dim(pp), " max |rowsum - 1|:",
      max(abs(apply(pp, 1:2, sum) - 1)), "\n")
  ep <- posterior_epred(d2, ndraws = 5)
  cat("dir2 epred dims:", dim(ep), dimnames(ep)[[3]], "\n")
}
b_dre <- readRDS(file.path(wt, "famsimplex-brms-dirre.rds"))
cat("brms dirre variables:", b_dre$variables, "\n")
fd <- frm(bf(Y ~ x + (1 | ID | g)), family = dirichlet(), data = b_dre$data)
dd <- smp("dirichlet (1 | ID | g)", fd)
dg <- smp("dirichlet with gamma(0.01, 0.01) on phi", fd,
          prior = set_prior("gamma(0.01, 0.01)", class = "phi"))
cat("\nfrm() MAP with gamma(0.01, 0.01) on phi:\n")
print(tryCatch(summary(frm(bf(Y ~ x), family = dirichlet(),
                           data = b_dre$data,
                           prior = set_prior("gamma(0.01, 0.01)",
                                             class = "phi")))$spec_pars,
               error = function(e) conditionMessage(e)))
