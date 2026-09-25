# Nothing outside items 3.3 to 3.5 may move. Fit the models that
# existed before this lane on the BASE build (rellib-r3) and on this
# lane's build, and compare with identical().
#
# Usage: Rscript dev/phase3b-bitwise.R base|new   then   ... compare
# Output: dev/phase3b-log/bitwise-<tag>.rds, bitwise-compare.txt
args <- commandArgs(trailingOnly = TRUE)
tag <- args[[1]]
if (tag == "compare") {
  a <- readRDS("dev/phase3b-log/bitwise-base.rds")
  b <- readRDS("dev/phase3b-log/bitwise-new.rds")
  out <- c(sprintf("base: %s", paste(a$versions, collapse = " ")),
           sprintf("new:  %s", paste(b$versions, collapse = " ")))
  for (nm in names(a$res)) {
    x <- a$res[[nm]]; y <- b$res[[nm]]
    for (q in names(x)) {
      out <- c(out, sprintf("%-28s %-10s identical=%s", nm, q,
                            identical(x[[q]], y[[q]])))
    }
  }
  n_same <- sum(grepl("identical=TRUE", out))
  n_all <- sum(grepl("identical=", out))
  out <- c(out, sprintf("IDENTICAL %d of %d", n_same, n_all))
  writeLines(out, "dev/phase3b-log/bitwise-compare.txt")
  cat(out, sep = "\n")
  quit(save = "no")
}
lib <- if (tag == "base") character(0) else
  unique(c(Sys.getenv("P3B_LIB", "C:/Users/adf44/source/r/phase3b-lib"),
           "C:/Users/adf44/source/r/phase3b-lib"))
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb); library(frmtmb.eam); library(frmtmb.learn)
})
versions <- c(eam = find.package("frmtmb.eam"),
              learn = find.package("frmtmb.learn"))
take <- function(fit, sim = TRUE) {
  r <- list(loglik = as.numeric(logLik(fit)),
            par = fit$opt$par,
            fixef = unlist(fixef_by_dpar(fit)))
  f <- tryCatch(as.numeric(fitted(fit)[, 1]), error = function(e) NULL)
  if (!is.null(f)) r$fitted <- f
  if (sim) {
    s <- tryCatch({set.seed(3); unlist(simulate(fit, nsim = 2))},
                  error = function(e) NULL)
    if (!is.null(s)) r$sim <- s
  }
  r
}
res <- list()
set.seed(2)
d <- ddm_simulate(400, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
d$g <- rep(1:4, each = 100)
d$x <- rnorm(400)
res$wiener_plain <- take(frm(bf(rt | dec(upper) ~ x, bias = 0.5),
                             family = wiener(), data = d))
res$wiener_bias <- take(frm(bf(rt | dec(upper) ~ 1, bs ~ 1, bias ~ 1),
                            family = wiener(), data = d))
res$wiener_sv <- take(frm(bf(rt | dec(upper) ~ 1, bias = 0.5),
                          family = wiener(variability = "sv"), data = d),
                      sim = FALSE)
res$wiener_ndtgroup <- take(frm(bf(rt | dec(upper) + ndt_group(g) ~ 1,
                                   ndt ~ 1 + (1 | g), bias = 0.5),
                                family = wiener(), data = d))
res$wiener_maxndt <- take(frm(bf(rt | dec(upper) ~ 1, bias = 0.5),
                              family = wiener(max_ndt = 0.25), data = d))
res$wiener_vint <- take(frm(bf(rt | vint(upper) ~ 1, bias = 0.5),
                            family = wiener(), data = d))
dm <- d
dm$rt[1:10] <- runif(10, 0.05, 0.2)
res$wiener_mixture <- take(frm(bf(rt | dec(upper) ~ 1, bias1 = 0.5),
                               family = mixture(
                                 wiener(max_ndt = 0.4,
                                        allow_unreachable = TRUE),
                                 lognormal()),
                               data = dm), sim = FALSE)
# punch round 1 touched lba(), rdm() and gddm() (the refusals) and
# shares helpers with wiener_gng(), so each gets a fit here
dlb <- lba_simulate(300, v = c(2.4, 1.6), A = 0.5, k = 0.4, ndt = 0.2)
res$lba <- take(frm(bf(rt | vint(choice) ~ 1), family = lba(2), data = dlb),
                sim = FALSE)
drd <- rdm_simulate(300, v = c(3, 2), A = 0.5, k = 0.5, ndt = 0.2)
drd$code <- as.integer(drd$rt > stats::quantile(drd$rt, 0.8))
drd$y <- pmin(drd$rt, stats::quantile(drd$rt, 0.8))
res$rdm_cens <- take(frm(bf(y | vint(choice) + cens(code) ~ 1),
                         family = rdm(2), data = drd), sim = FALSE)
dg <- d
dg$responded <- as.integer(dg$rt < 1.2)
dg$rt[dg$responded == 0] <- 1.2
res$gng <- take(frm(bf(rt | dec(responded) ~ 1, bias = 0.5),
                    family = wiener_gng(deadline = 1.2), data = dg),
                sim = FALSE)
res$wiener_lpdf_export <- list(v = wiener_lpdf(c(0.2, 0.5, 1.4), 0.7, 1.3,
                                               0.4, c(0, 1, 1)))
# learn: every family, single session
dl <- frm_task_design("bandit2arm", n_subject = 6, n_trial = 40, seed = 5)
dl$choice <- frm_task_simulate(bandit2arm_delta(subject = id, trial = trial),
                               dl, pars = list(alpha = 0.4, tau = 3),
                               seed = 5)[[1]]$choice
fl <- bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1)
res$learn_delta <- take(frm(fl, family = bandit2arm_delta(subject = id,
                                                          trial = trial),
                            data = dl))
res$learn_delta_re <- take(frm(bf(choice | reward(pay1, pay2) ~ 1 + (1 | id),
                                  tau ~ 1),
                               family = bandit2arm_delta(subject = id,
                                                         trial = trial),
                               data = dl))
res$learn_dual <- take(frm(bf(choice | reward(pay1, pay2) ~ 1, Apun ~ 1,
                              tau ~ 1),
                           family = bandit2arm_dual(subject = id,
                                                    trial = trial),
                           data = dl))
res$learn_prl <- take(frm(bf(choice | reward(pay1, pay2) ~ 1, bias ~ 1,
                             tau ~ 1),
                          family = prl_fictitious(subject = id,
                                                  trial = trial),
                          data = dl))
di <- frm_task_design("igt", n_subject = 5, n_trial = 50, seed = 7)
di$choice <- frm_task_simulate(igt_pvl_delta(subject = id, trial = trial), di,
                               pars = list(alpha = 0.3, shape = 0.5,
                                           lambda = 1, tau = 1),
                               seed = 7)[[1]]$choice
res$learn_pvl <- take(frm(bf(choice | payoff(pay1, pay2, pay3, pay4) ~ 1,
                             shape ~ 1,
                             lambda ~ 1, tau ~ 1),
                          family = igt_pvl_delta(subject = id, trial = trial),
                          data = di))
if (requireNamespace("RWiener", quietly = TRUE)) {
  dr <- frm_task_design("bandit2arm", n_subject = 5, n_trial = 30, seed = 62)
  dr <- frm_task_simulate(
    rlddm(subject = id, trial = trial), dr,
    pars = list(alpha = 0.4, drift = 3, bs = 1.6, ndt = 0.2, bias = 0.5),
    seed = 62)[[1]]
  res$learn_rlddm <- take(frm(bf(rt | dec(choice) + reward(pay1, pay2) ~ 1,
                                 drift ~ 1, bs ~ 1, ndt ~ 1, bias ~ 1),
                              family = rlddm(subject = id, trial = trial),
                              data = dr), sim = FALSE)
}
res$learn_trace <- list(tr = frm_value_trace(frm(
  fl, family = bandit2arm_delta(subject = id, trial = trial), data = dl)))
res$learn_tasksim <- list(y = frm_task_simulate(
  ts_par7(subject = id, trial = trial),
  frm_task_design("twostep", n_subject = 3, n_trial = 20, seed = 8),
  pars = list(w = 0.5, alpha1 = 0.4, tau1 = 2, alpha2 = 0.4, tau2 = 2,
              lambda = 0.5, pers = 0.1), seed = 8))
saveRDS(list(versions = versions, res = res),
        sprintf("dev/phase3b-log/bitwise-%s.rds", tag))
cat("done", tag, length(res), "\n")
