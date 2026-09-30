# Reviewer re-check (punch round 1), lane ordinal: B3's empty tau_raw path
# and the hidden disc in par_template(). Every level of thres(gr = h) at
# two categories under sum_to_zero (no threshold parameter), with and
# without disc; the ungrouped two-category sum_to_zero model; then every
# reader, the template round trip, and frmtmb.sample. Also start = and
# newparams = from par_template() on a flexible ordinal fit, whose template
# now carries the mapped disc. Data seed 20261015, sampler seed 5.
# Output: dev/ordinal-rev2-log-empty.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(20261015)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n),
                h = factor(sample(c("p", "q"), n, TRUE)))
d$yb <- 1L + (stats::rlogis(n) / exp(0.3 * d$z) + d$x > 0)
d$y <- 1L + findInterval(stats::rlogis(n) + d$x, c(-1, 0, 1))
run <- function(lab, expr) {
  cat("\n##", lab, "\n")
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("  WARNING:", substr(conditionMessage(w), 1, 160), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))
  invisible(r)
}
shapes <- list(
  gr_stz = quote(frm(yb | thres(gr = h) ~ x, data = d,
                     family = cumulative(threshold = "sum_to_zero"))),
  gr_stz_disc = quote(frm(bf(yb | thres(gr = h) ~ x, disc ~ 0 + z), data = d,
                          family = sratio(threshold = "sum_to_zero"))),
  ungr_stz = quote(frm(yb ~ x, data = d,
                       family = acat(threshold = "sum_to_zero"))))
for (nm in names(shapes)) {
  cat("\n==========", nm, "\n")
  f <- tryCatch(suppressWarnings(eval(shapes[[nm]])), error = function(e) e)
  if (inherits(f, "error")) { cat("ERROR:", conditionMessage(f), "\n"); next }
  run("par_template lengths", print(lengths(unclass(par_template(f)))))
  run("variables", print(variables(f)))
  run("fixef", print(fixef(f)))
  run("fixef(flatten = TRUE)", print(fixef(f, flatten = TRUE)))
  run("coef", print(coef(f)))
  run("confint rows", print(rownames(confint(f))))
  run("vcov rows", print(rownames(vcov(f))))
  run("brms_par_labels", print(brms_par_labels(f)))
  run("hypothesis b_Intercept[..,1] = 0 (class NULL)", {
    v <- grep("^b_Intercept", variables(f), value = TRUE)[1]
    print(hypothesis(f, paste0("`", v, "` = 0"), class = NULL)$hypothesis[, 1:3])
  })
  run("summary", print(summary(f)))
  run("refit from par_template(start =)", {
    g <- frm(formula(f), data = d, family = family(f), start = par_template(f))
    cat("  logLik identical:", identical(as.numeric(logLik(g)),
                                        as.numeric(logLik(f))),
        " rel diff", abs(as.numeric(logLik(g)) - as.numeric(logLik(f))) /
          abs(as.numeric(logLik(f))), "\n")
  })
  run("frm_simulate(newparams = par_template())", {
    s <- frm_simulate(f, newparams = par_template(f), nsim = 1, seed = 1)
    cat("  class", class(s)[1], "\n")
  })
  run("update() refit", cat("  logLik", as.numeric(logLik(update(f, data = d[1:250, ]))), "\n"))
  run("bootstrap", print(dim(frm_bootstrap(f, nsim = 3, seed = 2)$t)))
  run("frm_sample", {
    ds <- suppressMessages(frm_sample(f, chains = 1, iter = 100, refresh = 0,
                                      seed = 5))
    print(variables(ds))
    ep <- posterior_epred(ds, ndraws = 2)
    cat("  epred dim", dim(ep), " log_lik finite",
        all(is.finite(log_lik(ds))), "\n")
    print(fixef(ds))
  })
}
cat("\n========== flexible ordinal: template round trip with the mapped disc\n")
f <- frm(y ~ x, family = cumulative(), data = d)
run("par_template(f)", print(par_template(f)))
run("frm(start = par_template(f))", {
  g <- frm(y ~ x, family = cumulative(), data = d, start = par_template(f))
  cat("  logLik identical:", identical(as.numeric(logLik(g)),
                                      as.numeric(logLik(f))), "\n")
})
run("start = par_template with disc_(Intercept) edited to 0.5", {
  tp <- par_template(f); tp$betad[] <- 0.5
  g <- frm(y ~ x, family = cumulative(), data = d, start = tp)
  cat("  logLik", as.numeric(logLik(g)), "vs", as.numeric(logLik(f)), "\n")
})
