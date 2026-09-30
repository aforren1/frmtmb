# Reviewer, lane ordinal: the grouped equidistant MAP row of
# dev/ordinal-rev-map.R, where brms's gradient at frmtmb's MAP has an
# entry of 1. Which parameter, and why. Data seed 20261004 (as there).
# Output: dev/ordinal-rev-log-map-gr.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
sp <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad"
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = file.path(sp, "ordrev-stan-cache"))
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordinal"
env <- new.env(parent = asNamespace("frmtmb"))
sys.source(file.path(wt, "tests/testthat/helper-brms.R"), envir = env)
ctl <- frmtmb_control(grad_tol = 1e-7, restarts = 3)
set.seed(20261004)
n <- 400
d <- data.frame(x = rnorm(n, 0.7), z = rnorm(n),
                h = factor(sample(c("p", "q"), n, TRUE)))
u <- stats::rlogis(n) / exp(0.3 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -0.5) + (u > 0.4) + (u > 1.3) + (u > 2.3)
bform <- brms::bf(y | thres(gr = h) ~ x)
bfam <- brms::cumulative(threshold = "equidistant")
prior <- env$brms_flat_prior(bform, data = d, family = bfam)
print(prior)
setp <- function(prior, cls, grp, val) {
  i <- which(prior$class == cls & prior$group == grp & prior$coef == "")
  prior$prior[i] <- val; prior
}
arms <- list(
  both = list(fr = set_prior("normal(1, 0.1)", class = "delta", group = "p") +
                set_prior("normal(0, 0.5)", class = "Intercept", group = "q"),
              br = list(c("delta", "p", "normal(1, 0.1)"),
                        c("Intercept", "q", "normal(0, 0.5)"))),
  delta_only = list(fr = set_prior("normal(1, 0.1)", class = "delta", group = "p"),
                    br = list(c("delta", "p", "normal(1, 0.1)"))),
  intercept_only = list(fr = set_prior("normal(0, 0.5)", class = "Intercept",
                                       group = "q"),
                        br = list(c("Intercept", "q", "normal(0, 0.5)"))),
  intercept_p = list(fr = set_prior("normal(0, 0.5)", class = "Intercept",
                                    group = "p"),
                     br = list(c("Intercept", "p", "normal(0, 0.5)"))),
  none = list(fr = NULL, br = list()))
for (a in names(arms)) {
  cat("\n==", a, "==\n")
  p <- prior
  for (b in arms[[a]]$br) p <- setp(p, b[1], b[2], b[3])
  code <- brms::make_stancode(bform, data = d, family = bfam, prior = p)
  if (a == "both") {
    L <- strsplit(code, "\n")[[1]]
    cat(L[grep("^parameters", L):grep("^model", L)], sep = "\n")
    cat(grep("lprior \\+=", L, value = TRUE), sep = "\n")
    cat(grep("means_X|Xc", L, value = TRUE)[1:4], sep = "\n")
  }
  sdat <- env$brms_standata(bform, data = d, family = bfam, prior = p)
  sf <- suppressMessages(rstan::sampling(env$brms_stan_model(code),
                                         data = sdat, chains = 0))
  fit <- frm(y | thres(gr = h) ~ x, family = cumulative(threshold = "equidistant"),
             data = d, prior = arms[[a]]$fr, control = ctl)
  pars <- env$stan_pars_from_fit(fit, sdat, code)
  str(pars[setdiff(names(pars), character(0))])
  up <- rstan::unconstrain_pars(sf, pars)
  g <- rstan::grad_log_prob(sf, up, adjust_transform = TRUE)
  names(g) <- rstan:::unconstrained_param_names(sf)
  print(signif(g, 3))
  cat("frmtmb tau:", frmtmb:::ord_threshold_values(family(fit),
                                                   fit$estimates$tau_raw), "\n")
}
