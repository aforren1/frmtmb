# Reviewer, lane ordinal: (a) MAP with class delta and class Intercept
# priors against brms's compiled program: brms's gradient of the log
# posterior on the unconstrained scale (Jacobian included, which is
# frmtmb's "sd" placement) at frmtmb's MAP; (b) sum_to_zero: brms's own
# optimum, found by BFGS on its program from a perturbed start, against
# frmtmb's ML fit: thresholds, coefficients, logLik and fitted.
# Seeds: data 20261004, start perturbation 11.
# Output: dev/ordinal-rev-log-map.txt
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
cat("table(y):", table(d$y), "\n")

brms_setup <- function(bform, family, prior_mod) {
  prior <- env$brms_flat_prior(bform, data = d, family = family)
  for (pm in prior_mod) {
    i <- which(prior$class == pm$class & prior$coef == "" &
                 prior$group == (pm$group %||% "") & prior$dpar == "")
    stopifnot(length(i) == 1L)
    prior$prior[i] <- pm$prior
  }
  code <- brms::make_stancode(bform, data = d, family = family, prior = prior)
  sdat <- env$brms_standata(bform, data = d, family = family, prior = prior)
  mod <- env$brms_stan_model(code)
  sf <- suppressMessages(rstan::sampling(mod, data = sdat, chains = 0))
  list(code = code, sdat = sdat, sf = sf, prior = prior)
}

cat("\n#### (a) MAP with priors: brms gradient (Jacobian on) at frmtmb's MAP\n")
map_rows <- list(
  list(name = "cumulative equidistant, delta ~ N(1.5, .2), Intercept ~ N(-1, .5)",
       f = y ~ x, fam = cumulative(threshold = "equidistant"),
       bfam = brms::cumulative(threshold = "equidistant"),
       pl = set_prior("normal(1.5, 0.2)", class = "delta") +
         set_prior("normal(-1, 0.5)", class = "Intercept"),
       bp = list(list(class = "delta", prior = "normal(1.5, 0.2)"),
                 list(class = "Intercept", prior = "normal(-1, 0.5)"))),
  list(name = "sratio equidistant, delta ~ N(0.5, .1), Intercept ~ N(0, 2)",
       f = y ~ x, fam = sratio(threshold = "equidistant"),
       bfam = brms::sratio(threshold = "equidistant"),
       pl = set_prior("normal(0.5, 0.1)", class = "delta") +
         set_prior("normal(0, 2)", class = "Intercept"),
       bp = list(list(class = "delta", prior = "normal(0.5, 0.1)"),
                 list(class = "Intercept", prior = "normal(0, 2)"))),
  list(name = "cumulative equidistant thres(gr = h), delta|p ~ N(1, .1), Intercept|q ~ N(0, .5)",
       f = y | thres(gr = h) ~ x, fam = cumulative(threshold = "equidistant"),
       bfam = brms::cumulative(threshold = "equidistant"),
       pl = set_prior("normal(1, 0.1)", class = "delta", group = "p") +
         set_prior("normal(0, 0.5)", class = "Intercept", group = "q"),
       bp = list(list(class = "delta", group = "p", prior = "normal(1, 0.1)"),
                 list(class = "Intercept", group = "q",
                      prior = "normal(0, 0.5)"))),
  list(name = "acat equidistant, delta ~ student_t(3, 0.3, 0.1), disc ~ 0 + z",
       f = bf(y ~ x, disc ~ 0 + z), fam = acat(threshold = "equidistant"),
       bfam = brms::acat(threshold = "equidistant"),
       pl = set_prior("student_t(3, 0.3, 0.1)", class = "delta"),
       bp = list(list(class = "delta", prior = "student_t(3, 0.3, 0.1)")))
)
for (r in map_rows) {
  cat("\n==", r$name, "==\n")
  tryCatch({
    fit <- frm(r$f, family = r$fam, data = d, prior = r$pl, control = ctl)
    bform <- if (inherits(r$f, "formula")) brms::bf(r$f) else
      do.call(brms::bf, list(y ~ x, disc ~ 0 + z))
    b <- brms_setup(bform, r$bfam, r$bp)
    pars <- env$stan_pars_from_fit(fit, b$sdat, b$code)
    up <- rstan::unconstrain_pars(b$sf, pars)
    g <- rstan::grad_log_prob(b$sf, up, adjust_transform = TRUE)
    g0 <- rstan::grad_log_prob(b$sf, up, adjust_transform = FALSE)
    cat(sprintf("brms max|grad| Jacobian on %.3g, off %.3g; frmtmb max|grad| %.3g\n",
                max(abs(g)), max(abs(g0)),
                max(abs(fit$obj$gr(fit$opt$par)))))
    cat("grad (Jacobian on):", format(signif(g, 3)), "\n")
  }, error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
}

cat("\n#### (b) sum_to_zero: brms's own optimum against frmtmb's\n")
stz_rows <- list(
  list(name = "cumulative probit sum_to_zero", f = y ~ x,
       fam = cumulative("probit", threshold = "sum_to_zero"),
       bfam = brms::cumulative("probit", threshold = "sum_to_zero"),
       dens = "cumulative", link = "probit"),
  list(name = "acat sum_to_zero", f = y ~ x,
       fam = acat(threshold = "sum_to_zero"),
       bfam = brms::acat(threshold = "sum_to_zero"),
       dens = "acat", link = "logit"),
  list(name = "sratio sum_to_zero, disc ~ 0 + z", f = bf(y ~ x, disc ~ 0 + z),
       bf = brms::bf(y ~ x, disc ~ 0 + z),
       fam = sratio(threshold = "sum_to_zero"),
       bfam = brms::sratio(threshold = "sum_to_zero"),
       dens = "sratio", link = "logit")
)
for (r in stz_rows) {
  cat("\n==", r$name, "==\n")
  tryCatch({
    fit <- frm(r$f, family = r$fam, data = d, control = ctl)
    bform <- r[["bf"]] %||% brms::bf(r$f)
    b <- brms_setup(bform, r$bfam, list())
    pars <- env$stan_pars_from_fit(fit, b$sdat, b$code)
    up0 <- rstan::unconstrain_pars(b$sf, pars)
    set.seed(11)
    st <- up0 + rnorm(length(up0), 0, 0.3)
    fn <- function(p) -rstan::log_prob(b$sf, p, adjust_transform = FALSE)
    gr <- function(p) -rstan::grad_log_prob(b$sf, p, adjust_transform = FALSE)
    o <- optim(st, fn, gr, method = "BFGS",
               control = list(maxit = 5000, reltol = 1e-15))
    cp <- rstan::constrain_pars(b$sf, o$par)
    tb <- as.numeric(cp$b_Intercept)
    tf <- frmtmb:::ord_threshold_values(family(fit), fit$estimates$tau_raw)
    bb <- as.numeric(cp$b)
    bf_ <- unname(fixef(fit, flatten = TRUE)["x"])
    cat(sprintf("optim convergence %d; brms max|grad| at its optimum %.3g\n",
                o$convergence, max(abs(gr(o$par)))))
    cat("brms Intercept (uncentered, common location unplaced):",
        format(signif(as.numeric(cp$Intercept), 6)), "\n")
    cat("brms b_Intercept:", format(signif(tb, 10)), "\n")
    cat("frmtmb thresholds:", format(signif(tf, 10)), "\n")
    cat(sprintf("max rel diff thresholds %.3g; b: brms %.10g frmtmb %.10g rel %.3g\n",
                max(abs(tb - tf) / abs(tf)), bb[1], bf_, abs(bb[1] - bf_) / abs(bf_)))
    cat(sprintf("log-lik: brms optimum %.10g, frmtmb %.10g, rel %.3g\n",
                -o$value, as.numeric(logLik(fit)),
                abs(-o$value - as.numeric(logLik(fit))) /
                  abs(as.numeric(logLik(fit)))))
    # fitted at brms's optimum, through brms's own R-side density
    eta <- as.numeric(d$x * bb[1])
    disc <- if (!is.null(cp$b_disc)) exp(d$z * as.numeric(cp$b_disc)) else 1
    dens <- get(paste0("d", r$dens), asNamespace("brms"))
    P <- dens(1:5, eta = eta, thres = matrix(tb, n, 4, byrow = TRUE),
              disc = disc, link = r$link)
    Pf <- fitted(fit)[, "Estimate", ]
    cat(sprintf("fitted: max |brms-at-brms-optimum - frmtmb| %.3g\n",
                max(abs(P - Pf))))
  }, error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
}
