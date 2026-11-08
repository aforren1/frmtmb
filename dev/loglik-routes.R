# The first question of the lane brief: on each route that adds a prior,
# is "objective minus log prior at the estimates" the data-only objective
# at the same estimates?
#
# For each route, a MAP fit `fp` and the same model with no prior `f0`
# share one parameter layout, so `f0$obj$fn(fp$opt$par)` is the data-only
# objective of that route (same random set, same integrate=, same profile)
# evaluated at the MAP estimates. The split value is
#   split = -fp$opt$objective + nlp(fp$estimates)
# with nlp the fit's own neg_log_prior_fn() over its resolved entries.
# importance= is handled separately: its tape is anchored at the fit's
# own estimate, so the data-only value is the fit's own importance
# closure with the prior removed (see imp_data_only()).
#
# Arg 1: "before" (base build rellib-r18, internals only) or "after".
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "after"
libs <- c("C:/Users/adf44/source/r/wt-loglik-lib",
          "C:/Users/adf44/source/r/rellib-r18",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "before") libs <- libs[-1]
.libPaths(libs)
suppressMessages(library(frmtmb))
cat("arm:", arm, " lib:", dirname(find.package("frmtmb")),
    " version:", format(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb")

nlp_at <- function(fit) {
  ri <- ns$resolve_prior_input(list(frame = fit$frame, spec = fit$spec),
                               fit$prior)
  comps <- unique(unlist(lapply(ri$entries, function(e)
    c(e$comp, e$offset$comp))))
  list(nlp = ns$neg_log_prior_fn(ri$entries)(fit$estimates),
       comps = comps, n = length(ri$entries))
}

rows <- list()
report <- function(route, fp, f0, data_val = NULL) {
  np <- nlp_at(fp)
  split <- -fp$opt$objective + np$nlp
  if (is.null(data_val)) data_val <- -f0$obj$fn(fp$opt$par)
  rel <- abs(split - data_val) / max(1, abs(data_val))
  after <- if (arm == "after") {
    tryCatch(format(as.numeric(logLik(fp)), digits = 15),
             error = function(e) paste("REFUSED:",
                                       substr(conditionMessage(e), 1, 60)))
  } else "-"
  cat(sprintf(paste0("%-26s pen %.10f  log_prior %.10f  split %.10f  ",
                     "data %.10f  rel %.2e  comps [%s]  logLik() %s\n"),
              route, -fp$opt$objective, -np$nlp, split, data_val, rel,
              paste(np$comps, collapse = ","), after))
  rows[[route]] <<- data.frame(route = route, pen = -fp$opt$objective,
                               log_prior = -np$nlp, split = split,
                               data = data_val, rel = rel)
}

set.seed(4101)
n <- 200; ng <- 20
d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(ng), each = n / ng)))
eta <- 1 + 0.5 * d$x + rnorm(ng, 0, 0.7)[d$g]
d$y <- eta + rnorm(n)
d$yc <- rpois(n, exp(0.2 + 0.3 * d$x + rnorm(ng, 0, 0.5)[d$g]))
d$yb <- rbinom(n, 1, plogis(0.3 * d$x + rnorm(ng, 0, 0.8)[d$g]))
pb <- set_prior("normal(0, 0.2)", class = "b") +
  set_prior("student_t(3, 0, 1)", class = "sd") +
  set_prior("normal(2, 0.5)", class = "Intercept")
pb_sd <- set_prior("normal(0, 0.3)", class = "sd")

# A. Laplace ML, priors on b, Intercept (centered offset) and sd
f0 <- frm(y ~ x + (1 | g), data = d)
fp <- frm(y ~ x + (1 | g), data = d, prior = pb)
report("A laplace gaussian", fp, f0)
f0 <- frm(yc ~ x + (1 | g), data = d, family = poisson())
fp <- frm(yc ~ x + (1 | g), data = d, family = poisson(), prior = pb)
report("A laplace poisson", fp, f0)

# B. profile = TRUE (beta profiled inside the inner problem)
ctl <- frmtmb_control(profile = TRUE)
f0 <- frm(y ~ x + (1 | g), data = d, control = ctl)
fp <- frm(y ~ x + (1 | g), data = d, control = ctl, prior = pb)
# f0's tape re-maximizes beta (profiled), which is not "at the fit's
# estimates"; the non-profiled data tape is evaluated at (beta, theta)
fA <- frm(y ~ x + (1 | g), data = d)
pv <- unlist(fp$estimates[unique(names(fA$obj$par))])
report("B profile (tape re-max)", fp, f0)
report("B profile", fp, NULL, data_val = -fA$obj$fn(unname(pv)))

# C. quadrature = TRUE (adaptive Gauss-Kronrod per scalar effect)
f0 <- frm(yb ~ x + (1 | g), data = d, family = bernoulli(),
          quadrature = TRUE)
fp <- frm(yb ~ x + (1 | g), data = d, family = bernoulli(),
          quadrature = TRUE, prior = pb)
report("C quadrature", fp, f0)

# D. importance =: the fit's own corrected objective, prior removed. The
# closure the tape was built from is found on the fit's tape environment
# and evaluated with neg_log_prior = NULL at the same plan and offsets.
fp <- frm(yc ~ x + (1 | g), data = d, family = poisson(), prior = pb,
          importance = 200)
imp_data_only <- function(fit) {
  is_imp_fn <- function(v) {
    is.function(v) && identical(names(formals(v)), "pars") &&
      !is.null(environment(v)) &&
      exists("offs", environment(v), inherits = FALSE) &&
      exists("neg_log_prior", environment(v), inherits = FALSE) &&
      any(grepl("neg_log_prior", deparse(body(v)), fixed = TRUE))
  }
  # breadth-first over environments reachable from the tape's
  seen <- list()
  todo <- list(fit$obj$env)
  fn <- NULL
  while (length(todo) && is.null(fn) && length(seen) < 500) {
    e <- todo[[1L]]
    todo <- todo[-1L]
    if (any(vapply(seen, identical, TRUE, e))) next
    seen[[length(seen) + 1L]] <- e
    for (nm in ls(e, all.names = TRUE)) {
      v <- tryCatch(get(nm, e), error = function(err) NULL)
      if (is_imp_fn(v)) {
        fn <- v
        break
      }
      if (is.function(v) && is.environment(environment(v)) &&
            !isNamespace(environment(v)) &&
            !identical(environment(v), globalenv())) {
        todo[[length(todo) + 1L]] <- environment(v)
      }
      if (is.environment(v) && !isNamespace(v)) {
        todo[[length(todo) + 1L]] <- v
      }
    }
  }
  if (is.null(fn)) stop("importance closure not found on the tape env")
  e <- environment(fn)
  old <- e$neg_log_prior
  on.exit(assign("neg_log_prior", old, envir = e))
  assign("neg_log_prior", NULL, envir = e)
  pars <- fit$estimates
  pars[["b"]] <- NULL
  pars <- pars[names(fit$obj$env$parameters)]
  -fn(pars)
}
dv <- tryCatch(imp_data_only(fp), error = function(e) {
  cat("importance closure:", conditionMessage(e), "\n")
  NA_real_
})
report("D importance", fp, NULL, data_val = dv)

# E. REML: a prior on sd only reads outer parameters; a prior on b reads
# beta, which REML integrates
f0 <- frm(y ~ x + (1 | g), data = d, REML = TRUE)
fp <- frm(y ~ x + (1 | g), data = d, REML = TRUE, prior = pb_sd)
report("E REML, sd prior", fp, f0)
fp <- frm(y ~ x + (1 | g), data = d, REML = TRUE, prior = pb)
report("E REML, b prior", fp, f0)

# F. exact gp(): non-centered, z inside the Laplace integral
set.seed(4102)
dg <- data.frame(x = seq(0, 1, length.out = 40))
dg$y <- sin(4 * dg$x) + rnorm(40, 0, 0.3)
pg <- set_prior("normal(0, 1)", class = "sdgp") +
  set_prior("inv_gamma(3, 0.5)", class = "lscale")
f0 <- frm(y ~ gp(x), data = dg)
fp <- frm(y ~ gp(x), data = dg, prior = pg)
report("F exact gp", fp, f0)

# G. latent residual: ar() on a family with no residual
set.seed(4103)
nt <- 80
dt <- data.frame(t = seq_len(nt))
dt$yc <- rpois(nt, exp(0.5 + as.numeric(arima.sim(list(ar = 0.6), nt,
                                                   sd = 0.4))))
pa <- set_prior("normal(0, 0.5)", class = "ar") +
  set_prior("normal(0, 0.5)", class = "sderr")
f0 <- frm(yc ~ 1 + ar(time = t), data = dt, family = poisson())
fp <- frm(yc ~ 1 + ar(time = t), data = dt, family = poisson(), prior = pa)
report("G latent residual", fp, f0)

# H. me(): latent values in "miss"
set.seed(4104)
dm <- data.frame(xt = rnorm(100))
dm$sdx <- 0.3
dm$xo <- dm$xt + rnorm(100, 0, dm$sdx)
dm$y <- 1 + 0.7 * dm$xt + rnorm(100, 0, 0.5)
pm <- set_prior("normal(0, 0.3)", class = "b") +
  set_prior("normal(0, 1)", class = "sdme")
f0 <- frm(y ~ me(xo, sdx), data = dm)
fp <- frm(y ~ me(xo, sdx), data = dm, prior = pm)
report("H me()", fp, f0)

# I. mo() with a simo Dirichlet prior (the chart Jacobian)
set.seed(4105)
dmo <- data.frame(xo = sample(1:4, 120, TRUE))
dmo$y <- c(0, 0.2, 1.4, 1.6)[dmo$xo] + rnorm(120, 0, 0.5)
po <- set_prior("dirichlet(c(2, 2, 2))", class = "simo", coef = "moxo1")
f0 <- frm(y ~ mo(xo), data = dmo)
fp <- frm(y ~ mo(xo), data = dmo, prior = po)
report("I mo() simo", fp, f0)

# J. nonlinear with nlpar priors and a group-level term
set.seed(4106)
dn <- data.frame(x = runif(150, 0, 5), g = factor(rep(1:15, each = 10)))
dn$y <- 3 * (1 - exp(-0.8 * dn$x)) + rnorm(15, 0, 0.3)[dn$g] +
  rnorm(150, 0, 0.2)
nlb <- bf(y ~ a * (1 - exp(-k * x)), a ~ 1 + (1 | g), k ~ 1, nl = TRUE)
pn <- set_prior("normal(3, 1)", nlpar = "a") +
  set_prior("normal(1, 0.5)", nlpar = "k")
f0 <- frm(nlb, data = dn, start = list(beta = c(3, 1)))
fp <- frm(nlb, data = dn, prior = pn)
report("J nonlinear", fp, f0)

tab <- do.call(rbind, rows)
cat("\nmax rel over separable routes:",
    format(max(tab$rel[!tab$route %in% c("E REML, b prior", "B profile (tape re-max)")], na.rm = TRUE),
           digits = 3), "\n")
