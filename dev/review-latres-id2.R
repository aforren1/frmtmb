# Attack 3 and 10: the row-level identity with no translator, under cov = FALSE,
# at latents that are not modes: frm()'s joint env$f against brms's log_prob, at
# a feasible random zerr, per family. This closes the lane's NA stepped points
# (a point Stan rejects is resampled, not skipped). Also, on the post-link forms:
# is obj$fn a function of its argument, or of where the inner search started?
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
suppressPackageStartupMessages({library(rstan); library(brms)})
hel <- new.env(parent = asNamespace("frmtmb"))
hel$`%||%` <- `%||%`
sys.source("dev/equiv-lib.R", envir = hel)
set.seed(99)
d <- expand.grid(t = 1:4, g = factor(1:8))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
e <- rnorm(n, 0, .3)
d$cnt <- rpois(n, exp(1 + .3 * d$x + e))
d$bin <- rbinom(n, 1, plogis(.2 + .3 * d$x + e))
d$bi <- rbinom(n, 6, plogis(.2 + .3 * d$x + e))
d$tr <- 6L
d$pos <- rgamma(n, 4, 4 / exp(1 + .3 * d$x + e))
d$un <- pmin(pmax(rbeta(n, 4 * plogis(.3 * d$x + e), 4 * (1 - plogis(.3 * d$x + e))), 1e-3), 1 - 1e-3)
d$zi <- d$cnt * rbinom(n, 1, .8)
d$cn <- sample(c("none", "right"), n, TRUE, c(.8, .2))
d$hg <- d$pos * rbinom(n, 1, .8)
cases <- list(
  list("pois_log", cnt ~ x + ar(t, g), "poisson", list(), "poisson()"),
  list("pois_sqrt", cnt ~ x + ar(t, g), "poisson", list(), "poisson('sqrt')"),
  list("pois_ident", cnt ~ x + ar(t, g), "poisson", list(), "poisson('identity')"),
  list("pois_cens", cnt | cens(cn) ~ x + ar(t, g), "poisson", list(), "poisson()"),
  list("pois_trunc", cnt | trunc(ub = 500) ~ x + ar(t, g), "poisson", list(), "poisson()"),
  list("nb_log", cnt ~ x + ar(t, g), "negbinomial", list(shape = 3), "negbinomial()"),
  list("bern_logit", bin ~ x + ar(t, g), "bernoulli", list(), "bernoulli()"),
  list("bern_probit", bin ~ x + ar(t, g), "bernoulli", list(), "bernoulli('probit')"),
  list("bern_cloglog", bin ~ x + ar(t, g), "bernoulli", list(), "bernoulli('cloglog')"),
  list("binom_logit", bi | trials(tr) ~ x + ar(t, g), "binomial", list(), "binomial()"),
  list("binom_probit", bi | trials(tr) ~ x + ar(t, g), "binomial", list(), "binomial('probit')"),
  list("gamma_log", pos ~ x + ar(t, g), "Gamma", list(shape = 3), "Gamma(link = 'log')"),
  list("gamma_inv", pos ~ x + ar(t, g), "Gamma", list(shape = 3), "Gamma(link = 'inverse')"),
  list("beta_logit", un ~ x + ar(t, g), "beta", list(phi = 5), "Beta()"),
  list("lognormal", pos ~ x + ar(t, g), "lognormal", list(sigma = .5), "lognormal()"),
  list("weibull", pos ~ x + ar(t, g), "weibull", list(shape = 2), "weibull()"),
  list("exponential", pos ~ x + ar(t, g), "exponential", list(), "exponential()"),
  list("zip", zi ~ x + ar(t, g), "zero_inflated_poisson", list(zi = .2), "zero_inflated_poisson()"),
  list("hurdle_gamma", hg ~ x + ar(t, g), "hurdle_gamma", list(shape = 3, hu = .2), "hurdle_gamma()"),
  list("cumulative", NULL, "cumulative", list(), "cumulative()"),
  list("ar_T_probit", bin ~ x + ar(t, g, cov = TRUE), "bernoulli", list(), "bernoulli('probit')"))
d$ord <- factor(findInterval(.3 * d$x + e + rlogis(n), c(-1, 0, 1)) + 1L, ordered = TRUE)
cases[[20]][[2]] <- ord ~ x + ar(t, g)
out <- list()
for (cc in cases) {
  lab <- cc[[1]]
  form <- cc[[2]]
  fam_txt <- cc[[5]]
  fb <- eval(parse(text = fam_txt), asNamespace("brms"))
  ff <- eval(parse(text = fam_txt), asNamespace("frmtmb"))
  dp <- cc[[4]]
  pri <- c(set_prior("", class = "sderr"), set_prior("", class = "Intercept"))
  if (length(dp)) for (nm in names(dp)) if (nm != "zi" && nm != "hu") pri <- c(pri, set_prior("", class = nm))
  code <- brms::stancode(form, d, family = fb, prior = pri)
  sdat <- brms::standata(form, d, family = fb, prior = pri)
  m <- msg_of(hel$brms_stan_model(code))
  if (!is.null(m$err)) {
    cat(lab, "COMPILE ERR", substr(m$err, 1, 200), "\n")
    next
  }
  sf <- rstan::sampling(m$val, data = sdat, chains = 0, refresh = 0)
  fr <- msg_of(frm(form, data = d, family = ff, dry_run = "objective"))
  fit <- NULL
  mm <- msg_of(frm(form, data = d, family = ff))
  if (!is.null(mm$err)) {
    cat(lab, "frm ERR", substr(mm$err, 1, 200), "\n")
    next
  }
  fit <- mm$val
  env <- fit$obj$env
  bk <- fi("frame_latres")(fit$frame)[[1]]
  a <- bk$aux
  xm <- mean(d$x)
  ord <- identical(fb$family, "cumulative")
  res <- numeric(0)
  rej <- 0L
  for (pt in 1:6) {
    set.seed(1000 + pt)
    bc <- runif(1, -.5, .8)
    b0 <- switch(lab, pois_ident = 6, pois_sqrt = 3, gamma_inv = 2, bc + 1)
    b1 <- runif(1, -.2, .4)
    s <- runif(1, .1, .9)
    ar <- runif(1, -.8, .8)
    ok <- FALSE
    for (tries in 1:200) {
      z <- rnorm(n) * .7
      pl <- list(b = array(b1, 1), zerr = z, sderr = s, ar = array(ar, 1))
      if (ord) {
        pl$Intercept <- c(-1.2, 0, 1.1) * 1 + b1 * xm
      } else {
        pl$Intercept <- b0 + b1 * xm
      }
      for (nm in names(dp)) pl[[nm]] <- dp[[nm]]
      u <- tryCatch(rstan::unconstrain_pars(sf, pl), error = function(e) NULL)
      if (is.null(u)) stop(lab, " unconstrain failed")
      lp <- tryCatch(as.numeric(rstan::log_prob(sf, u, adjust_transform = FALSE)), error = function(e) NA_real_)
      if (is.finite(lp)) {
        ok <- TRUE
        break
      }
      rej <- rej + 1L
    }
    if (!ok) {
      cat(lab, "pt", pt, "no feasible point\n")
      next
    }
    full <- env$last.par.best
    nm <- names(full)
    full[nm == "beta"] <- if (ord) c(b1) else c(b0, b1)
    if (ord) {
      # the frame's internal thresholds
      th_idx <- which(nm == "beta" | nm == "thres" | grepl("thres", nm))
    }
    th <- c(log(s), if (a$kind == "condar") qlogis((ar + 1) / 2) else ar / sqrt(1 - ar^2))
    full[nm == "theta"] <- th
    full[nm == "b"] <- z
    if (length(dp)) {
      for (dn in names(dp)) {
        k <- which(nm == "betad")
        vals <- vapply(names(dp), function(q) switch(q, zi = , hu = qlogis(dp[[q]]), log(dp[[q]])), 0)
        full[k] <- vals
      }
    }
    if (ord) {
      res <- c(res, NA_real_)
      next
    }
    v <- -env$f(full)
    res <- c(res, v - lp)
  }
  cat(sprintf("%-14s on_response=%-5s kind=%-6s rejected=%3d  diff(frm - brms): %s\n", lab, isTRUE(a$on_response), a$kind, rej, paste(sprintf("%.2e", res), collapse = " ")))
  out[[lab]] <- res
}
