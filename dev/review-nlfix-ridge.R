# Item 1: frm() against brms's own optimum and posterior, on the ridge specs and
# on the default-start saddle. rstan::optimizing from 6 starts under flat priors;
# one short posterior for the identified saddle spec. Seeds as named.
source("dev/review-nlfix-common.R")
suppressPackageStartupMessages(library(testthat))
E <- new.env(parent = asNamespace("frmtmb"))
E$`%||%` <- function(a, b) if (is.null(a)) b else a
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) sys.source(f, envir = E)
opt <- function(fb, fam, d, k = 6) {
  prior <- E$brms_flat_prior(fb, data = d, family = fam)
  code <- brms::make_stancode(fb, data = d, family = fam, prior = prior)
  sdat <- brms::standata(fb, data = d, family = fam, prior = prior)
  m <- E$brms_stan_model(code)
  vals <- numeric(0)
  for (i in seq_len(k)) {
    r <- tryCatch(rstan::optimizing(m, data = sdat, seed = 100 + i, iter = 5000, hessian = FALSE, as_vector = FALSE, verbose = FALSE),
                  error = function(e) NULL)
    if (!is.null(r)) vals <- c(vals, r$value)
  }
  list(vals = vals, m = m, sdat = sdat)
}
report <- function(lab, ff, fb, fam, d, ...) {
  cat("\n===", lab, "\n")
  r <- tt(frm(ff, family = eval(parse(text = fam)), data = d, ...))
  if (iserr(r)) { cat("frm ERROR", r$msg, "\n"); return(invisible()) }
  cat(sprintf("frm ll=%.6f conv=%s msg=%s\n", as.numeric(logLik(r)), r$opt$convergence, r$opt$message))
  print(round(unlist(r$estimates[c("beta", "tau_raw")]), 5))
  cat(" warnings:\n"); for (w in attr(r, "w")) cat("   ", substr(w, 1, 700), "\n")
  g <- tt(max(abs(r$obj$gr(r$opt$par)))); cat(" max|grad|:", if (iserr(g)) g$msg else g, "\n")
  bf <- eval(parse(text = fam), envir = asNamespace("brms"))
  o <- opt(fb, bf, d)
  cat(" brms optimizing, 6 starts:", paste(format(sort(o$vals, decreasing = TRUE), digits = 9), collapse = " "), "\n")
  invisible(list(r = r, o = o))
}
for (k in 1:2) {
  fam <- c("cumulative()", "sratio()")[k]
  d <- E$eq_dat("ordinal", seed = 4930 + k)
  report(paste("ord nonlin spec", fam, "lane's start"), bf(y ~ b1 * exp(b2 * x), b1 ~ 1, b2 ~ 1 + z, nl = TRUE),
         brms::bf(y ~ b1 * exp(b2 * x), b1 ~ 1, b2 ~ 1 + z, nl = TRUE), fam, d,
         start = list(beta = c("b1_(Intercept)" = 0.5, "b2_(Intercept)" = 0.5)))
  report(paste("ord nonlin spec", fam, "default start"), bf(y ~ b1 * exp(b2 * x), b1 ~ 1, b2 ~ 1 + z, nl = TRUE),
         brms::bf(y ~ b1 * exp(b2 * x), b1 ~ 1, b2 ~ 1 + z, nl = TRUE), fam, d)
}
# an identified curve: the linear limit is not the maximum
set.seed(9202)
n <- 300
d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
d$y <- ordered(cut(1.5 * exp(0.6 * d$x) + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
x <- report("identified curve, cumulative, default start", bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE),
       brms::bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), "cumulative()", d)
report("identified curve, cumulative, start b1 = 1", bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE),
       brms::bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), "cumulative()", d, start = list(beta = c("b1_(Intercept)" = 1)))
report("identified curve, cumulative, start b2 = 0.1", bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE),
       brms::bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), "cumulative()", d, start = list(beta = c("b2_(Intercept)" = 0.1)))
cat("\nsummary at the default start\n")
print(summary(x$r))
cat("thresholds-only fit logLik:", as.numeric(logLik(frm(y ~ 1, family = cumulative(), data = d))), "\n")
cat("\nbrms posterior on the same spec, default inits, flat priors on b, 2 chains, 500 iter\n")
s <- suppressWarnings(rstan::sampling(x$o$m, data = x$o$sdat, chains = 2, iter = 500, seed = 11, refresh = 0))
ss <- rstan::summary(s, pars = c("b_b1", "b_b2", "lp__"))$summary
print(round(ss[, c("mean", "sd", "n_eff", "Rhat")], 3))
for (fm in c("sratio", "cratio", "acat", "hurdle_cumulative")) {
  dd <- d
  if (fm == "hurdle_cumulative") { set.seed(1); dd$y <- ordered(ifelse(runif(n) < 0.2, 0L, as.integer(d$y)), levels = 0:4) }
  r <- tt(frm(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), family = get(fm)(), data = dd))
  r2 <- tt(frm(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), family = get(fm)(), data = dd, start = list(beta = c("b1_(Intercept)" = 1, "b2_(Intercept)" = 0.3))))
  cat(fm, "default:", if (iserr(r)) r$msg else as.numeric(logLik(r)), " started:", if (iserr(r2)) r2$msg else as.numeric(logLik(r2)), " warn:", substr(paste(attr(r, "w"), collapse = "|"), 1, 110), "\n")
}
cat("\nprobit, equidistant, RE, other curves under cumulative from the default start\n")
for (sp in list(list(f = bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), fam = cumulative("probit")),
                list(f = bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), fam = cumulative(threshold = "equidistant")),
                list(f = bf(y ~ b1 * x + b2 * z, b1 + b2 ~ 1, nl = TRUE), fam = cumulative()),
                list(f = bf(y ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE), fam = cumulative()),
                list(f = bf(y ~ b1 / (1 + exp(b2 * x)), b1 + b2 ~ 1, nl = TRUE), fam = cumulative()),
                list(f = bf(y ~ inv_logit(b1 + b2 * x), b1 + b2 ~ 1, nl = TRUE), fam = cumulative()))) {
  r <- tt(frm(sp$f, family = sp$fam, data = d))
  cat(deparse1(sp$f$formula), ":", if (iserr(r)) substr(r$msg, 1, 120) else paste(format(as.numeric(logLik(r)), digits = 8), "beta", paste(format(r$estimates$beta, digits = 3), collapse = ",")), " warn:", length(attr(r, "w")), "\n")
}
