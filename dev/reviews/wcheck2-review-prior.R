# Review of lane wcheck2, item 4985 and the cost of the new rules.
#   Rscript dev/reviews/wcheck2-review-prior.R <base|lane> <part>
# part conds: refused models under valid and invalid priors; every
#   condition raised (class and first words) and the error.
# part ident: models the check lets through; logLik and parameters in
#   hex, to compare across the two arms.
# part time: check_weights_unbounded() alone on large frames.
a <- commandArgs(TRUE)
lib_arg <- a[1]
part <- a[2]
source("dev/wcheck2-common.R")
ns <- asNamespace("frmtmb")

conds <- function(lab, expr) {
  seen <- character()
  r <- tryCatch(withCallingHandlers(expr,
    warning = function(w) {
      seen <<- c(seen, paste0("W ", class(w)[1], ": ",
                              substr(gsub("\\s+", " ", conditionMessage(w)),
                                     1, 50)))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      seen <<- c(seen, paste0("M ", class(m)[1], ": ",
                              substr(gsub("\\s+", " ", conditionMessage(m)),
                                     1, 50)))
      invokeRestart("muffleMessage")
    }), error = function(e) e)
  res <- if (inherits(r, "error")) {
    m <- gsub("\\s+", " ", conditionMessage(r))
    paste(if (grepl("has no maximum|is improper", m)) "REFUSED" else "ERROR",
          substr(m, 1, 110))
  } else sprintf("fit logLik %.6g conv %d", as.numeric(logLik(r)),
                 r$opt$convergence)
  cat(sprintf("%-46s %s\n%-46s   conditions: %s\n", lab, res, "",
              if (length(seen)) paste(seen, collapse = " || ") else "none"))
}

set.seed(20261111)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n),
                g = factor(rep(1:16, each = 5)),
                fac = factor(rep(c("a", "b"), 40)))
d$y <- 0.5 + 0.4 * d$x + rnorm(n)
d$y2 <- d$y + rnorm(n)
d$w <- 1.1
d$yz <- ifelse(seq_len(n) %% 5 == 0, 0, exp(d$y / 2))
fo <- bf(y | weights(w) ~ x + (1 | obs))
fs <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + x)
S <- set_prior

if (part == "conds") {
  conds("no prior", frm(fo, data = d))
  conds("normal(0,5) b", frm(fo, data = d, prior = S("normal(0, 5)",
                                                       class = "b")))
  conds("normal(0,5) Intercept", frm(fo, data = d,
    prior = S("normal(0, 5)", class = "Intercept")))
  conds("student_t sd", frm(fo, data = d,
    prior = S("student_t(3, 0, 2.5)", class = "sd")))
  conds("normal sd lb = 0 (restated)", frm(fo, data = d,
    prior = S("normal(0, 1)", class = "sd", lb = 0)))
  conds("normal sd group obs", frm(fo, data = d,
    prior = S("normal(0, 1)", class = "sd", group = "obs")))
  conds("exponential(1) sigma", frm(fo, data = d,
    prior = S("exponential(1)", class = "sigma")))
  conds("normal b ub = 1 (bound warning)", frm(fo, data = d,
    prior = S("normal(0, 1)", class = "b", ub = 1)))
  conds("lognormal b (lower-bounded density)", frm(fo, data = d,
    prior = S("lognormal(0, 1)", class = "b")))
  if (requireNamespace("brms", quietly = TRUE)) {
    conds("brms::prior() object", frm(fo, data = d,
      prior = brms::prior(normal(0, 5), class = b)))
  }
  conds("quadrature, normal(0,5) b", frm(fo, data = d, quadrature = TRUE,
    prior = S("normal(0, 5)", class = "b")))
  # invalid priors
  conds("INVALID sd group zz", frm(fo, data = d,
    prior = S("normal(0, 1)", class = "sd", group = "zz")))
  conds("INVALID b coef zz", frm(fo, data = d,
    prior = S("normal(0, 1)", class = "b", coef = "zz")))
  conds("INVALID b coef Intercept, sigma ~ 1 + x", frm(fs, data = d,
    prior = S("normal(0, 1)", class = "b", coef = "Intercept",
              dpar = "sigma")))
  conds("INVALID class cor (no correlated block)", frm(fo, data = d,
    prior = S("lkj(2)", class = "cor")))
  conds("INVALID dpar zi on gaussian", frm(fo, data = d,
    prior = S("normal(0, 1)", class = "Intercept", dpar = "zi")))
  # a multivariate model refused on y, prior on y2 only
  fm <- mvbf(bf(y | weights(w) ~ x + (1 | obs)), bf(y2 ~ x))
  conds("mv, normal(0,1) sigma of y2", frm(fm, data = d,
    prior = S("normal(0, 1)", class = "sigma", resp = "y2")))
  conds("mv, INVALID resp zz", frm(fm, data = d,
    prior = S("normal(0, 1)", class = "sigma", resp = "zz")))
  # a hurdle family, now refused, with a prior on hu
  fh <- bf(yz | weights(w) ~ x + (1 | obs))
  d$w <- ifelse(d$yz == 0, 0.5, 1.2)
  conds("hurdle_lognormal, beta(1,1) hu", frm(fh, data = d,
    family = hurdle_lognormal(), prior = S("beta(1, 1)", class = "hu")))
}

if (part == "ident") {
  hx <- function(fit) paste(sprintf("%a", c(as.numeric(logLik(fit)),
                                            fit$opt$par)), collapse = " ")
  run <- function(lab, expr) {
    t0 <- proc.time()[["elapsed"]]
    f <- q(expr)
    cat(sprintf("%-40s %.2fs %s\n", lab, proc.time()[["elapsed"]] - t0, hx(f)))
  }
  d$w <- 0.9
  run("gaussian w 0.9 (1|obs)", frm(fo, data = d))
  d$w <- ifelse(d$fac == "a", 0.8, 1.05)
  run("sigma ~ fac, normal Intercept (R1)", frm(
    bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + fac), data = d,
    prior = S("normal(0, 1)", class = "Intercept", dpar = "sigma")))
  d$w <- 1.3
  run("gaussian w 1.3 (1|g)", frm(bf(y | weights(w) ~ x + (1 | g)),
                                  data = d))
  d$w <- ifelse(d$yz == 0, 3, 0.95)
  run("hurdle_lognormal E<0 (1|obs)", frm(bf(yz | weights(w) ~ x + (1 | obs)),
                                         data = d, family = hurdle_lognormal()))
  run("lognormal unweighted (1|g)", frm(bf(yz + 0.1 ~ x + (1 | g)), data = d,
                                        family = lognormal()))
}

if (part == "time") {
  for (cfg in list(c(2000, 100), c(6000, 300), c(10000, 1000))) {
    nn <- cfg[1]; L <- cfg[2]
    set.seed(3)
    dd <- data.frame(obs = factor(seq_len(nn)), x = rnorm(nn),
                     fac = factor(sample(sprintf("l%04d", seq_len(L)), nn,
                                         replace = TRUE)))
    dd$y <- rnorm(nn)
    dd$w <- 0.95
    fl <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + fac)
    fr <- q(frm(fl, data = dd, dry_run = "frame"))
    for (pp in list(NULL, S("normal(0, 1)", class = "Intercept",
                            dpar = "sigma"))) {
      t0 <- proc.time()[["elapsed"]]
      r <- tryCatch(ns$check_weights_unbounded(fr$spec, fr, pp),
                    error = function(e) "REFUSED")
      cat(sprintf("n %d, %d levels, prior %s: check %.2fs, %s\n", nn, L,
                  if (is.null(pp)) "none" else "Intercept",
                  proc.time()[["elapsed"]] - t0,
                  if (is.null(r)) "pass" else r))
    }
  }
}
cat("DONE\n")
