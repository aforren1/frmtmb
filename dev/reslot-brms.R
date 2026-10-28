# Lane reslot: what brms 2.23.0 does, run rather than recalled.
# default_prior(), validate_prior() and stancode() only; no sampling.
.libPaths(c("C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
source_data <- function() {
  set.seed(1)
  G <- 8L; n <- 10L
  g <- factor(rep(seq_len(G), each = n))
  x <- rnorm(G * n); z <- rnorm(G * n)
  y <- 1 + 0.5 * x + rnorm(G)[g] + rnorm(G * n)
  y2 <- 0.5 * y + rnorm(G * n)
  t <- rep(seq_len(n), G)
  data.frame(y, y2, x, z, g, t)
}
dd <- source_data()
show <- function(label, expr) {
  cat("\n--", label, "\n")
  r <- tryCatch(expr, error = function(e) e)
  if (inherits(r, "error")) cat("  ERROR:", conditionMessage(r), "\n") else
    print(r)
  invisible(r)
}
lines_with <- function(code, pat) {
  l <- strsplit(as.character(code), "\n")[[1]]
  cat(paste0("  | ", grep(pat, l, value = TRUE)), sep = "\n")
}
f_re <- bf(y ~ x + z + (1 + x + z | g))

show("default_prior (1 + x + z | g)", default_prior(f_re, data = dd))
show("validate_prior class theta coef theta_4",
     validate_prior(set_prior("normal(0, 1)", class = "theta",
                              coef = "theta_4"), f_re, data = dd))
show("validate_prior class theta, no coef",
     validate_prior(set_prior("normal(0, 1)", class = "theta"), f_re,
                    data = dd))
show("validate_prior cor with coef",
     validate_prior(set_prior("lkj(2)", class = "cor", coef = "x"), f_re,
                    data = dd))
show("set_prior cor with lb",
     validate_prior(set_prior("lkj(2)", class = "cor", lb = 0), f_re,
                    data = dd))
show("set_prior bounds-only cor", validate_prior(
  set_prior("", class = "cor", group = "g", lb = -0.5, ub = 0.5), f_re,
  data = dd))
show("set_prior rescor with lb", validate_prior(
  set_prior("lkj(2)", class = "rescor", lb = 0),
  bf(y ~ x) + bf(y2 ~ x) + set_rescor(TRUE), data = dd))

cat("\n== flat on brms's parameters: every default replaced by \"\"\n")
flat_re <- c(set_prior("", class = "Intercept"),
             set_prior("", class = "sd"),
             set_prior("", class = "sigma"),
             set_prior("", class = "cor"))
show("validate_prior flat RE", validate_prior(flat_re, f_re, data = dd))
lines_with(stancode(f_re, data = dd, prior = flat_re),
           "lprior|target|lkj|cholesky|<lower|sd_1|sigma;")
show("rescor flat", validate_prior(
  set_prior("", class = "rescor"), bf(y ~ x) + bf(y2 ~ x) +
    set_rescor(TRUE), data = dd))
lines_with(stancode(bf(y ~ x) + bf(y2 ~ x) + set_rescor(TRUE), data = dd,
                    prior = set_prior("", class = "rescor")),
           "lprior|Lrescor")

cat("\n== ar/ma tables\n")
show("default_prior ar(cov = TRUE)",
     default_prior(bf(y ~ x + ar(time = t, gr = g, cov = TRUE)), data = dd))
show("default_prior ma(cov = TRUE)",
     default_prior(bf(y ~ x + ma(time = t, gr = g, cov = TRUE)), data = dd))
show("default_prior arma(cov = TRUE)",
     default_prior(bf(y ~ x + arma(time = t, gr = g, cov = TRUE)), data = dd))
show("default_prior ar()", default_prior(bf(y ~ x + ar(time = t, gr = g)),
                                          data = dd))
show("default_prior cosy()",
     default_prior(bf(y ~ x + cosy(time = t, gr = g)), data = dd))
cat("\n-- stancode ar(cov = TRUE) with ar \"\" (flat):\n")
lines_with(stancode(bf(y ~ x + ar(time = t, gr = g, cov = TRUE)), data = dd,
                    prior = set_prior("", class = "ar")), "lprior|ar;|ar\\]")

cat("\n== other classes under \"\"\n")
lines_with(stancode(bf(y ~ x + gp(x)), data = dd,
                    prior = c(set_prior("", class = "sdgp"),
                              set_prior("", class = "lscale"))),
           "lprior|sdgp|lscale")
dd$m <- rep(1:4, length.out = nrow(dd))
show("mo simo flat", validate_prior(set_prior("", class = "simo",
                                              coef = "mom1"),
                                    bf(y ~ mo(m)), data = dd))
lines_with(stancode(bf(y ~ mo(m)), data = dd,
                    prior = set_prior("", class = "simo", coef = "mom1")),
           "lprior|simo")
dd$yp <- abs(dd$y) + 0.1
lines_with(stancode(bf(yp ~ x), data = dd, family = Gamma("log"),
                    prior = set_prior("", class = "shape")),
           "lprior|shape")
lines_with(stancode(bf(y ~ x), data = dd, family = student(),
                    prior = set_prior("", class = "nu")),
           "lprior|nu")
dd$ord <- cut(dd$y, 4, labels = FALSE)
lines_with(stancode(bf(ord ~ x), data = dd, family = cumulative(),
                    prior = set_prior("", class = "Intercept")),
           "lprior|Intercept;|ordered")
lines_with(stancode(bf(ord ~ x), data = dd,
                    family = cumulative(threshold = "equidistant"),
                    prior = c(set_prior("", class = "Intercept"),
                              set_prior("", class = "delta"))),
           "lprior|delta")
lines_with(stancode(bf(y ~ x), data = dd, family = mixture(gaussian,
                                                            gaussian),
                    prior = c(set_prior("", class = "theta"),
                              set_prior("normal(0, 5)", class = "Intercept",
                                        dpar = "mu1"),
                              set_prior("normal(1, 5)", class = "Intercept",
                                        dpar = "mu2"))),
           "lprior|theta")

cat("\n== item 4: other spellings of one parameter\n")
show("b coef Intercept on a centered model", validate_prior(
  set_prior("normal(0, 1)", class = "b", coef = "Intercept"), bf(y ~ x),
  data = dd))
show("Intercept with nlpar", validate_prior(
  set_prior("normal(0, 1)", class = "Intercept", nlpar = "a"),
  bf(y ~ a * x, a ~ 1, nl = TRUE), data = dd))
show("sigma class when sigma has a formula", validate_prior(
  set_prior("normal(0, 1)", class = "sigma"), bf(y ~ x, sigma ~ x),
  data = dd))
show("Intercept dpar sigma without a formula", validate_prior(
  set_prior("normal(0, 1)", class = "Intercept", dpar = "sigma"), bf(y ~ x),
  data = dd))
show("default_prior gp", default_prior(bf(y ~ gp(x)), data = dd))
show("default_prior s(x)", default_prior(bf(y ~ s(x)), data = dd))

cat("\n== ported rows priors:5 and priors:123\n")
show("priors:5", sort(default_prior(
  count ~ zBase * Trt + (1|patient) + (1+Trt|visit), data = epilepsy,
  family = "poisson")$class))
p1 <- validate_prior(prior(normal(0,10), class = b) +
                       prior(cauchy(0,2), class = sd),
                     count ~ zAge + zBase * Trt + (1|patient),
                     data = epilepsy, family = poisson())
show("priors:123 nrow", nrow(p1))
print(p1)
