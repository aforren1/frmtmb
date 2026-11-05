# Area 5: default starts. Other links and bodies, global maximum against a
# reference found from many starts. Seed 7501. FRM_ARM=before for the base.
source("dev/review-nlfix-common.R")
ns <- asNamespace("frmtmb")
msh <- function(z) if (length(z)) substr(z[1], 1, 150) else ""
one <- function(lab, f, fam, d, ref = NULL, ...) {
  r <- tt(frm(f, family = fam, data = d, ...))
  if (iserr(r)) {
    cat(sprintf("%-30s ERROR %s\n", lab, substr(r$msg, 1, 150)))
    return(invisible(NULL))
  }
  ll <- as.numeric(logLik(r))
  cat(sprintf("%-30s ll=%.5f conv=%s ref=%s gap=%s beta=%s\n    MSG[%s]\n    WARN[%s]\n", lab, ll,
              r$opt$convergence, if (is.null(ref)) "NA" else format(ref, digits = 9),
              if (is.null(ref)) "NA" else format(ref - ll, digits = 3),
              paste(format(r$estimates$beta, digits = 4), collapse = ","), msh(attr(r, "m")), msh(attr(r, "w"))))
  invisible(r)
}
# many random starts: the best maximum frm() can find for this body
best_ms <- function(f, fam, d, nms, k = 12, scale = 1, seed = 1) {
  set.seed(seed)
  b <- -Inf
  for (i in seq_len(k)) {
    st <- setNames(rnorm(length(nms), 0, scale), nms)
    r <- tryCatch(suppressWarnings(suppressMessages(frm(f, family = fam, data = d, start = list(beta = st)))), error = function(e) NULL)
    if (!is.null(r) && is.finite(as.numeric(logLik(r))) && as.numeric(logLik(r)) > b) b <- as.numeric(logLik(r))
  }
  b
}
set.seed(7501)
n <- 300
d <- data.frame(x = runif(n, 0.1, 3), z = rnorm(n), g = factor(rep(1:15, each = 20)))
d$xz <- c(0, d$x[-1])
d$y <- rgamma(n, shape = 4, rate = 4 / (2 * exp(0.3 * d$x)))
d$yig <- d$y
cat("\n--- Gamma inverse (default link)\n")
fam <- Gamma()
lin <- as.numeric(logLik(frm(y ~ x, family = fam, data = d)))
one("b0 + b1 * x", bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), fam, d, ref = lin)
ref <- best_ms(bf(y ~ b0 * exp(b1 * x), b0 + b1 ~ 1, nl = TRUE), fam, d, c("b0_(Intercept)", "b1_(Intercept)"), scale = 0.5)
one("b0 * exp(b1 * x)", bf(y ~ b0 * exp(b1 * x), b0 + b1 ~ 1, nl = TRUE), fam, d, ref = ref)
ref <- best_ms(bf(y ~ a * x^b, a + b ~ 1, nl = TRUE), fam, d, c("a_(Intercept)", "b_(Intercept)"), scale = 0.5)
one("a * x^b", bf(y ~ a * x^b, a + b ~ 1, nl = TRUE), fam, d, ref = ref)
one("a * xz^b (x = 0 in data)", bf(y ~ a * xz^b, a + b ~ 1, nl = TRUE), fam, d)
ref <- best_ms(bf(y ~ exp(a + b * x), a + b ~ 1, nl = TRUE), fam, d, c("a_(Intercept)", "b_(Intercept)"), scale = 0.5)
one("exp(a + b * x)", bf(y ~ exp(a + b * x), a + b ~ 1, nl = TRUE), fam, d, ref = ref)
one("a + b (two intercepts)", bf(y ~ a + b, a + b ~ 1, nl = TRUE), fam, d, ref = as.numeric(logLik(frm(y ~ 1, family = fam, data = d))))
one("a + b1 * b2 * x (saddle)", bf(y ~ a + b1 * b2 * x, a + b1 + b2 ~ 1, nl = TRUE), fam, d, ref = lin)
one("a * b (no unit slope)", bf(y ~ a * b, a + b ~ 1, nl = TRUE), fam, d, ref = as.numeric(logLik(frm(y ~ 1, family = fam, data = d))))
one("a + b * x, a ~ 0 + z", bf(y ~ a + b * x, a ~ 0 + z, b ~ 1, nl = TRUE), fam, d)
one("1 / (a + b * x), a body undef", bf(y ~ 1 / (a + b * x), a + b ~ 1, nl = TRUE), fam, d)
one("a - exp(b * x)", bf(y ~ a - exp(b * x), a + b ~ 1, nl = TRUE), fam, d)
one("a + b * z, z both signs", bf(y ~ a + b * z, a + b ~ 1, nl = TRUE), fam, d, ref = as.numeric(logLik(frm(y ~ z, family = fam, data = d))))
one("(a + b * x) ^ 2", bf(y ~ (a + b * x) ^ 2, a + b ~ 1, nl = TRUE), fam, d)
one("RE: a + b * x, a ~ 1 + (1|g)", bf(y ~ a + b * x, a ~ 1 + (1 | g), b ~ 1, nl = TRUE), fam, d,
    ref = as.numeric(logLik(frm(y ~ x + (1 | g), family = fam, data = d))))
one("RE: b ~ 1 + (1|g)", bf(y ~ a + b * x, a ~ 1, b ~ 1 + (1 | g), nl = TRUE), fam, d,
    ref = as.numeric(logLik(frm(y ~ x + (0 + x | g), family = fam, data = d))))
one("nl mu, shape ~ z", bf(y ~ a + b * x, a + b ~ 1, shape ~ z, nl = TRUE), fam, d,
    ref = as.numeric(logLik(frm(bf(y ~ x, shape ~ z), family = fam, data = d))))
one("nl mu, nl shape", bf(y ~ a + b * x, nlf(shape ~ exp(s0 + s1 * z)), a + b + s0 + s1 ~ 1, nl = TRUE), fam, d,
    ref = as.numeric(logLik(frm(bf(y ~ x, shape ~ z), family = fam, data = d))))
one("mean placed by prior", bf(y ~ a + b * x, a + b ~ 1, nl = TRUE), fam, d, prior = set_prior("normal(0.2, 1)", nlpar = "a"))
one("start on b only", bf(y ~ a + b * x, a + b ~ 1, nl = TRUE), fam, d, start = list(beta = c("b_(Intercept)" = 0.01)))
one("start on a = 0 stays", bf(y ~ a + b * x, a + b ~ 1, nl = TRUE), fam, d, start = list(beta = c("a_(Intercept)" = 0)))
cat("\n--- inverse.gaussian 1/mu^2\n")
fam <- inverse.gaussian()
lin <- as.numeric(logLik(frm(y ~ x, family = fam, data = d)))
one("b0 + b1 * x", bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), fam, d, ref = lin)
one("b0 * exp(b1 * x)", bf(y ~ b0 * exp(b1 * x), b0 + b1 ~ 1, nl = TRUE), fam, d)
one("a * x^b", bf(y ~ a * x^b, a + b ~ 1, nl = TRUE), fam, d)
cat("\n--- other links\n")
d$yc <- rpois(n, 2 + d$x)
d$yb <- rbinom(n, 1, 0.3 + 0.1 * d$x)
one("poisson identity b0+b1x", bf(yc ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), poisson("identity"), d, ref = as.numeric(logLik(frm(yc ~ x, family = poisson("identity"), data = d))))
one("poisson identity a * x^b", bf(yc ~ a * x^b, a + b ~ 1, nl = TRUE), poisson("identity"), d)
one("bernoulli identity b0+b1x", bf(yb ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), bernoulli("identity"), d, ref = tryCatch(as.numeric(logLik(frm(yb ~ x, family = bernoulli("identity"), data = d))), error = function(e) NULL))
one("Gamma identity", bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), Gamma("identity"), d, ref = as.numeric(logLik(frm(y ~ x, family = Gamma("identity"), data = d))))
one("weibull identity", bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), weibull("identity"), d, ref = tryCatch(as.numeric(logLik(frm(y ~ x, family = weibull("identity"), data = d))), error = function(e) NULL))
one("exponential identity", bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), exponential("identity"), d, ref = tryCatch(as.numeric(logLik(frm(y ~ x, family = exponential("identity"), data = d))), error = function(e) NULL))
d$yneg <- d$y - 5
one("gaussian inverse neg y", bf(yneg ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), gaussian("inverse"), d, ref = tryCatch(as.numeric(logLik(frm(yneg ~ x, family = gaussian("inverse"), data = d))), error = function(e) NULL))
one("student inverse", bf(yneg ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), student("inverse"), d)
cat("\n--- multi-dpar and mixtures under the inverse link\n")
one("mixture(Gamma, Gamma), inverse", bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), mixture(Gamma(), Gamma()), d)
one("hurdle_gamma inverse", bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), tryCatch(hurdle_gamma("inverse"), error = function(e) hurdle_gamma()), d)
one("Gamma mv", mvbf(bf(y ~ a + b * x, a + b ~ 1, nl = TRUE), bf(yig ~ c + e * x, c + e ~ 1, nl = TRUE), rescor = FALSE), Gamma(), d)
cat("\n--- saddles: the same default start under gaussian and under an ordinal family\n")
set.seed(7502)
d2 <- data.frame(x = rnorm(n), z = rnorm(n))
d2$yg <- 1 + 2 * exp(0.5 * d2$x) + rnorm(n, 0, 0.5)
d2$yo <- cut(0.8 * d2$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = 1:4, ordered_result = TRUE)
d2$yc <- rpois(n, exp(0.2 + 0.5 * d2$x))
ref <- best_ms(bf(yg ~ b0 + b1 * exp(b2 * x), b0 + b1 + b2 ~ 1, nl = TRUE), gaussian(), d2, c("b0_(Intercept)", "b1_(Intercept)", "b2_(Intercept)"))
one("gauss b0 + b1 exp(b2 x)", bf(yg ~ b0 + b1 * exp(b2 * x), b0 + b1 + b2 ~ 1, nl = TRUE), gaussian(), d2, ref = ref)
ref <- best_ms(bf(yg ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), gaussian(), d2, c("b1_(Intercept)", "b2_(Intercept)"))
one("gauss b1 exp(b2 x)", bf(yg ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), gaussian(), d2, ref = ref)
for (fm in c("cumulative", "sratio", "cratio", "acat")) {
  F <- get(fm)
  ref <- best_ms(bf(yo ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), F(), d2, c("b1_(Intercept)", "b2_(Intercept)"))
  one(paste(fm, "b1 exp(b2 x)"), bf(yo ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), F(), d2, ref = ref)
  ref <- best_ms(bf(yo ~ a * x^3 + b * exp(c * x), a + b + c ~ 1, nl = TRUE), F(), d2, c("a_(Intercept)", "b_(Intercept)", "c_(Intercept)"))
  one(paste(fm, "a x^3 + b exp(c x)"), bf(yo ~ a * x^3 + b * exp(c * x), a + b + c ~ 1, nl = TRUE), F(), d2, ref = ref)
}
one("cumul a * (x - c)^2", bf(yo ~ a * (x - c)^2, a + c ~ 1, nl = TRUE), cumulative(), d2, ref = as.numeric(logLik(frm(yo ~ x + I(x^2), family = cumulative(), data = d2))))
one("cumul exp(a) * x", bf(yo ~ exp(a) * x, a ~ 1, nl = TRUE), cumulative(), d2, ref = as.numeric(logLik(frm(yo ~ x, family = cumulative(), data = d2))))
one("cumul exp(a + b * x)", bf(yo ~ exp(a + b * x), a + b ~ 1, nl = TRUE), cumulative(), d2)
one("poisson exp(a) * exp(b x)", bf(yc ~ exp(a) * exp(b * x), a + b ~ 1, nl = TRUE), poisson("identity"), d2,
    ref = as.numeric(logLik(frm(yc ~ x, family = poisson(), data = d2))))
