# frmtmb's side of dev/priortab-brms.R, case for case.
# Before: PRIORTAB_LIB unset (rellib-r20). After: PRIORTAB_LIB = lane lib.
# Rscript dev/priortab-frm.R > dev/priortab-out/frm-<arm>.txt
lib <- Sys.getenv("PRIORTAB_LIB")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
cat("frmtmb.sample", dirname(find.package("frmtmb.sample")), "\n")
options(width = 150)
cols <- c("prior", "class", "coef", "group", "resp", "dpar", "nlpar", "lb",
          "ub", "source")
show <- function(label, expr) {
  cat("\n=== ", label, "\n", sep = "")
  r <- tryCatch(withCallingHandlers(force(expr),
         warning = function(w) {
           cat("WARN: ", conditionMessage(w), "\n")
           invokeRestart("muffleWarning")
         },
         message = function(m) {
           cat("MSG: ", conditionMessage(m))
           invokeRestart("muffleMessage")
         }),
       error = function(e) {
         cat("ERROR: ", conditionMessage(e), "\n")
         NULL
       })
  if (is.data.frame(r)) {
    print(as.data.frame(r)[, intersect(cols, names(r)), drop = FALSE],
          row.names = TRUE)
  } else if (!is.null(r)) {
    print(r)
  }
  invisible(r)
}
both <- function(label, expr) {
  e <- substitute(expr)
  for (rt in c("fit", "sample")) {
    show(paste0(label, " [", rt, "]"),
         eval(e, list(RT = rt), parent.frame()))
  }
}
set.seed(1)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 12)))
d$y <- 1 + d$z * 0.3 + rnorm(10)[d$g] * 0.3 + rnorm(n, sd = 0.3)
d$y2 <- d$y + rnorm(n)
d$Y_a <- d$y
d$Y_b <- d$y2
d$yp <- abs(d$y) + 0.1
d$f <- factor(sample(c("a", "b", "c"), n, TRUE))

both("1a gaussian y ~ z + (1 | g)",
     default_prior(y ~ z + (1 | g), data = d, route = RT))
both("1b student y ~ z",
     default_prior(y ~ z, data = d, family = student(), route = RT))
both("1c Gamma y ~ z",
     default_prior(yp ~ z, data = d, family = Gamma("log"), route = RT))
both("1d gaussian sigma ~ z",
     default_prior(bf(y ~ z, sigma ~ z), data = d, route = RT))
both("1e lognormal + (1 + z | g)",
     default_prior(yp ~ z + (1 + z | g), data = d, family = lognormal(),
                   route = RT))
both("1f zero_inflated_poisson",
     default_prior(yp ~ z, data = transform(d, yp = round(yp)),
                   family = zero_inflated_poisson(), route = RT))
both("1g beta", default_prior(pb ~ z, data = transform(d, pb = plogis(y)),
                              family = Beta(), route = RT))
both("1h negbinomial",
     default_prior(cnt ~ z, data = transform(d, cnt = rpois(n, 3)),
                   family = negbinomial(), route = RT))
both("1i skew_normal", default_prior(y ~ z, data = d, family = skew_normal(),
                                     route = RT))
both("1j hurdle_lognormal",
     default_prior(yp ~ z, data = d, family = hurdle_lognormal(),
                   route = RT))
both("1l shifted_lognormal",
     default_prior(yp ~ z, data = d, family = shifted_lognormal(),
                   route = RT))
both("1m exgaussian", default_prior(y ~ z, data = d, family = exgaussian(),
                                    route = RT))
both("1n von_mises", default_prior(vm ~ z,
  data = transform(d, vm = atan(y)), family = von_mises(), route = RT))
both("1o weibull", default_prior(yp ~ z, data = d, family = weibull(),
                                 route = RT))
both("1p gen_extreme_value", default_prior(y ~ z, data = d,
  family = gen_extreme_value(), route = RT))
both("1q asym_laplace", default_prior(y ~ z, data = d,
                                      family = asym_laplace(), route = RT))
both("1r zero_one_inflated_beta", default_prior(pb ~ z,
  data = transform(d, pb = plogis(y)), family = zero_one_inflated_beta(),
  route = RT))
both("1s frechet", default_prior(yp ~ z, data = d, family = frechet(),
                                 route = RT))
both("1t mv sigma rows", default_prior(bf(mvbind(y, y2) ~ z + (1 | g)) +
                                         set_rescor(TRUE), data = d,
                                       route = RT))

both("2a mixture(gaussian, gaussian) y ~ 1",
     default_prior(bf(y ~ 1), data = d, family = mixture(gaussian, gaussian),
                   route = RT))
both("2b mixture y ~ 1 + (1 | g)",
     default_prior(y ~ 1 + (1 | g), data = d,
                   family = mixture(gaussian, gaussian), route = RT))
both("2c mixture sigma1 ~ x",
     default_prior(bf(y ~ 1, sigma1 ~ x), data = d,
                   family = mixture(gaussian, gaussian), route = RT))
both("2d mixture of three",
     default_prior(bf(y ~ 1), data = d,
                   family = mixture(gaussian, gaussian, gaussian),
                   route = RT))
both("2e mixture theta1 ~ x",
     default_prior(bf(y ~ 1, theta1 ~ x), data = d,
                   family = mixture(gaussian, gaussian), route = RT))
both("2g validate dirichlet(2) on theta",
     validate_prior(set_prior("dirichlet(2)", class = "theta"),
                    bf(y ~ 1), data = d,
                    family = mixture(gaussian, gaussian), route = RT))
both("2i mixture(poisson, poisson)",
     default_prior(cnt ~ 1, data = transform(d, cnt = rpois(n, 3)),
                   family = mixture(poisson, poisson), route = RT))

p <- default_prior(y ~ z, data = d)
show("3a default_prior + set_prior", {
  q <- p + set_prior("normal(0,1)", class = "b"); print(class(q)); q })
show("3b set_prior + default_prior", {
  q <- set_prior("normal(0,1)", class = "b") + p; print(class(q)); q })
show("3c table + table", { q <- p + p; print(class(q)); length(q) })
show("3d validate_prior table + set_prior",
     validate_prior(p + set_prior("normal(0,1)", class = "b"), y ~ z,
                    data = d))
show("3e c(table, set_prior)", {
  q <- c(p, set_prior("normal(0,1)", class = "b")); print(class(q)); q })

fm <- bf(mvbind(y, y2) ~ z) + set_rescor(TRUE)
show("4a rescor coef", validate_prior(set_prior("lkj(2)", class = "rescor",
                                                coef = "y__y2"), fm, data = d))
show("4b rescor group", validate_prior(set_prior("lkj(2)", class = "rescor",
                                                 group = "g"), fm, data = d))
show("4c rescor resp", validate_prior(set_prior("lkj(2)", class = "rescor",
                                                resp = "y"), fm, data = d))

show("5a categorical cor dpar = mub",
     validate_prior(set_prior("lkj(2)", class = "cor", dpar = "mub"),
                    bf(f ~ x + (1 | p | g)), data = d,
                    family = categorical()))
show("5b categorical cor rows",
     default_prior(bf(f ~ x + (1 | p | g)), data = d, family = categorical()))
show("5c gaussian cor dpar = sigma",
     validate_prior(set_prior("lkj(2)", class = "cor", dpar = "sigma"),
                    bf(y ~ x + (1 + x | g), sigma ~ (1 + x | g)), data = d))
show("5d dpar block cor rows, no id",
     default_prior(bf(y ~ x + (1 + x | g), sigma ~ (1 + x | g)), data = d))
show("5e dpar block cor rows, id p",
     default_prior(bf(y ~ x + (1 | p | g), sigma ~ (1 | p | g)), data = d))
show("5f cor nlpar",
     validate_prior(set_prior("lkj(2)", class = "cor", nlpar = "a"),
                    bf(y ~ a + b * x, a ~ 1 + (1 + x | g), b ~ 1, nl = TRUE),
                    data = d))
show("5g nl cor rows",
     default_prior(bf(y ~ a + b * x, a ~ 1 + (1 + x | g), b ~ 1, nl = TRUE),
                   data = d))
show("5h mv cor rows", default_prior(bf(mvbind(y, y2) ~ x + (1 + x | g)),
                                     data = d))
show("5i mv cor resp = y2",
     validate_prior(set_prior("lkj(5)", class = "cor", group = "g",
                              resp = "y2"),
                    bf(mvbind(y, y2) ~ x + (1 + x | g)), data = d))
show("5j mv cor resp = y2, no group",
     validate_prior(set_prior("lkj(5)", class = "cor", resp = "y2"),
                    bf(mvbind(y, y2) ~ x + (1 + x | g)), data = d))
show("5k mv cor rows, shared id",
     default_prior(bf(mvbind(y, y2) ~ x + (1 + x | p | g)), data = d))
show("5l cor group only",
     validate_prior(set_prior("lkj(5)", class = "cor", group = "g"),
                    bf(y ~ x + (1 + x | g), sigma ~ (1 + x | g)), data = d))

show("6a mvbind(Y_a, Y_b)", default_prior(bf(mvbind(Y_a, Y_b) ~ z), data = d))
show("6b resp = Y_a", validate_prior(set_prior("normal(0, 1)", class = "b",
                                               resp = "Y_a"),
                                     bf(mvbind(Y_a, Y_b) ~ z), data = d))
show("6c resp = Ya", validate_prior(set_prior("normal(0, 1)", class = "b",
                                              resp = "Ya"),
                                    bf(mvbind(Y_a, Y_b) ~ z), data = d))

show("7a (1 | g) + (0 + x | g)",
     default_prior(y ~ (1 | g) + (0 + x | g), data = d))

show("9a ar(p = 2)", default_prior(y ~ z + ar(p = 2, gr = g), data = d))
show("9b ar(p = 1)", default_prior(y ~ z + ar(p = 1, gr = g), data = d))
show("9c uniform(0, 0.2) on b without bounds",
     validate_prior(set_prior("uniform(0, 0.2)", class = "b", coef = "z"),
                    y ~ z, data = d))
show("9d bound with coef", validate_prior(
  set_prior("", class = "b", coef = "z", lb = 0), y ~ z, data = d))
show("9e list(sp1) to validate_prior",
     validate_prior(list(set_prior("normal(0,1)", class = "b")), y ~ z,
                    data = d))
show("9f table column selection", {
  pe <- default_prior(y ~ z + x, data = d)
  validate_prior(pe[, "prior", drop = FALSE], y ~ z + x, data = d)
})
show("9g restated lb = 0 on sd with coef", validate_prior(
  set_prior("", class = "sd", group = "g", coef = "Intercept", lb = 0),
  y ~ z + (1 | g), data = d))
show("9h restated lb = 0 on sd class row", validate_prior(
  set_prior("", class = "sd", lb = 0), y ~ z + (1 | g), data = d))
show("9i restated lb = 0 on sigma", validate_prior(
  set_prior("", class = "sigma", lb = 0), y ~ z, data = d))
show("9j lb = 0 on sd group row", validate_prior(
  set_prior("", class = "sd", group = "g", lb = 0), y ~ z + (1 | g),
  data = d))
show("9l default table passed back", validate_prior(
  default_prior(y ~ z + (1 | g), data = d), y ~ z + (1 | g), data = d))
cat("\ndone\n")
