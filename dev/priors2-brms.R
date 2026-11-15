# brms 2.23.0 on the cases of lane priors2 (round of 2026-11-13): the
# mixture-weights row, its accepted spellings and its refusals; cor with
# resp; a bound on class Intercept.
# Run: Rscript dev/priors2-brms.R > dev/priors2-out/brms.txt
# No Stan compile: default_prior(), validate_prior() and stancode() only.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
cat("brms", as.character(packageVersion("brms")), "\n")
options(width = 150)
cols <- c("prior", "class", "coef", "group", "resp", "dpar", "nlpar", "lb",
          "ub", "source")
show <- function(label, expr) {
  cat("\n=== ", label, "\n", sep = "")
  r <- tryCatch(withCallingHandlers(force(expr),
         warning = function(w) {
           cat("WARN: ", conditionMessage(w), "\n")
           invokeRestart("muffleWarning")
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
grep_code <- function(code, pat) {
  cat(grep(pat, strsplit(code, "\n")[[1]], value = TRUE), sep = "\n")
}
set.seed(1)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 12)))
d$y <- 1 + d$z * 0.3 + rnorm(10)[d$g] * 0.3 + rnorm(n, sd = 0.3)
d$y2 <- d$y + rnorm(n)
d$cnt <- rpois(n, 3)
d$o <- cut(d$y, 3, labels = FALSE)
m2 <- mixture(gaussian, gaussian)
m3 <- mixture(gaussian, gaussian, gaussian)

# M: the mixture weights row
show("M1 default_prior, two gaussian", default_prior(y ~ 1, data = d,
                                                     family = m2))
show("M2 default_prior, three gaussian", default_prior(y ~ 1, data = d,
                                                       family = m3))
show("M3 default_prior, theta1 ~ x", default_prior(bf(y ~ 1, theta1 ~ x),
                                                   data = d, family = m2))
show("M4 default_prior, theta1 = 0.4 (fixed)",
     default_prior(bf(y ~ 1, theta1 = 0.4), data = d, family = m2))
show("M5 default_prior, ordinal mixture",
     default_prior(o ~ x, data = d, family = mixture(cumulative, cumulative)))
show("M6 default_prior, mv two mixtures",
     default_prior(mvbf(bf(y ~ 1), bf(y2 ~ 1)), data = d, family = m2))
show("M7 default_prior, mixture(poisson, poisson)",
     default_prior(cnt ~ 1, data = d, family = mixture(poisson, poisson)))
show("M8 stancode, dirichlet(2) on theta",
     grep_code(stancode(y ~ 1, data = d, family = m2,
                        prior = set_prior("dirichlet(2)", class = "theta")),
               "theta"))
show("M8b standata con_theta, dirichlet(2)",
     standata(y ~ 1, data = d, family = m2,
              prior = set_prior("dirichlet(2)", class = "theta"))$con_theta)
show("M9 stancode, dirichlet(c(1, 2, 3)) on theta, three",
     grep_code(stancode(y ~ 1, data = d, family = m3,
                        prior = set_prior("dirichlet(c(1, 2, 3))",
                                          class = "theta")), "theta"))
show("M9b standata con_theta, c(1, 2, 3)",
     standata(y ~ 1, data = d, family = m3,
              prior = set_prior("dirichlet(c(1, 2, 3))",
                                class = "theta"))$con_theta)
show("M10 dirichlet(c(1, 2)) on three: length mismatch",
     stancode(y ~ 1, data = d, family = m3,
              prior = set_prior("dirichlet(c(1, 2))", class = "theta")))
show("M11 normal(0, 1) on class theta",
     stancode(y ~ 1, data = d, family = m2,
              prior = set_prior("normal(0, 1)", class = "theta")))
show("M12 dirichlet(2) on class theta1",
     stancode(y ~ 1, data = d, family = m2,
              prior = set_prior("dirichlet(2)", class = "theta1")))
show("M13 normal(0, 1) on class theta1",
     stancode(y ~ 1, data = d, family = m2,
              prior = set_prior("normal(0, 1)", class = "theta1")))
show("M14 dirichlet(2) class theta, dpar theta1",
     stancode(y ~ 1, data = d, family = m2,
              prior = set_prior("dirichlet(2)", class = "theta",
                                dpar = "theta1")))
show("M15 dirichlet(2) class Intercept, dpar theta1",
     stancode(y ~ 1, data = d, family = m2,
              prior = set_prior("dirichlet(2)", class = "Intercept",
                                dpar = "theta1")))
show("M16 dirichlet(2) class theta, coef a",
     stancode(y ~ 1, data = d, family = m2,
              prior = set_prior("dirichlet(2)", class = "theta",
                                coef = "a")))
show("M17 dirichlet(2) class theta, group g",
     stancode(y ~ 1, data = d, family = m2,
              prior = set_prior("dirichlet(2)", class = "theta",
                                group = "g")))
show("M18 dirichlet(2) class theta, lb = 0",
     stancode(y ~ 1, data = d, family = m2,
              prior = set_prior("dirichlet(2)", class = "theta", lb = 0)))
show("M19 dirichlet(2) class theta when theta1 ~ x",
     stancode(bf(y ~ 1, theta1 ~ x), data = d, family = m2,
              prior = set_prior("dirichlet(2)", class = "theta")))
show("M20 dirichlet(2) class theta on a non-mixture",
     stancode(y ~ 1, data = d,
              prior = set_prior("dirichlet(2)", class = "theta")))
show("M21 mv: dirichlet(2) class theta resp y2",
     grep_code(stancode(mvbf(bf(y ~ 1), bf(y2 ~ 1)), data = d, family = m2,
                        prior = set_prior("dirichlet(2)", class = "theta",
                                          resp = "y2")), "theta"))
show("M22 mv: dirichlet(2) class theta no resp",
     grep_code(stancode(mvbf(bf(y ~ 1), bf(y2 ~ 1)), data = d, family = m2,
                        prior = set_prior("dirichlet(2)", class = "theta")),
               "theta"))
show("M23 validate_prior dirichlet(2) theta",
     validate_prior(set_prior("dirichlet(2)", class = "theta"), y ~ 1,
                    data = d, family = m2))
show("M24 prior(dirichlet(2), theta) via prior()",
     grep_code(stancode(y ~ 1, data = d, family = m2,
                        prior = prior(dirichlet(2), class = theta)), "theta"))
show("M25 dirichlet(2) class theta when theta1 = 0.4 fixed",
     stancode(bf(y ~ 1, theta1 = 0.4), data = d, family = m2,
              prior = set_prior("dirichlet(2)", class = "theta")))
show("M26 stancode, no prior: theta lines",
     grep_code(stancode(y ~ 1, data = d, family = m2), "theta"))
show("M27 ordinal mixture stancode dirichlet(3)",
     grep_code(stancode(o ~ x, data = d,
                        family = mixture(cumulative, cumulative),
                        prior = set_prior("dirichlet(3)", class = "theta")),
               "theta"))
show("M28 dirichlet(0.5) class theta",
     grep_code(stancode(y ~ 1, data = d, family = m2,
                        prior = set_prior("dirichlet(0.5)", class = "theta")),
               "theta"))
show("M29 empty prior on class theta (flat)",
     grep_code(stancode(y ~ 1, data = d, family = m2,
                        prior = set_prior("", class = "theta")), "theta"))
show("M30 class theta, nlpar a",
     stancode(y ~ 1, data = d, family = m2,
              prior = set_prior("dirichlet(2)", class = "theta",
                                nlpar = "a")))

# C: cor with resp (4767)
show("C1 default_prior mv cor rows",
     default_prior(mvbf(bf(y ~ 1 + (1 | p | g)), bf(y2 ~ 1 + (1 | p | g))),
                   data = d))
show("C2 cor resp y2",
     stancode(mvbf(bf(y ~ 1 + (1 + z | g)), bf(y2 ~ 1 + (1 + z | g))),
              data = d, prior = set_prior("lkj(5)", class = "cor",
                                          group = "g", resp = "y2")))
show("C3 cor resp y2 no group",
     stancode(mvbf(bf(y ~ 1 + (1 + z | g)), bf(y2 ~ 1 + (1 + z | g))),
              data = d, prior = set_prior("lkj(5)", class = "cor",
                                          resp = "y2")))
show("C4 cor resp y on univariate",
     stancode(y ~ 1 + (1 + z | g), data = d,
              prior = set_prior("lkj(5)", class = "cor", resp = "y")))
show("C5 cor resp y2, shared id",
     stancode(mvbf(bf(y ~ 1 + (1 | p | g)), bf(y2 ~ 1 + (1 | p | g))),
              data = d, prior = set_prior("lkj(5)", class = "cor",
                                          group = "g", resp = "y2")))
show("C6 cor no resp mv",
     grep_code(stancode(mvbf(bf(y ~ 1 + (1 + z | g)),
                             bf(y2 ~ 1 + (1 + z | g))),
                        data = d, prior = set_prior("lkj(5)", class = "cor",
                                                    group = "g")), "lkj"))

# I: a bound on class Intercept (5646)
d$xs <- d$z + 10
d$yi <- 1 + 0.41 * d$xs + rnorm(n, sd = 0.3)
show("I1 stancode lb = 3 on Intercept, y ~ xs",
     grep_code(stancode(yi ~ xs, data = d,
                        prior = set_prior("", class = "Intercept", lb = 3)),
               "Intercept"))
show("I2 stancode normal(0, 5) lb = 3 on Intercept",
     grep_code(stancode(yi ~ xs, data = d,
                        prior = set_prior("normal(0, 5)", class = "Intercept",
                                          lb = 3)), "Intercept"))
show("I3 stancode lb = 3 on Intercept, center = FALSE",
     grep_code(stancode(bf(yi ~ xs, center = FALSE), data = d,
                        prior = set_prior("", class = "b",
                                          coef = "Intercept", lb = 3)),
               "Intercept|b"))
show("I4 stancode lb on Intercept, sigma ~ xs",
     grep_code(stancode(bf(yi ~ 1, sigma ~ xs), data = d,
                        prior = set_prior("", class = "Intercept",
                                          dpar = "sigma", lb = -1)),
               "Intercept"))
show("I5 default_prior y ~ xs",
     default_prior(yi ~ xs, data = d))
show("I6 uniform(0, 0.5) on Intercept without lb",
     grep_code(stancode(yi ~ xs, data = d,
                        prior = set_prior("uniform(0, 0.5)",
                                          class = "Intercept")), "Intercept"))
show("I7 y ~ 0 + Intercept + xs, lb on b coef Intercept",
     grep_code(stancode(yi ~ 0 + Intercept + xs, data = d,
                        prior = set_prior("", class = "b",
                                          coef = "Intercept", lb = 3)),
               "Intercept|real.*b"))

# backlog 1132: a mixture whose sigma is modeled
show("M40 default_prior, mixture, sigma1 ~ x",
     default_prior(bf(y ~ 1, sigma1 ~ x), data = d, family = m2))
show("M41 standata con_theta, dirichlet(c(1, 2, 3)) and c(1, 2) on three",
     tryCatch(standata(y ~ 1, data = d, family = m3,
                       prior = set_prior("dirichlet(c(1, 2))",
                                         class = "theta"))$con_theta,
              error = function(e) conditionMessage(e)))
