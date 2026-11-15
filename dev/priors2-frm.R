# Lane priors2: frmtmb's side of every case, before (rellib-r21) and after
# (the lane library).
# PRIORS2_LIB unset = lane library; PRIORS2_LIB=none = rellib-r21 only.
# Rscript dev/priors2-frm.R > dev/priors2-out/frm-<arm>.txt
lib <- Sys.getenv("PRIORS2_LIB", "C:/Users/adf44/source/r/wt-priors2-lib")
if (identical(lib, "none")) lib <- ""
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "\n")
options(width = 150, digits = 10)
cols <- c("prior", "class", "coef", "group", "resp", "dpar", "nlpar", "lb",
          "ub", "source")
show <- function(label, expr) {
  cat("\n=== ", label, "\n", sep = "")
  r <- tryCatch(withCallingHandlers(force(expr),
         warning = function(w) {
           cat("WARN: ", substr(gsub("\n", " ", conditionMessage(w)), 1, 200),
               "\n")
           invokeRestart("muffleWarning")
         },
         message = function(m) {
           cat("MSG: ", substr(gsub("\n", " ", conditionMessage(m)), 1, 200),
               "\n")
           invokeRestart("muffleMessage")
         }),
       error = function(e) {
         cat("ERROR: ", substr(gsub("\n", " ", conditionMessage(e)), 1, 500),
             "\n")
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
set.seed(1)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 12)))
d$y <- 1 + d$z * 0.3 + rnorm(10)[d$g] * 0.3 + rnorm(n, sd = 0.3)
d$y2 <- d$y + rnorm(n)
d$yneg <- d$y - 0.8 * d$z
d$o <- cut(d$y, 3, labels = FALSE)
# two well-separated components, 30 percent in the first
set.seed(3)
dm <- data.frame(y = c(rnorm(48, 0, 1), rnorm(112, 4, 1)), x = rnorm(160))
dm$y2 <- dm$y + rnorm(160)
m2 <- mixture(gaussian(), gaussian())
m3 <- mixture(gaussian(), gaussian(), gaussian())
# the weights, the last component the reference (no weight has a formula)
wts <- function(f) {
  b <- f$estimates[["betad"]]
  th <- b[grepl("^(y2?_)?theta[0-9]+_", names(b))]
  if (!length(th)) return(numeric(0))
  v <- c(th, 0)
  round(exp(v) / sum(exp(v)), 6)
}

# M: the mixture weights row
show("M1 default_prior, two gaussian, fit route",
     default_prior(y ~ 1, data = dm, family = m2))
show("M1s default_prior, two gaussian, sample route",
     default_prior(y ~ 1, data = dm, family = m2, route = "sample"))
show("M2 default_prior, three gaussian",
     default_prior(y ~ 1, data = dm, family = m3))
show("M3 default_prior, theta1 ~ x",
     default_prior(bf(y ~ 1, theta1 ~ x), data = dm, family = m2))
show("M5 default_prior, ordinal mixture",
     default_prior(o ~ x, data = d,
                   family = mixture(cumulative(), cumulative())))
show("M6 default_prior, mv two mixtures",
     default_prior(mvbf(bf(y ~ 1), bf(y2 ~ 1)), data = dm, family = m2))
f0 <- frm(y ~ 1, data = dm, family = m2)
show("M8 frm, ML weights", wts(f0))
show("M8b frm, dirichlet(2) on theta", {
  f <- frm(y ~ 1, data = dm, family = m2,
           prior = set_prior("dirichlet(2)", class = "theta"))
  list(w = wts(f), obj_ml = f0$opt$objective, obj = f$opt$objective)
})
show("M8c frm, dirichlet(c(20, 1)) on theta", {
  f <- frm(y ~ 1, data = dm, family = m2,
           prior = set_prior("dirichlet(c(20, 1))", class = "theta"))
  list(w = wts(f), obj = f$opt$objective)
})
show("M9 frm three, dirichlet(c(1, 2, 3))", {
  f <- frm(y ~ 1, data = dm, family = m3,
           prior = set_prior("dirichlet(c(1, 2, 3))", class = "theta"))
  wts(f)
})
show("M10 dirichlet(c(1, 2)) on three: length",
     frm(y ~ 1, data = dm, family = m3,
         prior = set_prior("dirichlet(c(1, 2))", class = "theta")))
show("M11 normal(0, 1) on class theta of a mixture",
     frm(y ~ 1, data = dm, family = m2,
         prior = set_prior("normal(0, 1)", class = "theta")))
show("M12 dirichlet(2) on class theta1",
     frm(y ~ 1, data = dm, family = m2,
         prior = set_prior("dirichlet(2)", class = "theta1")))
show("M14 dirichlet(2) class theta, dpar theta1",
     frm(y ~ 1, data = dm, family = m2,
         prior = set_prior("dirichlet(2)", class = "theta", dpar = "theta1")))
show("M16 dirichlet(2) class theta, coef a",
     frm(y ~ 1, data = dm, family = m2,
         prior = set_prior("dirichlet(2)", class = "theta", coef = "a")))
show("M17 dirichlet(2) class theta, group g",
     frm(y ~ 1, data = dm, family = m2,
         prior = set_prior("dirichlet(2)", class = "theta", group = "g")))
show("M18 dirichlet(2) class theta, lb = 0",
     frm(y ~ 1, data = dm, family = m2,
         prior = set_prior("dirichlet(2)", class = "theta", lb = 0)))
show("M19 dirichlet(2) class theta when theta1 ~ x",
     frm(bf(y ~ 1, theta1 ~ x), data = dm, family = m2,
         prior = set_prior("dirichlet(2)", class = "theta")))
show("M20 dirichlet(2) class theta on a non-mixture",
     frm(y ~ 1, data = dm,
         prior = set_prior("dirichlet(2)", class = "theta")))
show("M21 mv: dirichlet(5) class theta resp y2", {
  f <- frm(mvbf(bf(y ~ 1), bf(y2 ~ 1)), data = dm, family = m2,
           prior = set_prior("dirichlet(5)", class = "theta", resp = "y2"))
  wts(f)
})
show("M22 mv: dirichlet(5) class theta no resp",
     frm(mvbf(bf(y ~ 1), bf(y2 ~ 1)), data = dm, family = m2,
         prior = set_prior("dirichlet(5)", class = "theta")))
show("M28 dirichlet(0.5) class theta on frm()",
     frm(y ~ 1, data = dm, family = m2,
         prior = set_prior("dirichlet(0.5)", class = "theta")))
show("M30 round trip: frm(prior = default_prior(fit)) == frm()", {
  tab <- default_prior(f0)
  f <- frm(y ~ 1, data = dm, family = m2, prior = tab)
  list(identical_par = identical(f$opt$par, f0$opt$par),
       obj = c(f0$opt$objective, f$opt$objective))
})
show("M31 brms prior object: prior(dirichlet(2), class = theta)", {
  if (requireNamespace("brms", quietly = TRUE)) {
    f <- frm(y ~ 1, data = dm, family = m2,
             prior = brms::prior(dirichlet(2), class = theta))
    wts(f)
  } else "brms absent"
})
show("M32 prior_summary of a dirichlet(2) fit", {
  f <- frm(y ~ 1, data = dm, family = m2,
           prior = set_prior("dirichlet(2)", class = "theta"))
  prior_summary(f)
})

# R: the covariance class, renamed
show("R1 default_prior toep(): class rows",
     default_prior(y ~ toep(0 + factor(rep(1:4, 30)) | g), data = d))
show("R2 set_prior class theta coef theta_1 on (1 | g)",
     frm(y ~ z + (1 | g), data = d,
         prior = set_prior("normal(0, 1)", class = "theta", coef = "theta_1")))
show("R3 set_prior class thetacov coef theta_1 on (1 | g)",
     frm(y ~ z + (1 | g), data = d,
         prior = set_prior("normal(0, 1)", class = "thetacov",
                           coef = "theta_1")))

# C: cor with resp (4767)
show("C2 frm cor resp y2",
     frm(mvbf(bf(y ~ 1 + (1 + z | g)), bf(y2 ~ 1 + (1 + z | g))), data = d,
         prior = set_prior("lkj(5)", class = "cor", group = "g",
                           resp = "y2")))
show("C3 frm cor resp y2, no group",
     frm(mvbf(bf(y ~ 1 + (1 + z | g)), bf(y2 ~ 1 + (1 + z | g))), data = d,
         prior = set_prior("lkj(5)", class = "cor", resp = "y2")))

# I: a bound on class Intercept (5646)
d$xs <- d$z + 10
set.seed(11)
d$yi <- 1 + 0.41 * d$xs + rnorm(n, sd = 0.3)
show("I0 lm", coef(lm(yi ~ xs, data = d)))
show("I1 frm lb = 3 on Intercept", {
  f <- frm(yi ~ xs, data = d,
           prior = set_prior("", class = "Intercept", lb = 3))
  e <- fixef(f)[, "Estimate"]
  c(e, centered = unname(e[1] + e[2] * mean(d$xs)), conv = f$opt$convergence)
})
show("I2 frm lb = 9 on Intercept (binds the centered one)", {
  f <- frm(yi ~ xs, data = d,
           prior = set_prior("", class = "Intercept", lb = 9))
  e <- fixef(f)[, "Estimate"]
  c(e, centered = unname(e[1] + e[2] * mean(d$xs)), conv = f$opt$convergence)
})
cat("\ndone\n")
