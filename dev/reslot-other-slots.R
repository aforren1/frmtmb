# Lane reslot, item 4: other parameters reachable through two class
# spellings. Usage: Rscript dev/reslot-other-slots.R before|after
# Seed 1 (re_data()). brms's answer for each spelling is in
# dev/reslot-out/brms.txt (dev/reslot-brms.R).
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
dd <- re_data(1L)
ri <- function(fit, pl) frmtmb:::resolve_prior_input(fit, pl)

cat("== 1. centered intercept: class Intercept and class b coef Intercept\n")
f1 <- frm(bf(y ~ x), family = gaussian(), data = dd)
a <- set_prior("normal(0, 1)", class = "Intercept")
b <- set_prior("normal(5, 1)", class = "b", coef = "Intercept")
try_show("Intercept alone", ri(f1, a))
try_show("b coef Intercept alone", ri(f1, b))
try_show("Intercept then b coef Intercept", ri(f1, a + b))
try_show("b coef Intercept then Intercept", ri(f1, b + a))

cat("\n== 2. nonlinear intercept: class b nlpar and class Intercept nlpar\n")
dd$yn <- 2 * exp(0.3 * dd$x) + rnorm(nrow(dd), 0, 0.3)
f2 <- frm(bf(yn ~ a * exp(k * x), a ~ 1, k ~ 1, nl = TRUE),
          family = gaussian(), data = dd,
          prior = set_prior("normal(2, 1)", nlpar = "a") +
            set_prior("normal(0, 1)", nlpar = "k"))
bn <- set_prior("normal(2, 1)", class = "b", nlpar = "a")
bi <- set_prior("normal(9, 1)", class = "Intercept", nlpar = "a")
bc <- set_prior("normal(7, 1)", class = "b", nlpar = "a", coef = "Intercept")
try_show("b nlpar a, then Intercept nlpar a", ri(f2, bn + bi))
try_show("Intercept nlpar a, then b nlpar a", ri(f2, bi + bn))
try_show("b nlpar a coef Intercept, then Intercept nlpar a",
         ri(f2, bc + bi))
try_show("Intercept nlpar a, then b nlpar a coef Intercept",
         ri(f2, bi + bc))
for (e in ri(f2, bc + bi)$entries) cat("  location", e$dist$location, "\n")
for (e in ri(f2, bi + bc)$entries) cat("  location", e$dist$location, "\n")

cat("\n== 3. sigma: own class and Intercept dpar sigma\n")
try_show("class sigma, no sigma formula",
         ri(f1, set_prior("normal(0, 1)", class = "sigma")))
try_show("Intercept dpar sigma, no sigma formula",
         ri(f1, set_prior("normal(0, 1)", class = "Intercept",
                          dpar = "sigma")))
f3 <- frm(bf(y ~ x, sigma ~ x), family = gaussian(), data = dd)
try_show("class sigma, sigma ~ x",
         ri(f3, set_prior("normal(0, 1)", class = "sigma")))

cat("\n== 4. gp(): sdgp / lscale and class theta\n")
f4 <- frm(bf(y ~ gp(x)), family = gaussian(), data = dd)
cat("theta:", length(f4$estimates$theta), "\n")
print(default_prior(bf(y ~ gp(x)), data = dd, family = gaussian()))
try_show("sdgp then theta_1",
         ri(f4, set_prior("student_t(3, 0, 1)", class = "sdgp") +
              set_prior("normal(0, 1)", class = "theta", coef = "theta_1")))

cat("\n== 5. smooth: class sd group s(x) and class theta\n")
f5 <- frm(bf(y ~ s(x)), family = gaussian(), data = dd)
print(default_prior(bf(y ~ s(x)), data = dd, family = gaussian()))
try_show("sd group s(x) then theta_1",
         ri(f5, set_prior("student_t(3, 0, 1)", class = "sd",
                          group = "s(x)") +
              set_prior("normal(0, 1)", class = "theta", coef = "theta_1")))

cat("\n== 6. a structure brms does not have: ar1(0 + factor(t) | g)\n")
dd$tf <- factor(dd$t)
f6 <- frm(bf(y ~ x + ar1(0 + tf | g)), family = gaussian(), data = dd)
cat("theta:", length(f6$estimates$theta), "\n")
print(default_prior(bf(y ~ x + ar1(0 + tf | g)), data = dd,
                    family = gaussian()))
try_show("cor then theta_2",
         ri(f6, set_prior("lkj(2)", class = "cor") +
              set_prior("normal(0, 1)", class = "theta", coef = "theta_2")))

cat("\n== 7. toep(): no LKJ, theta only\n")
f7 <- suppressWarnings(frm(bf(y ~ x + toep(0 + tf | g)), family = gaussian(),
                           data = dd))
print(default_prior(bf(y ~ x + toep(0 + tf | g)), data = dd,
                    family = gaussian()))
try_show("theta_3 on toep (a standard deviation)",
         ri(f7, set_prior("normal(0, 1)", class = "theta",
                          coef = "theta_3")))
try_show("class theta, no coef, on toep",
         ri(f7, set_prior("normal(0, 2)", class = "theta")))
try_show("theta_11, a toep correlation",
         ri(f7, set_prior("normal(0, 1)", class = "theta",
                          coef = "theta_11")))

cat("\n== 2b. the nonlinear intercept's two spellings as frm() reads them\n")
print(as.data.frame(bi))
r2 <- tryCatch(frm(bf(yn ~ a * exp(k * x), a ~ 1, k ~ 1, nl = TRUE),
                   family = gaussian(), data = dd,
                   prior = bc + bi + set_prior("normal(0, 1)", nlpar = "k")),
               error = function(e) conditionMessage(e))
cat("  frm(b nlpar a coef Intercept + Intercept nlpar a):",
    if (is.character(r2)) r2 else "fitted", "\n")
