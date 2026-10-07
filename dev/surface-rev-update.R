# Reviewer of lane surface, claim 3: update() and the stored prior.
#   Rscript dev/surface-rev-update.R lane|base
source("dev/surface-rev-env.R")
arm <- rev_env(commandArgs(TRUE)[1])
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "\n")
set.seed(11)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n),
                g = factor(sample(letters[1:12], n, TRUE)),
                f = factor(sample(c("A", "B", "C"), n, TRUE)))
u <- rnorm(12, 0, 0.7)[d$g]
d$y <- 1 + 0.8 * d$x - 0.5 * d$z + 0.3 * d$w + u + rnorm(n)
ll <- function(f) as.numeric(logLik(f))

## a: no prior in the original fit
f0 <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
cat("a: class(f0$prior):", class(f0$prior), "| length",
    length(f0$prior), "\n")
ua <- rev_show("a update(no prior, + z)", update(f0, ~ . + z))
cat("   logLik update - direct:",
    ll(ua) - ll(frm(bf(y ~ x + z + (1 | g)), family = gaussian(),
                    data = d)), "\n")

## c: a typo in the ORIGINAL fit's prior
rev_show("c frm(prior coef typo 'xx')",
         frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d,
             prior = set_prior("normal(0, 1)", class = "b", coef = "xx")))
rev_show("c frm(prior on class cor, no correlation)",
         frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d,
             prior = set_prior("lkj(2)", class = "cor")))
rev_show("c frm(prior group typo 'gg')",
         frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d,
             prior = set_prior("normal(0, 1)", class = "sd", group = "gg")))

## d: coef priors, one matches after the update
pr <- c(set_prior("normal(0, 0.1)", class = "b", coef = "x"),
        set_prior("normal(0, 0.1)", class = "b", coef = "z"),
        set_prior("student_t(3, 0, 1)", class = "sd", group = "g"))
fd <- frm(bf(y ~ x + z + (1 | g)), family = gaussian(), data = d,
          prior = pr)
ud <- rev_show("d update(- z): coef z prior unmatched", update(fd, ~ . - z))
if (!inherits(ud, "rev_err")) {
  fd2 <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d,
             prior = c(set_prior("normal(0, 0.1)", class = "b", coef = "x"), set_prior("student_t(3, 0, 1)", class = "sd", group = "g")))
  cat("   logLik update - direct(kept priors):", ll(ud) - ll(fd2),
      "| b_x", fixef(ud)["x", 1], fixef(fd2)["x", 1], "\n")
}

## e: group-level priors, coef in a group
pe <- c(set_prior("lkj(2)", class = "cor", group = "g"),
        set_prior("normal(0, 1)", class = "sd", group = "g"))
fe <- rev_show("e fit (1 + x | g) with cor g, sd g",
               frm(bf(y ~ x + (1 + x | g)), family = gaussian(), data = d,
                   prior = pe))
if (!inherits(fe, "rev_err")) {
  rev_show("e update to (1 | g)",
           update(fe, ~ . - (1 + x | g) + (1 | g)))
  rev_show("e update to (1 | f) (group gone)",
           update(fe, ~ . - (1 + x | g) + (1 | f)))
}

## f: nlpar priors after the nl formula changes
d$yn <- 2 * exp(0.3 * d$x) + rnorm(n, 0, 0.3)
pn <- c(set_prior("normal(2, 1)", nlpar = "a"),
        set_prior("normal(0, 1)", nlpar = "b"))
fn <- rev_show("f nl fit a * exp(b x)",
               frm(bf(yn ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE),
                   family = gaussian(), data = d, prior = pn))
if (!inherits(fn, "rev_err")) {
  rev_show("f update(new nl body a + c x, nlpar b gone)",
           update(fn, bf(yn ~ a + c * x, a ~ 1, c ~ 1, nl = TRUE)))
  rev_show("f update(newdata)", update(fn, newdata = d[1:200, ]))
}

## h: a factor level vanishes from newdata
ph <- set_prior("normal(0, 0.1)", class = "b", coef = "fC")
fh <- frm(bf(y ~ x + f), family = gaussian(), data = d, prior = ph)
rev_show("h update(newdata without level C)",
         update(fh, newdata = droplevels(d[d$f != "C", ])))

## i: family change, sigma prior unmatched
fi <- frm(bf(y ~ x), family = gaussian(), data = d,
          prior = set_prior("exponential(1)", class = "sigma"))
d$cnt <- rpois(n, exp(0.2 + 0.3 * d$x))
rev_show("i update(family = poisson, y = cnt)",
         update(fi, cnt ~ x, family = poisson()))

## g: update(prior =) replaces the stored prior; brms merges
pg <- c(set_prior("normal(0, 0.05)", class = "b", coef = "x"),
        set_prior("normal(0, 0.5)", class = "sd", group = "g"))
fg <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d, prior = pg)
newp <- set_prior("exponential(1)", class = "sigma")
ug <- rev_show("g update(prior = sigma prior only)",
               update(fg, prior = newp))
if (!inherits(ug, "rev_err")) {
  merged <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d,
                prior = c(newp, pg))
  cat("   prior of update():", nrow(as.data.frame(prior_summary(ug))),
      "rows\n")
  print(prior_summary(ug))
  cat("   b_x: original", signif(fixef(fg)["x", 1], 5),
      "| update(prior = newp)", signif(fixef(ug)["x", 1], 5),
      "| brms-style merge", signif(fixef(merged)["x", 1], 5), "\n")
  cat("   sd_g: update", signif(VarCorr(ug)$g$sd[1], 5),
      "| merge", signif(VarCorr(merged)$g$sd[1], 5), "\n")
}
