# Reviewer, punch round 2: (a) the re_formula text in the reuse key,
# spellings of one formula and different formulas; (b) the stop on an
# unseen smooth-factor level: correct calls under re_formula = NA (and
# NULL) must not be refused.
#   Rscript dev/postfit2-rev-p2-keysmooth.R
.libPaths(c("C:/Users/adf44/source/r/wt-postfit2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
say <- function(...) cat(sprintf(...), "\n", sep = "")
try_it <- function(label, expr) {
  r <- tryCatch(suppressMessages(expr), error = function(e) e)
  say("%s: %s", label, if (inherits(r, "error"))
    paste("REFUSED:", substr(gsub("\n", " ", conditionMessage(r)), 1, 170)) else "ok")
  invisible(r)
}
## (a) keys
set.seed(62)
d <- data.frame(x = rnorm(240), g = factor(rep(1:12, 20)),
                h = factor(rep(1:8, 30)))
d$y <- rnorm(240, 1 + 0.5 * d$x + rnorm(12)[d$g] + rnorm(8, 0, 0.5)[d$h], 0.5)
f <- frm(bf(y ~ x + (1 | g) + (1 | h)), family = gaussian(), data = d)
cn <- list(g = "3", h = "2")
ceb <- function(re, boot) suppressMessages(conditional_effects(
  f, "x", resolution = 3, conditions = cn, re_formula = re, band = "boot",
  boot = boot, seed = 5))
a <- ceb(~ (1 | g), 10)
try_it("~(1|g) boot reused by ~ ( 1|g ) (spaces)", ceb(~ ( 1|g ), attr(a, "boot")))
ab <- ceb(~ (1 | g) + (1 | h), 10)
try_it("~(1|g)+(1|h) boot reused by ~(1|h)+(1|g) (term order)",
       ceb(~ (1 | h) + (1 | g), attr(ab, "boot")))
try_it("~(1|g)+(1|h) boot reused by NULL (the same terms)",
       ceb(NULL, attr(ab, "boot")))
try_it("~(1|g) boot reused by ~(1|h) (different formula)",
       ceb(~ (1 | h), attr(a, "boot")))
p <- suppressMessages(conditional_effects(f, "x", resolution = 3,
                                          conditions = cn, band = "boot",
                                          boot = 10, seed = 5))
try_it("NA boot reused by ~0 (population, other spelling)",
       ceb(~ 0, attr(p, "boot")))
say("key text of a long formula: %s", frmtmb:::ce_re_key(
  ~ (1 | g) + (1 | h) + (1 | g) + (1 | h) + (1 | g) + (1 | h) + (1 | g) +
    (1 | h) + (1 | g) + (1 | h) + (1 | g) + (1 | h) + (1 | g) + (1 | h)))

## (b) smooth-factor stop, correct calls
set.seed(63)
ds <- data.frame(x = runif(200), z = rnorm(200), g = factor(rep(1:5, 40)),
                 gc = rep(c("p", "q", "r", "s", "t"), 40))
ds$o <- factor(ds$g, ordered = TRUE)
ds$y <- sin(2 * pi * ds$x) + rnorm(5, 0, 0.5)[ds$g] + 0.3 * ds$z +
  rnorm(200, 0, 0.3)
fs <- frm(bf(y ~ z + s(x, g, bs = "fs", k = 5)), family = gaussian(), data = ds)
try_it("fs, default call (re_formula NA)", conditional_effects(fs))
try_it("fs, effects x:g (NA)", conditional_effects(fs, "x:g"))
try_it("fs, conditions g = 3 (NA)", conditional_effects(fs, "x", conditions = list(g = "3")))
try_it("fs, conditions g = factor 3 with other level order (NA)",
       conditional_effects(fs, "x", conditions = data.frame(g = factor("3", levels = c("5", "3", "1")))))
try_it("fs, conditions g = 3 numeric (NA)", conditional_effects(fs, "x", conditions = list(g = 3)))
try_it("fs, re_formula NULL, nothing set", conditional_effects(fs, "x", re_formula = NULL))
try_it("fs, make_conditions(g) (NA)", conditional_effects(fs, "x", conditions = make_conditions(fs, "g")))
fb <- frm(bf(y ~ g + s(x, by = g, k = 5)), family = gaussian(), data = ds)
try_it("by-factor smooth, default (NA)", conditional_effects(fb))
try_it("by-factor smooth, g = 4 (NA)", conditional_effects(fb, "x", conditions = list(g = "4")))
fo <- frm(bf(y ~ o + s(x, by = o, k = 5)), family = gaussian(), data = ds)
try_it("ordered by-factor smooth, default (NA)", conditional_effects(fo))
fc <- tryCatch(frm(bf(y ~ gc + s(x, by = gc, k = 5)), family = gaussian(), data = ds), error = function(e) { say("character by-variable fit: ERROR %s (mgcv; not this lane)", conditionMessage(e)); NULL })
if (!is.null(fc)) try_it("character by-variable smooth, default (NA)", conditional_effects(fc))
if (!is.null(fc)) try_it("character by-variable smooth, gc = r (NA)",
       conditional_effects(fc, "x", conditions = list(gc = "r")))
fz <- frm(bf(y ~ s(x, by = z, k = 5) + z), family = gaussian(), data = ds)
try_it("numeric by-variable smooth, default (NA)", conditional_effects(fz))
fsr <- frm(bf(y ~ z + s(x, g, bs = "fs", k = 5) + (1 | gc)), family = gaussian(), data = ds)
try_it("fs + (1 | gc), re_formula NULL, nothing set", conditional_effects(fsr, "x", re_formula = NULL))
try_it("fs + (1 | gc), re_formula NA, nothing set", conditional_effects(fsr, "x"))
fre <- frm(bf(y ~ z + s(g, bs = "re")), family = gaussian(), data = ds)
try_it("s(g, bs = 're'), default (NA)", conditional_effects(fre))
try_it("s(g, bs = 're'), re_formula NULL", conditional_effects(fre, "z", re_formula = NULL))
say("done")
