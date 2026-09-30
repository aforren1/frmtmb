# Reviewer: each new refusal with its condition ABSENT, and a few
# option combinations. Lane library only.
#   Rscript dev/postfit2-rev-guards.R
.libPaths(c("C:/Users/adf44/source/r/wt-postfit2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
say <- function(...) cat(sprintf(...), "\n", sep = "")
f3 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
try_it <- function(label, expr) {
  r <- tryCatch(suppressMessages(expr), error = function(e) e)
  if (inherits(r, "error")) {
    say("%s: REFUSED: %s", label, substr(gsub("\n", " ",
                                              conditionMessage(r)), 1, 200))
  } else {
    say("%s: ok", label)
  }
  invisible(r)
}
set.seed(4)
d <- data.frame(z = runif(200), x = rnorm(200), g = factor(rep(1:10, 20)))
d$y <- 2 * exp(0.5 * d$z) + rnorm(10, 0, 0.3)[d$g] + rnorm(200, 0, 0.3)

## nonlinear Wald band, group terms
fa <- frm(bf(y ~ a * exp(b * z), a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
          family = gaussian(), data = d)
r <- try_it("nl, group in nlpar a, observed g=3, re_formula NULL, wald",
            conditional_effects(fa, "z", resolution = 5, re_formula = NULL,
                                conditions = data.frame(
                                  g = factor(3, levels = 1:10))))
if (!inherits(r, "error")) {
  nd <- r$z
  ref <- as.vector(frm_linpred(fa, newdata = nd, re_formula = NULL,
                               type = "response"))
  say("  estimate vs frm_linpred(re_formula = NULL): max rel gap %.3g",
      max(abs(r$z$estimate__ - ref)) / max(abs(ref)))
  pop <- conditional_effects(fa, "z", resolution = 5)$z
  say("  se observed g %s | se population %s", f3(r$z$se__), f3(pop$se__))
}
try_it("nl, group in nlpar a, new group (no conditions), wald",
       conditional_effects(fa, "z", resolution = 5, re_formula = NULL))
fs <- frm(bf(y ~ a * exp(b * z), a ~ 1, b ~ 1, sigma ~ (1 | g), nl = TRUE),
          family = gaussian(), data = d)
try_it("nl mu with NO group term; group only in sigma; re_formula NULL",
       conditional_effects(fs, "z", resolution = 5, re_formula = NULL))
try_it("same, the population curve", conditional_effects(fs, "z",
                                                          resolution = 5))
fl <- frm(bf(y ~ z + (1 | g), sigma ~ (1 | g)), family = gaussian(),
          data = d)
r1 <- try_it("linear mu, same shape, re_formula NULL (control)",
             conditional_effects(fl, "z", resolution = 5, re_formula = NULL))

## surface with the profile band, too_far = 0 (guard absent)
fx <- frm(bf(y ~ x * z), family = gaussian(), data = d)
rp <- try_it("surface + profile, too_far 0",
             conditional_effects(fx, "x:z", surface = TRUE, band = "profile",
                                 resolution = 4))
rw <- conditional_effects(fx, "x:z", surface = TRUE, resolution = 4)
if (!inherits(rp, "error")) {
  say("  profile vs wald on the surface: max rel gap lower %.3g",
      max(abs(rp[[1]]$lower__ - rw[[1]]$lower__)) /
        max(abs(rw[[1]]$lower__)))
}
try_it("surface + profile, too_far 0.1 (guard present)",
       conditional_effects(fx, "x:z", surface = TRUE, band = "profile",
                           resolution = 4, too_far = 0.1))
try_it("x:g surface (second is a factor), spaghetti boot",
       conditional_effects(fx, "x", surface = TRUE, spaghetti = TRUE,
                           band = "boot", boot = 5, seed = 1,
                           resolution = 4))
rs <- try_it("surface on a factor second predictor",
             conditional_effects(fl, "z:g", surface = TRUE, resolution = 4))
if (!inherits(rs, "error")) say("  surface attr %s, rows %d",
                                attr(rs[[1]], "surface"), nrow(rs[[1]]))
## surface with several condition rows
rc <- try_it("surface with two condition rows",
             conditional_effects(fx, "x:z", surface = TRUE, resolution = 4,
                                 conditions = data.frame(g = c("1", "2"))))
if (!inherits(rc, "error")) say("  rows %d, cond__ %s", nrow(rc[[1]]),
                                paste(unique(rc[[1]]$cond__), collapse = ","))
grDevices::pdf(NULL)
try_it("plot surface with two condition rows", plot(rc, ask = FALSE))
grDevices::dev.off()

## spaghetti on a boot band with conditions of two rows and a factor pair
fb <- frm(bf(y ~ x + g), family = gaussian(), data = d)
sp <- try_it("spaghetti boot, x:g, two condition rows? (none)",
             conditional_effects(fb, "x:g", spaghetti = TRUE, band = "boot",
                                 boot = 4, seed = 1, resolution = 3))
if (!inherits(sp, "error")) {
  s <- attr(sp[[1]], "spaghetti")
  say("  spaghetti rows %d = 4 draws x %d grid rows: %s", nrow(s),
      nrow(sp[[1]]), nrow(s) == 4 * nrow(sp[[1]]))
}
fo <- frm(bf(y ~ x + z), family = gaussian(), data = d)
cd <- data.frame(z = c(0.2, 0.8))
sp2 <- try_it("spaghetti boot with two condition rows",
              conditional_effects(fo, "x", spaghetti = TRUE, band = "boot",
                                  boot = 4, seed = 1, resolution = 3,
                                  conditions = cd))
if (!inherits(sp2, "error")) {
  s <- attr(sp2[[1]], "spaghetti")
  say("  spaghetti %s rows, frame %d rows, cond__ in spaghetti: %s",
      if (is.null(s)) "NULL" else nrow(s), nrow(sp2[[1]]),
      !is.null(s) && !is.null(s$cond__))
  if (!is.null(s)) {
    ok <- all(vapply(split(s, s$sample__), function(p) {
      all(tapply(p$estimate__, p$cond__, length) == 3)
    }, NA))
    say("  every draw has 3 points per condition: %s", ok)
  }
}
## conditional_smooths edge cases
d$f <- factor(rep(c("a", "b"), 100))
fsm <- frm(bf(y ~ f + s(z, by = f) + s(x, g, bs = "fs", k = 4),
              sigma ~ s(x)), family = gaussian(), data = d)
cs <- try_it("conditional_smooths: by, fs, sigma smooth",
             conditional_smooths(fsm, resolution = 5))
if (!inherits(cs, "error")) say("  keys: %s", paste(names(cs), collapse = " | "))
try_it("conditional_smooths smooths = 's(z,by=f)' (no spaces)",
       conditional_smooths(fsm, smooths = "s(z,by=f)", resolution = 5))
try_it("conditional_smooths smooths = 's(z, by=f)'",
       conditional_smooths(fsm, smooths = "s(z, by=f)", resolution = 5))
try_it("conditional_smooths ndraws on a fit (refusal by design)",
       conditional_smooths(fsm, ndraws = 5))
try_it("conditional_smooths on a model without smooths",
       conditional_smooths(fo))
try_it("make_conditions on a fit with a factor and numeric",
       make_conditions(fsm, c("x", "f")))
try_it("update_adterms on an mvbf (refusal by design)",
       update_adterms(mvbf(bf(y ~ x), bf(z ~ x)), ~ weights(w)))
try_it("update_adterms on a bf() with sigma", print(update_adterms(
  bf(y ~ x, sigma ~ x), ~ weights(w))))
say("done")
