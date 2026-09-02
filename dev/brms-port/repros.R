# Minimal repro for every FAILS-NEW item in dev/brms-vignette-port.md.
# Each block prints PASS (works today) or the exact error. Rerun after a
# fix to see an item flip.
suppressMessages(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-audit", quiet = TRUE))
set.seed(1)

show <- function(tag, expr) {
  r <- tryCatch({ force(expr); "PASS" },
                error = function(e) paste0("ERROR: ", conditionMessage(e)))
  cat(sprintf("%-6s %s\n", tag, substr(gsub("\n", " ", r), 1, 200)))
}

d <- data.frame(x = rnorm(60), g = gl(6, 10))
d$y <- 1 + 2 * d$x + rnorm(60)
d$o <- factor(sample(1:4, 60, TRUE), ordered = TRUE)
d$k <- as.integer(cut(d$y, 3))

show("F1", frm(y ~ x, data = d))
show("F2", frm(y ~ 1 + (1 | mm(g, g)), data = d, family = gaussian()))

fnl <- frm(bf(y ~ a * x + b, a ~ 1, b ~ 1, nl = TRUE), data = d,
           family = gaussian(), start = list(beta = c(1, 1)))
show("F3a", conditional_effects(fnl))
fmo <- frm(y ~ mo(o), data = d, family = gaussian())
show("F3b", conditional_effects(fmo))

d2 <- d; d2$z <- rnorm(60)
fsurf <- frm(y ~ t2(x, z), data = d2, family = gaussian())
show("F4", conditional_effects(fsurf, surface = TRUE))

flin <- frm(y ~ x, data = d, family = gaussian())
show("F5", hypothesis(flin, "x > 0"))
show("F6", frm(k ~ x, data = d, family = cumulative))
show("F7", frm(k ~ x, data = d, family = cumulative(), threshold = "equidistant"))
show("F8", frm(bf(y ~ a * x + b, a + b ~ 1, nl = TRUE), data = d,
               family = gaussian(), start = list(beta = c(1, 1))))
show("F9a", update(flin, formula. = ~ . + I(x^2)))
show("F9b", update(flin, formula = y ~ x, newdata = d))
show("F10", lf(sigma ~ x))

dmi <- d; dmi$x[1:5] <- NA
fmi <- frm(bf(y | mi() ~ mi(x)) + bf(x | mi() ~ 1) + set_rescor(FALSE),
           data = dmi, family = gaussian())
show("F13", conditional_effects(fmi, "x", resp = "y"))

show("F14", conditional_smooths(fsurf))
show("F15a", add_criterion(flin, "loo"))
show("F15b", bayes_R2(flin))

# frm_multiple post-processing (F11) and multivariate pp_check (F12)
imps <- lapply(1:3, function(i) { z <- d; z$y <- z$y + rnorm(60, 0, .1); z })
fm <- frm_multiple(y ~ x, data = imps, family = gaussian())
show("F11a", plot(fm, variable = "^b", regex = TRUE))
show("F11b", conditional_effects(fm, "x"))
d$y2 <- d$y + rnorm(60)
fmv <- frm(bf(mvbind(y, y2) ~ x) + set_rescor(TRUE), data = d, family = gaussian())
show("F12", pp_check(fmv, resp = "y"))
