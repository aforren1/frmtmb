# Copy of dev/postfit2-rev-p1-attack.R for the ceplot review, with the library
# paths of this round (wt-ceplot-lib, rellib-r4); otherwise unchanged.
# Reviewer, punch round 1: the placeholder logic, the by-value reuse
# key, and allow_new_levels = TRUE on every non-population display.
# For each shape: the Wald band (which carries a new level's variance
# analytically), the boot band (fit) and the draws band (hand-built
# draws, the ML estimate plus N(0, 0.05) noise, 200 draws, seed 1).
#   Rscript dev/postfit2-rev-p1-attack.R <base|lane>
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "lane"
libs <- c("C:/Users/adf44/source/r/wt-ceplot-lib",
          "C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
say <- function(...) cat(sprintf(...), "\n", sep = "")
f3 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
say("ARM %s from %s", arm, find.package("frmtmb"))
hand <- function(fit, n = 200, seed = 1) {
  tpl <- fit$frame$par_template
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  set.seed(seed)
  M <- matrix(rep(est, each = n) + rnorm(n * length(est), 0, 0.05), n,
              dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
  structure(list(stanfit = NULL,
                 draws = cbind(frmtmb.sample:::draws_to_natural(M, fit),
                               lp__ = 0), fit = fit),
            class = "frmtmb_draws")
}
run <- function(label, expr) {
  r <- tryCatch(suppressWarnings(suppressMessages(expr)),
                error = function(e) e)
  if (inherits(r, "error")) {
    say("%s: ERROR %s", label, substr(gsub("\n", " ", conditionMessage(r)),
                                      1, 230))
    return(invisible(NULL))
  }
  d <- r[[1]]
  for (cv in unique(as.character(d$cond__))) {
    s <- d[as.character(d$cond__) == cv, ]
    say("%s [%s]: est %s | width %s", label, cv, f3(s$estimate__),
        f3(s$upper__ - s$lower__))
  }
  invisible(r)
}
three <- function(label, fit, ds, ...) {
  run(paste(label, "| wald"), conditional_effects(fit, "x", resolution = 3,
                                                  re_formula = NULL, ...))
  run(paste(label, "| boot40"), conditional_effects(
    fit, "x", resolution = 3, re_formula = NULL, band = "boot", boot = 40,
    seed = 5, ...))
  if (!is.null(ds)) {
    run(paste(label, "| draws"), conditional_effects(
      ds, "x", resolution = 3, re_formula = NULL, seed = 1, ...))
  }
}

## P1. a grouping variable that is also a fixed predictor:
## y ~ x + trt + (1 | trt:subj)
set.seed(41)
d1 <- expand.grid(rep = 1:4, trt = factor(c("a", "b")), subj = factor(1:12))
d1$x <- rnorm(nrow(d1))
d1$y <- rnorm(nrow(d1), 1 + 0.5 * d1$x + (d1$trt == "b") +
                rnorm(24, 0, 1)[interaction(d1$trt, d1$subj)], 0.5)
f1 <- frm(bf(y ~ x + trt + (1 | trt:subj)), family = gaussian(), data = d1)
ds1 <- hand(f1)
three("P1 trt:subj, nothing set", f1, ds1)
three("P1 trt:subj, trt=b", f1, ds1, conditions = list(trt = "b"))
three("P1 trt:subj, trt=b subj=3", f1, ds1,
      conditions = list(trt = "b", subj = "3"))

## P2. a grouping variable that is also a smooth's by variable
set.seed(42)
d2 <- data.frame(x = runif(240), g = factor(rep(1:6, 40)))
d2$y <- sin(2 * pi * d2$x) * (1 + 0.3 * as.integer(d2$g)) +
  rnorm(6, 0, 0.5)[d2$g] + rnorm(240, 0, 0.3)
f2 <- tryCatch(frm(bf(y ~ g + s(x, by = g, k = 5) + (1 | g)),
                   family = gaussian(), data = d2), error = function(e) e)
if (inherits(f2, "error")) say("P2 fit: ERROR %s", conditionMessage(f2)) else {
  three("P2 s(x, by = g) + (1 | g), nothing set", f2, hand(f2))
  three("P2 ..., g = 2", f2, NULL, conditions = list(g = "2"))
}

## P3. mm() membership term
set.seed(43)
d3 <- data.frame(x = rnorm(200), g1 = factor(sample(1:10, 200, TRUE)),
                 g2 = factor(sample(1:10, 200, TRUE)))
u <- rnorm(10, 0, 1)
d3$y <- rnorm(200, 1 + 0.5 * d3$x + 0.5 * (u[d3$g1] + u[d3$g2]), 0.5)
f3m <- tryCatch(frm(bf(y ~ x + (1 | mm(g1, g2))), family = gaussian(),
                    data = d3), error = function(e) e)
if (inherits(f3m, "error")) say("P3 fit: ERROR %s", conditionMessage(f3m)) else {
  say("P3 VarCorr sd %s", f3(unlist(lapply(VarCorr(f3m), function(v) v$sd[, 1]))))
  run("P3 mm, population | wald", conditional_effects(f3m, "x", resolution = 3))
  three("P3 mm(g1, g2), nothing set", f3m, hand(f3m))
  three("P3 mm, g1 = g2 = 2", f3m, NULL,
        conditions = list(g1 = "2", g2 = "2"))
}

## P4. gr(g, by = f)
set.seed(44)
d4 <- data.frame(x = rnorm(240), g = factor(rep(1:12, 20)),
                 f = factor(rep(c("a", "b"), each = 120)))
d4$y <- rnorm(240, 1 + 0.5 * d4$x + rnorm(24, 0, 1)[interaction(d4$g, d4$f)],
              0.5)
f4 <- tryCatch(frm(bf(y ~ x + f + (1 | gr(g, by = f))), family = gaussian(),
                   data = d4), error = function(e) e)
if (inherits(f4, "error")) say("P4 fit: ERROR %s", conditionMessage(f4)) else {
  say("P4 block group_name %s, levels %s ...",
      f4$frame$re_blocks[[1]]$group_name,
      paste(head(f4$frame$re_blocks[[1]]$levels, 3), collapse = ","))
  three("P4 gr(g, by = f), f = b", f4, hand(f4), conditions = list(f = "b"))
  three("P4 gr(g, by = f), f = b, g = 3", f4, NULL,
        conditions = list(f = "b", g = "3"))
}

## P5. crossed with an interaction term, both main levels set, their
## combination unseen: brms reads g and h and draws a new g:h
set.seed(45)
d5 <- expand.grid(rep = 1:3, g = factor(1:8), h = factor(1:8))
d5 <- d5[!(d5$g == "1" & d5$h == "2"), ]
d5$x <- rnorm(nrow(d5))
d5$y <- rnorm(nrow(d5), 1 + 0.5 * d5$x + rnorm(8, 0, 1)[d5$g] +
                rnorm(8, 0, 1)[d5$h] +
                rnorm(64, 0, 0.7)[interaction(d5$g, d5$h)], 0.5)
f5 <- frm(bf(y ~ x + (1 | g) + (1 | h) + (1 | g:h)), family = gaussian(),
          data = d5)
ds5 <- hand(f5)
three("P5 g=1, h=2 (g:h unseen)", f5, ds5,
      conditions = list(g = "1", h = "2"))
three("P5 g=1, h=3 (all observed)", f5, ds5,
      conditions = list(g = "1", h = "3"))

## P6. multivariate, one block shared by |p| across responses
set.seed(46)
d6 <- data.frame(x = rnorm(200), g = factor(rep(1:20, 10)))
u1 <- rnorm(20, 0, 1)
d6$y1 <- rnorm(200, 1 + 0.5 * d6$x + u1[d6$g], 0.5)
d6$y2 <- rnorm(200, -1 + 0.3 * d6$x + 0.8 * u1[d6$g], 0.5)
f6 <- frm(mvbf(bf(y1 ~ x + (1 | p | g)), bf(y2 ~ x + (1 | p | g))),
          family = gaussian(), data = d6)
ds6 <- hand(f6)
for (r in c("y1", "y2")) {
  three(paste("P6 mv |p|, nothing set, resp", r), f6, ds6, resp = r)
  three(paste("P6 mv |p|, g = 4, resp", r), f6, ds6, resp = r,
        conditions = list(g = "4"))
}

## P7. a factor-smooth's factor set to a level the fit never saw
set.seed(47)
d7 <- data.frame(x = runif(200), g = factor(rep(1:5, 40)))
d7$y <- sin(2 * pi * d7$x) + rnorm(5, 0, 0.5)[d7$g] + rnorm(200, 0, 0.3)
f7 <- frm(bf(y ~ s(x, g, bs = "fs", k = 5)), family = gaussian(), data = d7)
for (cnd in list(list(g = "2"), list(g = "99"))) {
  run(sprintf("P7 fs, g = %s, re_formula NULL | wald", cnd$g),
      conditional_effects(f7, "x", resolution = 3, re_formula = NULL,
                          conditions = cnd))
  run(sprintf("P7 fs, g = %s, re_formula NA | wald", cnd$g),
      conditional_effects(f7, "x", resolution = 3, conditions = cnd))
}

## P8. the by-value reuse key: (1 + x | g), g set, re_formula NULL and
## ~ (1 | g) read the same grid and hold the same block, and predict
## different curves
set.seed(48)
d8 <- data.frame(x = rnorm(240), g = factor(rep(1:12, 20)))
b0 <- rnorm(12, 0, 1)
b1 <- rnorm(12, 0, 1)
d8$y <- rnorm(240, 1 + b0[d8$g] + (0.5 + b1[d8$g]) * d8$x, 0.5)
f8 <- frm(bf(y ~ x + (1 + x | g)), family = gaussian(), data = d8)
c8 <- list(g = "3")
a <- suppressMessages(conditional_effects(f8, "x", resolution = 3,
                                          re_formula = NULL, conditions = c8,
                                          band = "boot", boot = 20, seed = 5))
b <- tryCatch(suppressMessages(conditional_effects(
  f8, "x", resolution = 3, re_formula = ~ (1 | g), conditions = c8,
  band = "boot", boot = attr(a, "boot"))), error = function(e) e)
bf <- suppressMessages(conditional_effects(
  f8, "x", resolution = 3, re_formula = ~ (1 | g), conditions = c8,
  band = "boot", boot = 20, seed = 5))
say("P8 NULL estimate %s | ~(1|g) estimate %s", f3(a$x$estimate__),
    f3(bf$x$estimate__))
if (inherits(b, "error")) {
  say("P8 reuse of the NULL boot for ~(1|g): REFUSED %s",
      substr(conditionMessage(b), 1, 120))
} else {
  say("P8 reuse of the NULL boot for ~(1|g): ACCEPTED; band %s..%s vs its own boot %s..%s",
      f3(b$x$lower__), f3(b$x$upper__), f3(bf$x$lower__), f3(bf$x$upper__))
  say("P8 reused band contains its estimate: %s",
      all(b$x$lower__ <= b$x$estimate__ & b$x$estimate__ <= b$x$upper__))
}
## factor level order and numeric type in the grid
d9 <- d8
d9$g2 <- as.integer(d9$g)
f9 <- frm(bf(y ~ x + (1 | g2)), family = gaussian(), data = d9)
a9 <- suppressMessages(conditional_effects(f9, "x", resolution = 3,
                                           re_formula = NULL,
                                           conditions = list(g2 = 3L),
                                           band = "boot", boot = 10, seed = 5))
b9 <- tryCatch(suppressMessages(conditional_effects(
  f9, "x", resolution = 3, re_formula = NULL, conditions = list(g2 = 3),
  band = "boot", boot = attr(a9, "boot"))), error = function(e) e)
say("P8b integer 3L boot reused for double 3: %s",
    if (inherits(b9, "error")) "REFUSED" else "ACCEPTED")
say("done")
