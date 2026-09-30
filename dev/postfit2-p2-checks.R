# Punch round 2: the numbers of dev/postfit2-findings.md section 7d
# that are not brms comparisons (those are dev/postfit2-p2-brms.R).
#   Rscript dev/postfit2-p2-checks.R <lane|r1> > dev/postfit2-log/p2-checks-<arm>.txt
# Seeds: data 44 (gr by), 43 (mm), 41 (trt:subj), 48 ((1 + x | g)),
# 47 (fs), 49 (crossed); boot seed 5, draws seed 1.
arm <- commandArgs(TRUE)[1]
lib1 <- if (identical(arm, "r1")) "C:/Users/adf44/source/r/wt-postfit2-r1lib" else
  "C:/Users/adf44/source/r/wt-postfit2-lib"
.libPaths(c(lib1, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
say <- function(...) cat(sprintf(...), "\n", sep = "")
f4 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
wd <- function(d) d$upper__ - d$lower__
hand <- function(fit, n = 200, seed = 1) {
  tpl <- fit$frame$par_template
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  set.seed(seed)
  M <- matrix(rep(est, each = n) + rnorm(n * length(est), 0, 0.05), n,
              dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
  structure(list(stanfit = NULL,
                 draws = cbind(frmtmb.sample:::draws_to_natural(M, fit),
                               lp__ = 0),
                 fit = fit), class = "frmtmb_draws")
}
run <- function(label, expr) {
  r <- tryCatch(suppressMessages(expr), error = function(e) e)
  if (inherits(r, "error")) {
    return(say("%s: ERROR %s", label, substr(conditionMessage(r), 1, 160)))
  }
  x <- r[[1]]
  say("%s: est %s | width %s", label, f4(x$estimate__), f4(wd(x)))
  invisible(x)
}
ce <- function(o, ...) {
  conditional_effects(o, "x", resolution = 3, re_formula = NULL, ...)
}
say("ARM %s: frmtmb from %s", arm, find.package("frmtmb"))

## P1-B1: gr(g, by = f)
set.seed(44)
d4 <- data.frame(x = rnorm(240), g = factor(rep(1:12, 20)))
d4$f <- factor(ifelse(as.integer(d4$g) <= 6, "a", "b"))
d4$y <- rnorm(240, 1 + 0.5 * d4$x +
                rnorm(12, 0, ifelse(1:12 <= 6, 0.3, 1.5))[d4$g], 0.5)
f4m <- frm(bf(y ~ x + f + (1 | gr(g, by = f))), family = gaussian(),
           data = d4)
ds4 <- hand(f4m)
for (cnd in list(list(f = "a", g = "3"), list(f = "b", g = "9"),
                 list(f = "a"), list(f = "b"), list(f = "a", g = "99"))) {
  lab <- paste0("B1 ", paste(names(cnd), unlist(cnd), sep = " = ",
                             collapse = ", "))
  run(paste(lab, "| wald"), ce(f4m, conditions = cnd))
  run(paste(lab, "| boot40"), ce(f4m, conditions = cnd, band = "boot",
                                 boot = 40, seed = 5))
  run(paste(lab, "| draws"), ce(ds4, conditions = cnd, seed = 1))
}

## P1-B2: y ~ x + trt + (1 | trt:subj), trt = "b"
set.seed(41)
d1 <- expand.grid(rep = 1:4, trt = factor(c("a", "b")), subj = factor(1:12))
d1$x <- rnorm(nrow(d1))
d1$y <- rnorm(nrow(d1), 1 + 0.5 * d1$x + (d1$trt == "b") +
                rnorm(24, 0, 1)[interaction(d1$trt, d1$subj)], 0.5)
f1 <- frm(bf(y ~ x + trt + (1 | trt:subj)), family = gaussian(), data = d1)
b <- suppressMessages(ce(f1, conditions = list(trt = "b"), band = "boot",
                         boot = 40, seed = 5))
say("B2 trt = b: wald est %s | boot mean of refits %s",
    f4(suppressMessages(ce(f1, conditions = list(trt = "b")))$x$estimate__),
    f4(colMeans(attr(b, "boot")$t, na.rm = TRUE)))

## P1-M1: mm(g1, g2)
set.seed(43)
d3 <- data.frame(x = rnorm(200), g1 = factor(sample(1:10, 200, TRUE)),
                 g2 = factor(sample(1:10, 200, TRUE)))
u <- rnorm(10, 0, 1)
d3$y <- rnorm(200, 1 + 0.5 * d3$x + 0.5 * (u[d3$g1] + u[d3$g2]), 0.5)
f3 <- frm(bf(y ~ x + (1 | mm(g1, g2))), family = gaussian(), data = d3)
ds3 <- hand(f3)
say("M1 VarCorr sd %s", f4(VarCorr(f3)[[1]]$sd[, 1]))
x <- run("M1 population (re_formula NA) | wald",
         conditional_effects(f3, "x", resolution = 3))
x <- suppressMessages(ce(f3, conditions = list()))$x
say("M1 nothing set: frame g1 %s g2 %s", as.character(x$g1[1]),
    as.character(x$g2[1]))
for (cnd in list(list(), list(g1 = "2", g2 = "2"), list(g1 = "2"),
                 list(g1 = "99", g2 = "98"))) {
  lab <- paste0("M1 ", if (length(cnd)) paste(names(cnd), unlist(cnd),
                                              sep = " = ", collapse = ", ")
                else "nothing set")
  w <- run(paste(lab, "| wald"), ce(f3, conditions = cnd))
  bo <- run(paste(lab, "| boot60"), ce(f3, conditions = cnd, band = "boot",
                                       boot = 60, seed = 5))
  if (!is.null(w) && !is.null(bo)) {
    say("%s: boot width / wald width %s", lab, f4(wd(bo) / wd(w)))
  }
  run(paste(lab, "| draws"), ce(ds3, conditions = cnd, seed = 1))
}

## P1-M3: the reuse key carries re_formula
set.seed(48)
d8 <- data.frame(x = rnorm(240), g = factor(rep(1:12, 20)))
b0 <- rnorm(12, 0, 1)
b1 <- rnorm(12, 0, 1)
d8$y <- rnorm(240, 1 + b0[d8$g] + (0.5 + b1[d8$g]) * d8$x, 0.5)
f8 <- frm(bf(y ~ x + (1 + x | g)), family = gaussian(), data = d8)
a <- suppressMessages(ce(f8, conditions = list(g = "3"), band = "boot",
                         boot = 20, seed = 5))
run("M3 re_formula NULL | boot20", list(a$x))
run("M3 ~ (1 | g) reusing the NULL bootstrap",
    conditional_effects(f8, "x", resolution = 3, re_formula = ~ (1 | g),
                        conditions = list(g = "3"), band = "boot",
                        boot = attr(a, "boot")))
run("M3 ~ (1 | g) | wald",
    conditional_effects(f8, "x", resolution = 3, re_formula = ~ (1 | g),
                        conditions = list(g = "3")))

## P1-M2: crossed (1 | g) + (1 | h) + (1 | g:h), an unseen g:h
set.seed(49)
dc <- expand.grid(g = factor(1:6), h = factor(1:5), r = 1:4)
dc <- dc[!(dc$g == "1" & dc$h == "1"), ]
dc$x <- rnorm(nrow(dc))
dc$y <- rnorm(nrow(dc), 1 + 0.5 * dc$x + rnorm(6)[dc$g] + rnorm(5)[dc$h] +
                rnorm(30, 0, 0.7)[as.integer(interaction(dc$g, dc$h))], 0.5)
fc <- frm(bf(y ~ x + (1 | g) + (1 | h) + (1 | g:h)), family = gaussian(),
          data = dc)
dsc <- hand(fc)
cc <- list(g = "1", h = "1")
run("M2 g = 1, h = 1 (never together) | wald", ce(fc, conditions = cc))
run("M2 g = 1, h = 2 (observed together) | wald",
    ce(fc, conditions = list(g = "1", h = "2")))
run("M2 g = 1, h = 1 | boot20", ce(fc, conditions = cc, band = "boot",
                                   boot = 20, seed = 5))
run("M2 g = 1, h = 1 | draws", ce(dsc, conditions = cc, seed = 1))

## minor: fs smooth at an unseen level
set.seed(47)
d7 <- data.frame(x = runif(200), g = factor(rep(1:5, 40)))
d7$y <- sin(2 * pi * d7$x) + rnorm(5, 0, 0.5)[d7$g] + rnorm(200, 0, 0.3)
f7 <- frm(bf(y ~ s(x, g, bs = "fs", k = 5)), family = gaussian(), data = d7)
run("fs g = 99, re_formula NULL | wald", ce(f7, conditions = list(g = "99")))
run("fs g = 99, re_formula NA | wald",
    conditional_effects(f7, "x", resolution = 3, conditions = list(g = "99")))
say("done")
