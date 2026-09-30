# Reviewer: adversarial shapes for the observed-group fix, fit method.
#   Rscript dev/postfit2-rev-adv-fit.R <base|lane>
# Data seeds and boot seeds are printed beside each number.
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "lane"
libs <- c("C:/Users/adf44/source/r/wt-postfit2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressMessages(library(frmtmb))
say <- function(...) cat(sprintf(...), "\n", sep = "")
say("ARM %s: frmtmb %s from %s", arm, format(packageVersion("frmtmb")),
    find.package("frmtmb"))
f3 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
try_ce <- function(label, expr) {
  r <- tryCatch(suppressMessages(expr), error = function(e) e,
                warning = function(w) w)
  if (inherits(r, "condition")) {
    say("%s: %s: %s", label, class(r)[1],
        substr(gsub("\n", " ", conditionMessage(r)), 1, 220))
    return(invisible(NULL))
  }
  r
}
show <- function(label, ce, key = 1) {
  if (is.null(ce)) return(invisible(NULL))
  d <- ce[[key]]
  say("%s: est %s | lo %s | hi %s", label, f3(d$estimate__), f3(d$lower__),
      f3(d$upper__))
  invisible(d)
}

## A. y ~ x + (1 | g), group SD 2, data seed 21 (the worker's second fit)
set.seed(21)
d <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
d$y <- rnorm(80, 1 + 0.5 * d$x + rnorm(8, 0, 2)[d$g], 1)
fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
re <- ranef(fit)$g
re <- if (length(dim(re)) == 3L) re[, 1, 1] else re[, 1]
say("A modes: %s", f3(re))
pop <- show("A pop wald", conditional_effects(fit, "x", resolution = 3))
say("A1: a level conditions sets that the fit never saw ('99')")
for (spec in list(list(nm = "list(g='99')", c = list(g = "99")),
                  list(nm = "df(g=factor('99'))",
                       c = data.frame(g = factor("99"))))) {
  for (anl in c(FALSE, TRUE)) {
    lab <- sprintf("A1 %s anl=%s", spec$nm, anl)
    show(paste(lab, "wald"), try_ce(lab, conditional_effects(
      fit, "x", resolution = 3, re_formula = NULL, conditions = spec$c,
      allow_new_levels = anl)))
    show(paste(lab, "boot40 seed5"), try_ce(lab, conditional_effects(
      fit, "x", resolution = 3, re_formula = NULL, conditions = spec$c,
      allow_new_levels = anl, band = "boot", boot = 40, seed = 5)))
  }
}
say("A0: no conditions (a new group), for the width a new level has")
show("A0 wald", conditional_effects(fit, "x", resolution = 3,
                                    re_formula = NULL))
show("A0 boot40 seed5", suppressMessages(conditional_effects(
  fit, "x", resolution = 3, re_formula = NULL, band = "boot", boot = 40,
  seed = 5)))
say("A2: rows mixing a level and NA")
cm <- data.frame(g = factor(c("2", NA), levels = 1:8))
rownames(cm) <- c("lev2", "unset")
r <- try_ce("A2 wald", conditional_effects(fit, "x", resolution = 3,
                                           re_formula = NULL,
                                           conditions = cm))
if (!is.null(r)) {
  dd <- r$x
  for (cv in unique(as.character(dd$cond__))) {
    s <- dd[dd$cond__ == cv, ]
    say("A2 wald cond %s: est %s | lo %s | hi %s", cv, f3(s$estimate__),
        f3(s$lower__), f3(s$upper__))
  }
}
r <- try_ce("A2 boot", conditional_effects(fit, "x", resolution = 3,
                                           re_formula = NULL,
                                           conditions = cm, band = "boot",
                                           boot = 40, seed = 5))
if (!is.null(r)) {
  dd <- r$x
  for (cv in unique(as.character(dd$cond__))) {
    s <- dd[dd$cond__ == cv, ]
    say("A2 boot cond %s: est %s | lo %s | hi %s", cv, f3(s$estimate__),
        f3(s$lower__), f3(s$upper__))
  }
}
show("A2 ref: lev2 alone wald", conditional_effects(
  fit, "x", resolution = 3, re_formula = NULL,
  conditions = data.frame(g = factor("2", levels = 1:8))))

say("A3: character grouping column")
dc <- d
dc$g <- as.character(dc$g)
fitc <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dc)
show("A3 chr g='2' wald", try_ce("A3 wald", conditional_effects(
  fitc, "x", resolution = 3, re_formula = NULL,
  conditions = list(g = "2"))))
show("A3 chr g='2' boot", try_ce("A3 boot", conditional_effects(
  fitc, "x", resolution = 3, re_formula = NULL, conditions = list(g = "2"),
  band = "boot", boot = 40, seed = 5)))
show("A3 fct g='2' boot", conditional_effects(
  fit, "x", resolution = 3, re_formula = NULL,
  conditions = list(g = factor("2", levels = 1:8)), band = "boot",
  boot = 40, seed = 5))

## B. two grouping factors, only one set
set.seed(22)
db <- data.frame(x = rnorm(160), g = factor(rep(1:8, 20)),
                 h = factor(rep(1:10, each = 16)))
db$y <- rnorm(160, 1 + 0.5 * db$x + rnorm(8, 0, 2)[db$g] +
                rnorm(10, 0, 0.5)[db$h], 1)
fb <- frm(bf(y ~ x + (1 | g) + (1 | h)), family = gaussian(), data = db)
rg <- ranef(fb)$g
rg <- if (length(dim(rg)) == 3L) rg[, 1, 1] else rg[, 1]
lv <- which.max(abs(rg))
say("B modes g: %s; level %d", f3(rg), lv)
cb <- data.frame(g = factor(lv, levels = 1:8))
w <- show("B g set, h unset: wald", conditional_effects(
  fb, "x", resolution = 3, re_formula = NULL, conditions = cb))
bo <- suppressMessages(conditional_effects(
  fb, "x", resolution = 3, re_formula = NULL, conditions = cb,
  band = "boot", boot = 60, seed = 5))
b <- show("B g set, h unset: boot60 seed5", bo)
m <- colMeans(attr(bo, "boot")$t, na.rm = TRUE)
say("B boot mean minus wald estimate: %s; over |mode| %.3f",
    f3(m - w$estimate__), max(abs(m - w$estimate__)) / abs(rg[lv]))
say("B boot width / wald width: %s",
    f3((b$upper__ - b$lower__) / (w$upper__ - w$lower__)))
cb2 <- data.frame(g = factor(lv, levels = 1:8), h = factor(1, levels = 1:10))
w2 <- show("B both set: wald", conditional_effects(
  fb, "x", resolution = 3, re_formula = NULL, conditions = cb2))
bo2 <- suppressMessages(conditional_effects(
  fb, "x", resolution = 3, re_formula = NULL, conditions = cb2,
  band = "boot", boot = 60, seed = 5))
b2 <- show("B both set: boot60 seed5", bo2)
m2 <- colMeans(attr(bo2, "boot")$t, na.rm = TRUE)
say("B both set: boot mean minus wald estimate: %s; over |mode| %.3f",
    f3(m2 - w2$estimate__), max(abs(m2 - w2$estimate__)) / abs(rg[lv]))
say("B both set: boot width / wald width: %s",
    f3((b2$upper__ - b2$lower__) / (w2$upper__ - w2$lower__)))

## C. nested (1 | g / h): g set, g:h not
set.seed(23)
dn <- data.frame(x = rnorm(160), g = factor(rep(1:8, 20)),
                 h = factor(rep(1:4, each = 2, length.out = 160)))
dn$y <- rnorm(160, 1 + 0.5 * dn$x + rnorm(8, 0, 2)[dn$g] +
                rnorm(32, 0, 0.7)[interaction(dn$g, dn$h)], 1)
fn <- frm(bf(y ~ x + (1 | g / h)), family = gaussian(), data = dn)
rgn <- ranef(fn)$g
rgn <- if (length(dim(rgn)) == 3L) rgn[, 1, 1] else rgn[, 1]
lvn <- which.max(abs(rgn))
say("C modes g: %s; level %d", f3(rgn), lvn)
cn <- data.frame(g = factor(lvn, levels = 1:8))
wn <- show("C g set: wald", try_ce("C wald", conditional_effects(
  fn, "x", resolution = 3, re_formula = NULL, conditions = cn)))
bon <- try_ce("C boot", conditional_effects(
  fn, "x", resolution = 3, re_formula = NULL, conditions = cn,
  band = "boot", boot = 60, seed = 5))
bn <- show("C g set: boot60 seed5", bon)
if (!is.null(bon) && !is.null(wn)) {
  mn <- colMeans(attr(bon, "boot")$t, na.rm = TRUE)
  say("C boot mean minus wald estimate: %s; over |mode| %.3f",
      f3(mn - wn$estimate__), max(abs(mn - wn$estimate__)) / abs(rgn[lvn]))
}

## D. no grouping factor at all, a smooth: re_formula = NULL vs NA
set.seed(24)
ds <- data.frame(x = runif(150))
ds$y <- sin(2 * pi * ds$x) + rnorm(150, 0, 0.3)
fs <- frm(bf(y ~ s(x)), family = gaussian(), data = ds)
dna <- show("D s(x) re_formula=NA boot60 seed5", suppressMessages(
  conditional_effects(fs, "x", resolution = 5, band = "boot", boot = 60,
                      seed = 5)))
dnu <- show("D s(x) re_formula=NULL boot60 seed5", suppressMessages(
  conditional_effects(fs, "x", resolution = 5, band = "boot", boot = 60,
                      seed = 5, re_formula = NULL)))
say("D width NULL / NA: %s", f3((dnu$upper__ - dnu$lower__) /
                                  (dna$upper__ - dna$lower__)))
dw <- show("D s(x) wald", conditional_effects(fs, "x", resolution = 5))
say("D width NULL-boot / wald: %s; NA-boot / wald: %s",
    f3((dnu$upper__ - dnu$lower__) / (dw$upper__ - dw$lower__)),
    f3((dna$upper__ - dna$lower__) / (dw$upper__ - dw$lower__)))
bna <- suppressMessages(conditional_effects(fs, "x", resolution = 5,
                                            band = "boot", boot = 10,
                                            seed = 5))
try_ce("D reuse NA boot for NULL call", conditional_effects(
  fs, "x", resolution = 5, band = "boot", boot = attr(bna, "boot"),
  re_formula = NULL))

## E. boot reuse key: a population boot and an observed-level boot on a
## byte-identical grid (level 1 is the reference level)
say("E: reuse")
pb <- suppressMessages(conditional_effects(fit, "x", resolution = 3,
                                           band = "boot", boot = 10,
                                           seed = 5))
c1 <- data.frame(g = factor("1", levels = 1:8))
ob <- suppressMessages(conditional_effects(fit, "x", resolution = 3,
                                           band = "boot", boot = 10,
                                           seed = 5, re_formula = NULL,
                                           conditions = c1))
say("E pop grid g column: %s; obs grid g column: %s",
    paste(pb$x$g, collapse = ","), paste(ob$x$g, collapse = ","))
r1 <- try_ce("E pop boot into observed-level call", conditional_effects(
  fit, "x", resolution = 3, band = "boot", boot = attr(pb, "boot"),
  re_formula = NULL, conditions = c1))
if (!is.null(r1)) say("E pop boot into observed-level call: ACCEPTED")
r2 <- try_ce("E observed boot into pop call (with conditions g=1)",
             conditional_effects(fit, "x", resolution = 3, band = "boot",
                                 boot = attr(ob, "boot"), conditions = c1))
if (!is.null(r2)) say("E observed boot into pop call: ACCEPTED")
r3 <- try_ce("E observed boot into pop call (no conditions)",
             conditional_effects(fit, "x", resolution = 3, band = "boot",
                                 boot = attr(ob, "boot")))
if (!is.null(r3)) say("E observed boot into pop call (no cond): ACCEPTED")
r4 <- try_ce("E observed boot reused by the same call",
             conditional_effects(fit, "x", resolution = 3, band = "boot",
                                 boot = attr(ob, "boot"), re_formula = NULL,
                                 conditions = c1))
if (!is.null(r4)) say("E same-call reuse: ACCEPTED, identical band %s",
                      identical(r4$x$lower__, ob$x$lower__))
nb <- suppressMessages(conditional_effects(fit, "x", resolution = 3,
                                           band = "boot", boot = 10,
                                           seed = 5, re_formula = NULL))
r5 <- try_ce("E new-level boot into observed-level call",
             conditional_effects(fit, "x", resolution = 3, band = "boot",
                                 boot = attr(nb, "boot"), re_formula = NULL,
                                 conditions = c1))
if (!is.null(r5)) say("E new-level boot into observed call: ACCEPTED")
## a formula re_formula with g set
show("E2 re_formula=~(1|g) g=2 boot", try_ce("E2", conditional_effects(
  fit, "x", resolution = 3, band = "boot", boot = 40, seed = 5,
  re_formula = ~ (1 | g),
  conditions = data.frame(g = factor("2", levels = 1:8)))))
show("E2 re_formula=NULL g=2 boot", conditional_effects(
  fit, "x", resolution = 3, band = "boot", boot = 40, seed = 5,
  re_formula = NULL,
  conditions = data.frame(g = factor("2", levels = 1:8))))
say("done")
