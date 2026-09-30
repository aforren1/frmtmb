# Reviewer: adversarial shapes for the observed-group fix, draws method.
# Hand-built draws (the lane's construction in test-postfit-draws.R):
# the ML estimate plus N(0, 0.05) noise, 200 draws, seed 1; no sampler.
#   Rscript dev/postfit2-rev-adv-draws.R <base|lane>
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "lane"
libs <- c("C:/Users/adf44/source/r/wt-postfit2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
say <- function(...) cat(sprintf(...), "\n", sep = "")
say("ARM %s: frmtmb %s, frmtmb.sample %s from %s", arm,
    format(packageVersion("frmtmb")), format(packageVersion("frmtmb.sample")),
    find.package("frmtmb.sample"))
f3 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
hand_draws <- function(fit, n = 200, sd = 0.05, seed = 1) {
  tpl <- fit$frame$par_template
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  set.seed(seed)
  M <- matrix(rep(est, each = n) + rnorm(n * length(est), 0, sd), n,
              dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  structure(list(stanfit = NULL, draws = M, fit = fit),
            class = "frmtmb_draws")
}
run <- function(label, expr) {
  r <- tryCatch(suppressMessages(expr), error = function(e) e)
  if (inherits(r, "condition")) {
    say("%s: %s: %s", label, class(r)[1],
        substr(gsub("\n", " ", conditionMessage(r)), 1, 200))
    return(invisible(NULL))
  }
  for (cv in unique(as.character(r[[1]]$cond__))) {
    s <- r[[1]][as.character(r[[1]]$cond__) == cv, ]
    say("%s [cond %s]: est %s | lo %s | hi %s | width %s", label, cv,
        f3(s$estimate__), f3(s$lower__), f3(s$upper__),
        f3(s$upper__ - s$lower__))
  }
  invisible(r)
}
set.seed(9)
dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
ds <- hand_draws(fit)
say("sd_g at the estimate %.4f", exp(fit$estimates$theta[1]))
run("pop", conditional_effects(ds, "x", resolution = 3, seed = 1))
run("new group (no conditions)", conditional_effects(
  ds, "x", resolution = 3, re_formula = NULL, seed = 1))
run("observed g=2", conditional_effects(
  ds, "x", resolution = 3, re_formula = NULL, seed = 1,
  conditions = data.frame(g = factor(2, levels = 1:6))))
for (anl in c(FALSE, TRUE)) {
  run(sprintf("unseen g='99' anl=%s", anl), conditional_effects(
    ds, "x", resolution = 3, re_formula = NULL, seed = 1,
    conditions = list(g = "99"), allow_new_levels = anl))
}
cm <- data.frame(g = factor(c("2", NA), levels = 1:6))
rownames(cm) <- c("lev2", "unset")
run("mixed rows g=2 / NA", conditional_effects(
  ds, "x", resolution = 3, re_formula = NULL, seed = 1, conditions = cm))
run("mixed rows g=2 / NA, anl", conditional_effects(
  ds, "x", resolution = 3, re_formula = NULL, seed = 1, conditions = cm,
  allow_new_levels = TRUE))

## crossed: g set, h unset
set.seed(22)
db <- data.frame(x = rnorm(160), g = factor(rep(1:8, 20)),
                 h = factor(rep(1:10, each = 16)))
db$y <- rnorm(160, 1 + 0.5 * db$x + rnorm(8, 0, 2)[db$g] +
                rnorm(10, 0, 0.5)[db$h], 1)
fb <- frm(bf(y ~ x + (1 | g) + (1 | h)), family = gaussian(), data = db)
dsb <- hand_draws(fb)
run("crossed g=7 set, h unset", conditional_effects(
  dsb, "x", resolution = 3, re_formula = NULL, seed = 1,
  conditions = data.frame(g = factor(7, levels = 1:8))))
run("crossed g=7, h=1 both set", conditional_effects(
  dsb, "x", resolution = 3, re_formula = NULL, seed = 1,
  conditions = data.frame(g = factor(7, levels = 1:8),
                          h = factor(1, levels = 1:10))))
run("crossed neither set", conditional_effects(
  dsb, "x", resolution = 3, re_formula = NULL, seed = 1))
## a moved draw column shows which level's draws the curve reads
mv <- dsb
mv$draws[, "r_g[7,Intercept]"] <- mv$draws[, "r_g[7,Intercept]"] + 10
a0 <- suppressMessages(conditional_effects(
  dsb, "x", resolution = 3, re_formula = NULL, seed = 1,
  conditions = data.frame(g = factor(7, levels = 1:8))))[[1]]
a1 <- suppressMessages(conditional_effects(
  mv, "x", resolution = 3, re_formula = NULL, seed = 1,
  conditions = data.frame(g = factor(7, levels = 1:8))))[[1]]
say("crossed g=7 set: curve moves by %s when r_g[7,] moves by 10",
    f3(a1$estimate__ - a0$estimate__))

## nested (1 | g / h): g set
set.seed(23)
dn <- data.frame(x = rnorm(160), g = factor(rep(1:8, 20)),
                 h = factor(rep(1:4, each = 2, length.out = 160)))
dn$y <- rnorm(160, 1 + 0.5 * dn$x + rnorm(8, 0, 2)[dn$g] +
                rnorm(32, 0, 0.7)[interaction(dn$g, dn$h)], 1)
fn <- frm(bf(y ~ x + (1 | g / h)), family = gaussian(), data = dn)
dsn <- hand_draws(fn)
run("nested g=5 set", conditional_effects(
  dsn, "x", resolution = 3, re_formula = NULL, seed = 1,
  conditions = data.frame(g = factor(5, levels = 1:8))))
run("nested g=5, h=1 set", conditional_effects(
  dsn, "x", resolution = 3, re_formula = NULL, seed = 1,
  conditions = data.frame(g = factor(5, levels = 1:8),
                          h = factor(1, levels = 1:4))))
run("nested neither set", conditional_effects(
  dsn, "x", resolution = 3, re_formula = NULL, seed = 1))
say("done")
