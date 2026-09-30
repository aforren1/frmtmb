# Lane ceplot: what the base build (0.66.0, rellib-r4) does on each item.
#   Rscript dev/ceplot-probe-base.R [base|lane] > dev/ceplot-log/probe-<arm>.txt
arm <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(arm)) arm <- "base"
libs <- c("C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ceplot-lib", libs)
.libPaths(libs)
suppressMessages(library(frmtmb))
cat("ARM", arm, "frmtmb from", find.package("frmtmb"), "\n")
grDevices::pdf(NULL)
try_ <- function(label, expr) {
  w <- character()
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  if (inherits(r, "error")) {
    cat(label, ": ERROR ", conditionMessage(r), "\n", sep = "")
  } else {
    cat(label, ": OK class ", paste(class(r), collapse = "/"), "\n", sep = "")
  }
  for (x in w) cat("   warning: ", x, "\n", sep = "")
  invisible(r)
}

set.seed(1)
d <- data.frame(x = rnorm(100), f = factor(rep(c("a", "b"), 50)))
d$y <- rnorm(100, 1 + 0.5 * d$x + (d$f == "b"))
fit <- frm(bf(y ~ x * f), family = gaussian(), data = d)
ce <- conditional_effects(fit, effects = c("x", "f"))

## item 1: plot() of conditional effects, brms's arguments
try_("1 plot(rug, plot = FALSE)",
     plot(ce, points = TRUE, rug = TRUE, ask = FALSE, plot = FALSE))
try_("1 plot(stype, plot = FALSE)", plot(ce, stype = "raster", plot = FALSE))
try_("1 plot(theme = NULL)", plot(ce, theme = NULL, ask = FALSE))
try_("1 plot(mean = FALSE)", plot(ce, mean = FALSE, ask = FALSE))

## item 2: plot() of a hypothesis
h <- hypothesis(fit, "x > 0")
try_("2 plot(hyp, plot = FALSE)", plot(h, plot = FALSE))
try_("2 plot(hyp, ignore_prior, plot = FALSE)",
     plot(h, ignore_prior = TRUE, plot = FALSE))
try_("2 plot(hyp, N = 5, ask = FALSE)", plot(h, N = 5, ask = FALSE))

## item 3a: one invalid effect among valid ones
try_("3a effects = c(xx, x)", conditional_effects(fit, effects = c("xx", "x")))
try_("3a effects = xx", conditional_effects(fit, effects = "xx"))

## item 3b: ordinal
set.seed(2)
do <- data.frame(x = rnorm(200))
do$y <- factor(cut(do$x + rlogis(200), c(-Inf, -1, 0, 1, Inf)),
               ordered = TRUE)
fo <- frm(bf(y ~ x), family = cumulative(), data = do)
try_("3b ordinal default", conditional_effects(fo))
try_("3b ordinal categorical = FALSE",
     conditional_effects(fo, categorical = FALSE))
try_("3b ordinal dpar = mu", conditional_effects(fo, dpar = "mu"))

## item 3c: mi(x, idx = )
set.seed(26)
n <- 120
dm <- data.frame(g1 = sample(seq(1, n - 1, 2), n, TRUE), g2 = seq_len(n),
                 s = rep(c(TRUE, FALSE), n / 2), w = rnorm(n))
dm$x <- rnorm(n)
dm$y <- 1 + 0.5 * dm$x[match(dm$g1, dm$g2)] + 0.3 * dm$w +
  rnorm(n, sd = 0.5)
mv <- frm(bf(y ~ mi(x, idx = g1) + w) +
            bf(x | mi() + index(g2) + subset(s) ~ 1),
          data = dm, family = gaussian())
try_("3c mi idx, default", conditional_effects(mv, resp = "y"))
try_("3c mi idx, w", conditional_effects(mv, "w", resp = "y"))

## item 3d: crossed, unseen g:h
set.seed(49)
dc <- expand.grid(g = factor(1:6), h = factor(1:5), r = 1:4)
dc <- dc[!(dc$g == "1" & dc$h == "1"), ]
dc$x <- rnorm(nrow(dc))
dc$y <- rnorm(nrow(dc), 1 + 0.5 * dc$x + rnorm(6)[dc$g] + rnorm(5)[dc$h] +
                rnorm(30, 0, 0.7)[as.integer(interaction(dc$g, dc$h))], 0.5)
fc <- frm(bf(y ~ x + (1 | g) + (1 | h) + (1 | g:h)), family = gaussian(),
          data = dc)
try_("3d crossed wald", conditional_effects(fc, "x", resolution = 3,
     re_formula = NULL, conditions = list(g = "1", h = "1")))
try_("3d crossed boot", conditional_effects(fc, "x", resolution = 3,
     re_formula = NULL, conditions = list(g = "1", h = "1"), band = "boot",
     boot = 5, seed = 1))

## item 4
set.seed(3)
dg <- data.frame(x = rnorm(120), g = factor(rep(1:12, 10)))
dg$y <- rnorm(120, 1 + 0.5 * dg$x + rnorm(12)[dg$g])
fg <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dg)
nd <- data.frame(x = c(0, 1), g = factor(c("99", "1")))
try_("4a fitted old_levels", fitted(fg, newdata = nd, allow_new_levels = TRUE,
                                    sample_new_levels = "old_levels"))
try_("4a predict old_levels", predict(fg, newdata = nd,
                                      allow_new_levels = TRUE,
                                      sample_new_levels = "old_levels",
                                      ndraws = 20))
try_("4c parnames(fit)", parnames(fg))
cat("done\n")
