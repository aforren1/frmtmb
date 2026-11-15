# Short check of parity2's fix round: F1 (old_levels on a fit), F2, F3
# names on the frmtmb side. Arg "base" drops the lane library.
#   Rscript dev/reviews/parity2-short-a.R [base] > <out> 2>&1
args <- commandArgs(TRUE)
base <- length(args) && identical(args[1], "base")
.libPaths(c(if (!base) "C:/Users/adf44/source/r/rellib-r22" else
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
try_msg <- function(expr) {
  r <- tryCatch(expr, error = function(e) e)
  if (inherits(r, "error")) paste("ERROR:", substr(conditionMessage(r), 1, 200))
  else r
}

cat("\n#### F1/F2 data: 30 groups of 2, sd_g 1, sigma 2, seed 100\n")
set.seed(100)
d <- data.frame(x = rnorm(60), g = factor(rep(1:30, each = 2)))
d$y <- 1 + 0.5 * d$x + rnorm(30)[d$g] + rnorm(60, sd = 2)
fb <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
ic <- list(x = c(-1, 0, 1))
wid <- function(ce) ce$x$upper__ - ce$x$lower__

cat("\n## sample_new_levels on a fit's CE, each band\n")
for (s in c("gaussian", "uncertainty", "old_levels")) {
  for (b in c("wald", "boot", "profile")) {
    r <- try_msg(suppressMessages(conditional_effects(
      fb, "x", re_formula = NULL, band = b, boot = 4, seed = 2,
      int_conditions = ic, sample_new_levels = s)))
    cat(s, b, ":", if (is.character(r)) r else
      paste("answered, width", paste(format(wid(r), digits = 4),
                                     collapse = " ")), "\n")
  }
}
cat("\n## re_formula = NA with sample_new_levels (brms ignores it)\n")
for (s in c("gaussian", "uncertainty")) {
  r <- try_msg(conditional_effects(fb, "x", re_formula = NA,
                                   int_conditions = ic,
                                   sample_new_levels = s))
  cat(s, "NA wald:", if (is.character(r)) r else "answered", "\n")
}

if (!base) {
  cw <- conditional_effects(fb, "x", re_formula = NULL, int_conditions = ic)
  cg <- conditional_effects(fb, "x", re_formula = NULL, int_conditions = ic,
                            sample_new_levels = "gaussian")
  cat("\nF2 gaussian on Wald identical to default:", identical(cw, cg), "\n")
  bd <- suppressMessages(conditional_effects(
    fb, "x", re_formula = NULL, band = "boot", boot = 20, seed = 2,
    int_conditions = ic))
  bg <- suppressMessages(conditional_effects(
    fb, "x", re_formula = NULL, band = "boot", boot = 20, seed = 2,
    int_conditions = ic, sample_new_levels = "gaussian"))
  cat("F1 gaussian on boot identical to default:", identical(bd, bg), "\n")
}

cat("\n## fitted() at a new level, old_levels vs gaussian vs observed\n")
nd <- data.frame(x = c(-1, 0, 1), g = "new")
fg <- try_msg(fitted(fb, newdata = nd, allow_new_levels = TRUE))
print(fg)
set.seed(7)
fo <- try_msg(fitted(fb, newdata = nd, allow_new_levels = TRUE,
                     sample_new_levels = "old_levels"))
print(fo)
obs <- do.call(rbind, lapply(levels(d$g), function(l) {
  f <- fitted(fb, newdata = data.frame(x = c(-1, 0, 1), g = l))
  data.frame(g = l, x = c(-1, 0, 1), f[, 1:4])
}))
if (!is.character(fo)) {
  obs0 <- obs[obs$x == 0, ]
  hit <- obs0[abs(obs0$Estimate - fo[2, "Estimate"]) < 1e-10, ]
  cat("old_levels x=0 matches observed level:\n"); print(hit)
  l <- hit$g[1]
  cat("observed level", l, "fitted rows:\n"); print(obs[obs$g == l, ])
  cat("old_levels width (Q97.5 - Q2.5):",
      format(fo[, 4] - fo[, 3], digits = 4), "\n")
  cat("observed level width:",
      format(obs$Q97.5[obs$g == l] - obs$Q2.5[obs$g == l], digits = 4), "\n")
  cw_l <- conditional_effects(fb, "x", re_formula = NULL,
                              int_conditions = ic,
                              conditions = data.frame(g = l))
  cat("CE Wald at conditions g =", l, "estimate",
      format(cw_l$x$estimate__, digits = 6), "width",
      format(wid(cw_l), digits = 4), "se", format(cw_l$x$se__, digits = 4),
      "\n")
  cat("fitted Est.Error at that level:",
      format(obs$Est.Error[obs$g == l], digits = 4), "\n")
  cat("old_levels Est.Error:", format(fo[, 2], digits = 4), "\n")
}
cat("observed levels' width at x=0, range:",
    format(range(obs$Q97.5[obs$x == 0] - obs$Q2.5[obs$x == 0]), digits = 4),
    "\n")
cat("observed levels' Est.Error at x=0, range:",
    format(range(obs$Est.Error[obs$x == 0]), digits = 4), "\n")

cat("\n#### F3 default displays, frmtmb\n")
set.seed(21)
n <- 200
e <- data.frame(x = rnorm(n), z = rnorm(n),
                f = factor(sample(c("a", "b", "c"), n, TRUE)),
                h = factor(sample(c("u", "v"), n, TRUE)),
                grp = factor(sample(letters[1:8], n, TRUE)),
                kk = sample(1:6, n, TRUE))
e$ch <- as.character(e$f)
e$lg <- e$h == "u"
e$y <- 1 + e$x + 0.5 * e$z + as.integer(e$f) * 0.3 + rnorm(8)[e$grp] +
  rnorm(6)[e$kk] + rnorm(n)
shapes <- list(
  fx = y ~ f * x, xf = y ~ x * f, fh = y ~ f * h, xz = y ~ x * z,
  fx_re = y ~ f * x + (1 | grp), sby = y ~ s(x, by = f),
  chx = y ~ ch * x, lgx = y ~ lg * x, kkx_grp = y ~ kk * x + (1 | kk))
for (nm in names(shapes)) {
  fit <- try_msg(suppressMessages(frm(bf(shapes[[nm]]),
                                      family = gaussian(), data = e)))
  r <- if (is.character(fit)) fit else try_msg(names(suppressMessages(
    conditional_effects(fit, resolution = 5))))
  cat(nm, ":", paste(r, collapse = " "), "\n")
  if (nm %in% c("fx", "xf", "chx", "lgx") && !is.character(fit)) {
    for (ef in list("f:x", "x:f", c("x:f", "x"), c("x", "f:x"))) {
      ef2 <- if (nm == "chx") sub("f", "ch", ef) else
        if (nm == "lgx") sub("f", "lg", ef) else ef
      r2 <- try_msg(names(suppressMessages(conditional_effects(
        fit, effects = ef2, resolution = 5))))
      cat("  effects =", paste(ef2, collapse = ","), ":",
          paste(r2, collapse = " "), "\n")
    }
    # the display's axis: the first column varies over resolution points
    ce <- suppressMessages(conditional_effects(fit, resolution = 5))
    k <- grep(":", names(ce), value = TRUE)
    if (length(k)) {
      cat("  ", k, "effect1__ class", class(ce[[k]]$effect1__)[1],
          "n distinct", length(unique(ce[[k]]$effect1__)),
          "effect2__ class", class(ce[[k]]$effect2__)[1], "\n")
    }
  }
}
cat("done\n")
