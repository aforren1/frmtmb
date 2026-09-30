# Reviewer probes (lane ceplot):
#  1. crossed unseen g:h with CHARACTER grouping columns (fit, boot, draws)
#  2. y ~ x + trt + (1 | trt:subj) with nothing set: the NA label path
#     (base refuses boot and draws by name; see pf2/p1-attack-*.txt)
#  3. the "Valid effects" list on a nonlinear model
#  4. an invalid effect with data = (the NEWS "flat curve" claim)
#   Rscript dev/ceplot-rev-misc.R lane|base
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/wt-ceplot-lib",
          "C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
source("C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-rev-shapes.R")
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
f4 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
run <- function(label, expr) {
  w <- character(0)
  r <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }), error = function(e) paste("ERROR", conditionMessage(e)))
  if (is.character(r)) {
    cat(label, ":", substr(r, 1, 220), "\n")
  } else {
    cat(label, ": est", f4(r$estimate__), "| width",
        f4(r$upper__ - r$lower__), "\n")
  }
  for (x in unique(w)) cat("   warning:", substr(x, 1, 200), "\n")
}
## 1
dA <- crossed_data(49)
dch <- dA
dch$g <- as.character(dch$g)
dch$h <- as.character(dch$h)
fch <- frm(bf(y ~ x + (1 | g) + (1 | h) + (1 | g:h)), family = gaussian(),
           data = dch)
ffa <- frm(bf(y ~ x + (1 | g) + (1 | h) + (1 | g:h)), family = gaussian(),
           data = dA)
cat("logLik character", logLik(fch), "factor", logLik(ffa), "\n")
ce1 <- function(o, cond, ...) {
  conditional_effects(o, "x", resolution = 3, re_formula = NULL,
                      conditions = cond, ...)$x
}
for (cc in list(list(g = "1", h = "1"), list(g = "1", h = "2"))) {
  lab <- paste0("g=", cc$g, ",h=", cc$h)
  run(paste("1 character", lab, "wald"), ce1(fch, cc))
  run(paste("1 factor   ", lab, "wald"), ce1(ffa, cc))
  run(paste("1 character", lab, "boot40"), ce1(fch, cc, band = "boot",
                                               boot = 40, seed = 5))
  run(paste("1 factor   ", lab, "boot40"), ce1(ffa, cc, band = "boot",
                                               boot = 40, seed = 5))
  run(paste("1 character", lab, "draws"), ce1(hand(fch, 100), cc, seed = 1))
  run(paste("1 factor   ", lab, "draws"), ce1(hand(ffa, 100), cc, seed = 1))
}
## 2
set.seed(41)
dt <- expand.grid(subj = factor(1:12), trt = factor(c("a", "b")), r = 1:5)
dt$x <- rnorm(nrow(dt))
dt$y <- rnorm(nrow(dt), 1 + 0.5 * dt$x + (dt$trt == "b") +
                rnorm(24)[as.integer(interaction(dt$trt, dt$subj))], 0.5)
ft <- frm(bf(y ~ x + trt + (1 | trt:subj)), family = gaussian(), data = dt)
run("2 trt:subj nothing set, wald", ce1(ft, list()))
run("2 trt:subj nothing set, boot", ce1(ft, list(), band = "boot", boot = 10,
                                        seed = 1))
run("2 trt:subj nothing set, draws", ce1(hand(ft, 50), list(), seed = 1))
run("2 trt:subj trt = b, subj = 99, boot",
    ce1(ft, list(trt = "b", subj = "99"), band = "boot", boot = 10, seed = 1))
## 3
set.seed(11)
dn <- data.frame(x = rnorm(200), z = rnorm(200), extra = rnorm(200))
dn$yp <- 2 * exp(0.3 * dn$x) + 0.2 * dn$z + rnorm(200, 0, 0.3)
fn <- frm(bf(yp ~ a * exp(b * x), a ~ 1 + z, b ~ 1, nl = TRUE),
          family = gaussian(), data = dn)
w <- tryCatch(withCallingHandlers(
  conditional_effects(fn, effects = c("x", "extra"), resolution = 3,
                      data = dn),
  warning = function(x) {
    cat("3 warning:", conditionMessage(x), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) cat("3 ERROR", conditionMessage(e), "\n"))
## 4
set.seed(12)
d4 <- data.frame(x = rnorm(100), extra = rnorm(100))
d4$y <- rnorm(100, 0.5 * d4$x)
f4fit <- frm(bf(y ~ x), family = gaussian(), data = d4)
r4 <- tryCatch(withCallingHandlers(
  conditional_effects(f4fit, effects = "extra", resolution = 3, data = d4),
  warning = function(x) invokeRestart("muffleWarning")),
  error = function(e) paste("ERROR", conditionMessage(e)))
if (is.character(r4)) cat("4 effects = 'extra' with data =:", substr(r4, 1, 150), "\n") else
  cat("4 effects = 'extra' with data =: est", f4(r4$extra$estimate__), "\n")
