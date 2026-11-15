# Review of lane slopescale, attack items 1 and 5: ordinary fits with
# random slopes at spread near 1, base (rellib-r21) against the lane, in
# ONE process (frmtmb has no compiled code, so the namespace is unloaded
# and reloaded from the other library; rtmb-pitfalls item 21). Arms run
# base, lane, base: the two base runs show the process is deterministic.
# Compares opt$par, logLik, vcov, VarCorr Est.Error, conditions bitwise.
# Usage: Rscript slopescale-review-ordinary.R
LANE <- "C:/Users/adf44/source/r/rellib-r22"
BASE <- "C:/Users/adf44/source/r/rellib-r21"
.libPaths(c(BASE, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
use <- function(lib) {
  if ("frmtmb" %in% loadedNamespaces()) unloadNamespace("frmtmb")
  suppressPackageStartupMessages(library(frmtmb, lib.loc = lib))
  cat("lib:", dirname(getNamespaceInfo("frmtmb", "path")), as.character(packageVersion("frmtmb", lib.loc = lib)), "
")
  cat("loaded", find.package("frmtmb"), "\n")
}
data("sleepstudy", package = "lme4")
mk <- function(seed, fam = "gaussian", sd1 = 0.5, xfun = function(n) rnorm(n)) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  x <- xfun(300)
  tf <- factor(rep(1:3, length.out = 300))
  eta <- 0.3 + rnorm(20, 0, 0.5)[g] + rnorm(20, 0, sd1)[g] * x + 0.2 * x
  y <- switch(fam,
              gaussian = eta + rnorm(300),
              poisson = rpois(300, exp(eta)),
              bernoulli = rbinom(300, 1, plogis(eta)))
  data.frame(g = g, x = x, y = y, tf = tf)
}
fam_of <- function(f) switch(f, gaussian = gaussian(), poisson = poisson(),
                             bernoulli = bernoulli())
cases <- list()
add <- function(nm, fo, d, fam = "gaussian") {
  cases[[nm]] <<- list(fo = fo, d = d, fam = fam)
}
ss <- sleepstudy
add("sleep_DaysDays", Reaction ~ Days + (Days | Subject), ss)
add("sleep_1Days", Reaction ~ 1 + (Days | Subject), ss)
add("sleep_0Days", Reaction ~ Days + (0 + Days | Subject), ss)
for (s in 1:3) add(paste0("gau_1x_", s), y ~ 1 + (1 + x | g), mk(s))
for (s in 1:2) add(paste0("gau_xx_", s), y ~ x + (1 + x | g), mk(s))
for (s in 1:3) add(paste0("pois_1x_", s), y ~ 1 + (1 + x | g),
                   mk(s, "poisson"), "poisson")
for (s in 1:3) add(paste0("bern_xx_", s), y ~ x + (1 + x | g),
                   mk(s, "bernoulli"), "bernoulli")
add("gau_cs", y ~ 1 + cs(1 + x | g), mk(4))
add("gau_ar1", y ~ 1 + (1 | g) + ar1(0 + tf:x | g), mk(5))
add("gau_homdiag", y ~ 1 + homdiag(1 + x | g), mk(6))
add("gau_unif01", y ~ 1 + (1 + x | g), mk(7, xfun = function(n) runif(n)))
add("gau_x5", y ~ 1 + (1 + x | g), mk(8, xfun = function(n) 5 + rnorm(n)))
# singular (true slope sd 0) at spread 1
add("gau_sing", y ~ 1 + (1 + x | g), mk(9, sd1 = 0))
add("pois_sing", y ~ 1 + (0 + x | g), mk(10, "poisson", sd1 = 0), "poisson")

run1 <- function(cs) {
  cond <- character(0)
  f <- withCallingHandlers(
    tryCatch(frm(cs$fo, family = fam_of(cs$fam), data = cs$d),
             error = function(e) e),
    warning = function(w) {
      cond <<- c(cond, paste0("W:", class(w)[1], ":", conditionMessage(w)))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      cond <<- c(cond, paste0("M:", class(m)[1], ":", conditionMessage(m)))
      invokeRestart("muffleMessage")
    })
  if (inherits(f, "error")) return(list(err = conditionMessage(f), cond = cond))
  vc <- tryCatch(VarCorr(f), error = function(e) NULL)
  ee <- unlist(lapply(vc, function(v) v$sd))
  list(par = f$opt$par, obj = f$opt$objective, conv = f$opt$convergence,
       ll = as.numeric(logLik(f)),
       V = tryCatch(vcov(f), error = function(e) conditionMessage(e)),
       vc = ee, autoscaled = isTRUE(f$autoscaled), cond = cond)
}
arm <- function(lib) {
  use(lib)
  lapply(cases, run1)
}
rB1 <- arm(BASE)
rL <- arm(LANE)
rB2 <- arm(BASE)
same <- function(a, b) identical(a, b)
rel <- function(a, b) {
  if (!is.numeric(a) || !is.numeric(b) || length(a) != length(b)) return(NA)
  max(abs(a - b) / pmax(abs(b), 1e-300))
}
options(width = 200)
out <- do.call(rbind, lapply(names(cases), function(nm) {
  b1 <- rB1[[nm]]; l <- rL[[nm]]; b2 <- rB2[[nm]]
  data.frame(case = nm, base_det = same(b1, b2),
             par = same(b1$par, l$par), ll = same(b1$ll, l$ll),
             V = same(b1$V, l$V), vc = same(b1$vc, l$vc),
             cond = same(b1$cond, l$cond),
             all = same(b1, l), rel_par = rel(l$par, b1$par),
             as_lane = isTRUE(l$autoscaled), ncond_b = length(b1$cond),
             ncond_l = length(l$cond))
}))
print(out)
cat("\nbase deterministic:", sum(out$base_det), "of", nrow(out),
    "; lane identical to base (all fields):", sum(out$all), "of", nrow(out),
    "\n")
for (nm in names(cases)) {
  if (!identical(rB1[[nm]]$cond, rL[[nm]]$cond)) {
    cat("\n== conditions differ:", nm, "\n base:\n")
    cat(paste0("  ", substr(rB1[[nm]]$cond, 1, 300)), sep = "\n")
    cat(" lane:\n")
    cat(paste0("  ", substr(rL[[nm]]$cond, 1, 300)), sep = "\n")
  }
}
