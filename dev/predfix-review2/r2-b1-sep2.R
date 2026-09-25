source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-helpers.R")
# Separation under an engaged default, part 2: (a) the estimates and
# standard errors each arm reports, seed 511, design of r2-b1-sep.R;
# (b) the ordinary-data trigger, a random slope on a column spread 0.03,
# seeds 531..540; (c) is the warning the PRE-FIT's? Count the warnings
# the pre-fit raises by tracing autoscale_prefit's fit_assembled call.
run <- function(...) {
  w <- character(0)
  f <- withCallingHandlers(tryCatch(frm(...), error = function(e) e),
         warning = function(x) { w <<- c(w, conditionMessage(x))
                                 invokeRestart("muffleWarning") })
  list(f = f, n = length(w), w = w)
}
set.seed(511)
n <- 240
d <- data.frame(x = rnorm(n) * 1e-4, z = rnorm(n))
d$yb <- as.integer(d$z > 0)
a <- run(yb ~ z + x, family = bernoulli(), data = d)
b <- run(yb ~ z + x, family = bernoulli(), data = d,
         control = frmtmb_control(autoscale = FALSE))
cat("(a) default: code", a$f$opt$convergence, "warnings", a$n, "\n")
print(signif(fixef(a$f)[, 1:2], 4))
cat("    FALSE: code", b$f$opt$convergence, "warnings", b$n, "\n")
print(signif(fixef(b$f)[, 1:2], 4))
dg <- tryCatch(capture.output(diagnose(a$f)), error = function(e) conditionMessage(e))
cat("    diagnose(default) lines mentioning separation/boundary/large:",
    sum(grepl("separat|bound|large|infinite|non-finite", dg, ignore.case = TRUE)),
    "of", length(dg), "\n")
cat(head(dg, 25), sep = "\n")
gl <- withCallingHandlers(glm(yb ~ z + x, binomial, d),
        warning = function(w) { cat("    glm warns:", conditionMessage(w), "\n")
                                invokeRestart("muffleWarning") })

cat("\n(b) random slope on a column spread 0.03, separated z\n")
tab <- NULL
for (seed in 531:540) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 12))
  d2 <- data.frame(g, x = rnorm(240) * 0.03, z = rnorm(240))
  d2$yb <- as.integer(d2$z > 0)
  a2 <- run(yb ~ z + x + (1 + x | g), family = bernoulli(), data = d2)
  b2 <- run(yb ~ z + x + (1 + x | g), family = bernoulli(), data = d2,
            control = frmtmb_control(autoscale = FALSE))
  ok <- function(r) !inherits(r$f, "error")
  if (!ok(a2) || !ok(b2)) cat("seed", seed, "ERROR default:",
    if (!ok(a2)) conditionMessage(a2$f) else "none", "| FALSE:",
    if (!ok(b2)) conditionMessage(b2$f) else "none", "
")
  tab <- rbind(tab, data.frame(seed,
    engaged = ok(a2) && !is.null(a2$f$par_units),
    warnD = a2$n, warnF = b2$n,
    codeD = if (ok(a2)) a2$f$opt$convergence else NA,
    codeF = if (ok(b2)) b2$f$opt$convergence else NA,
    bz_D = if (ok(a2)) fixef(a2$f)["z", 1] else NA,
    bz_F = if (ok(b2)) fixef(b2$f)["z", 1] else NA))
}
print(tab)
cat("warned under FALSE, silent under default:",
    sum(tab$warnF > 0 & tab$warnD == 0), "of", nrow(tab), "\n")
