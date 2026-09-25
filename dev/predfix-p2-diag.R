# Punch round 2, requirement (c): on the reviewer's two designs the
# default must warn at least wherever autoscale = FALSE warns, and never
# error where FALSE returns a fit. Also records whether it engaged or
# fell back, and TRUE's behavior beside it.
#   PREDFIX_ARM=lane Rscript dev/predfix-p2-diag.R > dev/predfix-log/p2-diag.txt
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
run <- function(fo, fam, d, a) {
  k <- 0L
  f <- withCallingHandlers(
    tryCatch(frm(fo, family = fam, data = d,
                 control = frmtmb_control(autoscale = a)),
             error = function(e) e),
    warning = function(w) {
      k <<- k + 1L
      invokeRestart("muffleWarning")
    })
  list(err = inherits(f, "error"), warn = k,
       engaged = !inherits(f, "error") && !is.null(f$par_units),
       ll = if (inherits(f, "error")) NA else as.numeric(logLik(f)),
       msg = if (inherits(f, "error")) conditionMessage(f) else "")
}
rows <- list()
add <- function(design, seed, fo, fam, d) {
  a <- run(fo, fam, d, NULL)
  b <- run(fo, fam, d, FALSE)
  t <- run(fo, fam, d, TRUE)
  rows[[length(rows) + 1L]] <<- data.frame(design, seed,
    errD = a$err, errF = b$err, errT = t$err, warnD = a$warn,
    warnF = b$warn, warnT = t$warn, engagedD = a$engaged,
    llD_minus_llF = a$ll - b$ll,
    T_error_names_FALSE = t$err && grepl("autoscale = FALSE", t$msg,
                                         fixed = TRUE))
}
for (seed in 511:520) {
  set.seed(seed)
  n <- 240
  d <- data.frame(x = rnorm(n) * 1e-4, z = rnorm(n))
  d$yb <- as.integer(d$z > 0)
  add("A separated", seed, yb ~ z + x, bernoulli(), d)
}
for (seed in 531:540) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 12))
  d2 <- data.frame(g, x = rnorm(240) * 0.03, z = rnorm(240))
  d2$yb <- as.integer(d2$z > 0)
  add("B slope 0.03", seed, yb ~ z + x + (1 + x | g), bernoulli(), d2)
}
r <- do.call(rbind, rows)
options(width = 200)
print(r)
cat("\ndefault errors where FALSE returns:", sum(r$errD & !r$errF),
    "; default warns less than FALSE:", sum(!r$errD & r$warnD < r$warnF),
    "; default engaged:", sum(r$engagedD),
    "; TRUE errors:", sum(r$errT), "of which name autoscale = FALSE:",
    sum(r$T_error_names_FALSE), "\n")
