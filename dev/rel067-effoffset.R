# conditional_effects(effects = ) naming a variable that only an
# offset() reads, or brms's reserved `intercept`: brms's effect check
# reads brms:::get_all_effects(), which holds neither. Run per library:
#
#   Rscript dev/rel067-effoffset.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", lib, as.character(packageVersion("frmtmb")), "\n")
set.seed(20260930)
n <- 120
d <- data.frame(x = rnorm(n), time = runif(n, 1, 5))
d$yc <- rpois(n, d$time * exp(0.3 + 0.4 * d$x))
f1 <- yc ~ x + offset(log(time))
f2 <- yc ~ 0 + intercept + x
cat("brms get_all_effects, offset:",
    vapply(brms:::get_all_effects(brms::brmsterms(brms::bf(f1))), paste, ""),
    "\n")
cat("brms get_all_effects, 0 + intercept:",
    vapply(brms:::get_all_effects(brms::brmsterms(brms::bf(f2))), paste, ""),
    "\n")
one <- function(lab, fo, eff) {
  w <- character(0)
  r <- withCallingHandlers(tryCatch({
    fit <- frm(bf(fo), data = d, family = poisson())
    ce <- conditional_effects(fit, eff)
    paste("answers:", paste(names(ce), collapse = " "))
  }, error = function(e) paste("ERROR:", conditionMessage(e))),
  warning = function(cw) {
    w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning")
  })
  cat(sprintf("%-26s %s\n", lab, gsub("\n", " ", r)))
  for (x in w) cat("   warning:", gsub("\n", " ", x), "\n")
}
one("offset, 'time'", f1, "time")
one("offset, c('x', 'time')", f1, c("x", "time"))
one("intercept, 'intercept'", f2, "intercept")
