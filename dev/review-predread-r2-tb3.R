# Short check of fix round 2: where fitted(newdata = ) on 500 rows of a me() fit stops, and
# whether the same row count stops on a plain fit, on the other readers, and on rellib-r16.
# Seed 7974. REV_LIB picks the build.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7974)
n <- 300
o <- data.frame(lat = rnorm(n), sx = 0.4, z = rnorm(n))
o$xo <- o$lat + rnorm(n, 0, 0.4)
o$cnt <- rpois(n, exp(0.3 + 0.5 * o$lat))
fit <- quiet(frm(bf(cnt ~ me(xo, sx) + z), family = poisson(), data = o))
pl <- quiet(frm(bf(cnt ~ xo + z), family = poisson(), data = o))
big <- data.frame(xo = rnorm(500L), sx = 0.4, z = 0)
sh <- function(lab, e) {
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(paste(dim(as.matrix(suppressWarnings(suppressMessages(e)))), collapse = "x"),
                error = function(err) paste("ERROR:", substr(conditionMessage(err), 1, 70)))
  cat(sprintf("%-30s %6.2f s  %s\n", lab, proc.time()[["elapsed"]] - t0, r))
}
sh("plain fitted", fitted(pl, newdata = big))
sh("me fitted", fitted(fit, newdata = big))
sh("me frm_linpred", frm_linpred(fit, newdata = big))
sh("me frm_linpred se", frm_linpred(fit, newdata = big, se.fit = TRUE)$se.fit)
sh("me predict", predict(fit, newdata = big, ndraws = 5))
sh("me simulate", simulate(fit, nsim = 2, newdata = big))
sh("me frm_curve", frm_curve(fit, newdata = big, simultaneous = FALSE)$.se)
sh("me get_predict", marginaleffects::get_predict(fit, newdata = big))
options(error = NULL)
e <- tryCatch(fitted(fit, newdata = big), error = function(e) e)
print(conditionCall(e))
tr <- NULL
withCallingHandlers(
  try(suppressWarnings(fitted(fit, newdata = big)), silent = TRUE),
  error = function(e) tr <<- sys.calls())
cat(paste(tail(vapply(tr, function(cl) substr(paste(deparse(cl), collapse = " "), 1, 90), ""), 14), collapse = "\n"), "\n")
