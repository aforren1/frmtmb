# Short check of fix round 2, item 3: where fitted() on a cumulative me() fit stops with
# "variable names are limited to 10000 bytes", the row count at which it starts, by family,
# and the cost of the band where it answers. Seed 7973. REV_LIB picks the build.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7973)
n <- 300
o <- data.frame(lat = rnorm(n), sx = 0.4)
o$xo <- o$lat + rnorm(n, 0, 0.4)
o$ord <- as.integer(cut(0.4 + 0.7 * o$lat + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf), labels = FALSE))
o$cnt <- rpois(n, exp(0.3 + 0.5 * o$lat))
fits <- list(
  cumulative = quiet(frm(bf(ord ~ me(xo, sx)), family = cumulative(), data = o)),
  categorical = quiet(frm(bf(ord ~ me(xo, sx)), family = categorical(), data = o)),
  poisson = quiet(frm(bf(cnt ~ me(xo, sx)), family = poisson(), data = o)))
for (nm in names(fits)) {
  for (m in c(20L, 100L, 200L, 400L, 500L, 2000L, 20000L)) {
    big <- data.frame(xo = rnorm(m), sx = 0.4)
    t0 <- proc.time()[["elapsed"]]
    r <- tryCatch({
      a <- suppressWarnings(fitted(fits[[nm]], newdata = big))
      paste(dim(a), collapse = "x")
    }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 70)))
    cat(sprintf("%-12s %6d rows %6.2f s  %s\n", nm, m, proc.time()[["elapsed"]] - t0, r))
  }
}
cat("traceback of the stop:\n")
big <- data.frame(xo = rnorm(500L), sx = 0.4)
withCallingHandlers(
  tryCatch(fitted(fits$cumulative, newdata = big), error = function(e) NULL),
  error = function(e) {
    cs <- sys.calls()
    cat(paste(tail(vapply(cs, function(cl) substr(paste(deparse(cl), collapse = " "), 1, 100), ""), 12), collapse = "\n"), "\n")
  })
fn <- frmtmb:::response_mean_gh
print(names(formals(fn)))
