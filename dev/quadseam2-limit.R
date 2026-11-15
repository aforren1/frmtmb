# Item 4867: a coefficient that runs off under one-block quadrature.
# Poisson, a factor level whose every count is zero, plus (1 | g)
# (dev/reviews/quadseam-review-limit.R), Laplace and quadrature; plus the
# absent case (the same data without the zero level) on both.
#   Rscript quadseam2-limit.R <libpath|base>
a <- commandArgs(TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
run <- function(lbl, d, fo, fam, quad) {
  ms <- character(); w <- character()
  t0 <- proc.time()[["elapsed"]]
  f <- withCallingHandlers(
    frm(fo, family = fam, data = d, quadrature = quad),
    warning = function(x) {
      w <<- c(w, paste(class(x)[1], substr(conditionMessage(x), 1, 160)))
      invokeRestart("muffleWarning")
    },
    message = function(x) {
      ms <<- c(ms, paste(class(x)[1], substr(conditionMessage(x), 1, 160)))
      invokeRestart("muffleMessage")
    })
  v <- withCallingHandlers(
    tryCatch(sqrt(diag(vcov(f))), error = function(e) NA),
    warning = function(x) {
      w <<- c(w, paste(class(x)[1], substr(conditionMessage(x), 1, 160)))
      invokeRestart("muffleWarning")
    },
    message = function(x) {
      ms <<- c(ms, paste(class(x)[1], substr(conditionMessage(x), 1, 160)))
      invokeRestart("muffleMessage")
    })
  t1 <- proc.time()[["elapsed"]] - t0
  b <- f$estimates$beta
  cat(sprintf("%s QUAD %s ll %.8f conv %d secs %.1f\n", lbl, quad,
              as.numeric(logLik(f)), f$opt$convergence, t1))
  cat("   est", paste(sprintf("%s=%.4g", names(b), b), collapse = " "), "\n")
  cat("   se ", paste(sprintf("%s=%s", names(v), format(v, digits = 4)),
                      collapse = " "), "\n")
  sl <- f$cache$se_lost %||% attr(f$cache$sdr, "se_lost")
  for (x in c(w, ms)) cat("   ", gsub("\\s+", " ", x), "\n")
  invisible(f)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
set.seed(5)
d <- data.frame(g = factor(rep(1:20, each = 9)),
                f = factor(rep(c("a", "b", "c"), 60)), x = rnorm(180))
d$y <- rpois(180, exp(0.5 + 0.3 * d$x + rnorm(20, 0, 0.5)[d$g]))
d0 <- d
d$y[d$f == "c"] <- 0L
for (quad in c(FALSE, TRUE)) {
  run("runoff", d, y ~ f + x + (1 | g), poisson(), quad)
}
# absent case: the same design with the level's counts kept
for (quad in c(FALSE, TRUE)) {
  run("absent", d0, y ~ f + x + (1 | g), poisson(), quad)
}
# bernoulli: a level whose every response is 0 (separation on a factor)
set.seed(7)
db <- data.frame(g = factor(rep(1:30, each = 6)),
                 f = factor(rep(c("a", "b", "c"), 60)), x = rnorm(180))
db$y <- rbinom(180, 1, plogis(0.2 + 0.5 * db$x + rnorm(30, 0, 0.8)[db$g]))
db0 <- db
db$y[db$f == "c"] <- 0L
for (quad in c(FALSE, TRUE)) {
  run("bernsep", db, y ~ f + x + (1 | g), bernoulli(), quad)
}
for (quad in c(FALSE, TRUE)) {
  run("bernabs", db0, y ~ f + x + (1 | g), bernoulli(), quad)
}
cat("END\n")
