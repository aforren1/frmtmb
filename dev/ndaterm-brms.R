# Lane ndaterm: what brms 2.23.0 does with addition terms on newdata.
# binomial trials(N) and gaussian se(s), the same data as
# dev/ndaterm-probe.R (seed 7). Log dev/ndaterm-out/brms.txt.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
try_it <- function(label, expr) {
  r <- tryCatch(expr, error = function(e) {
    structure(paste("ERROR:", conditionMessage(e)), class = "err")
  })
  if (inherits(r, "err")) cat(label, ":", r, "\n")
  invisible(r)
}
set.seed(7)
n <- 200
db <- data.frame(x = rnorm(n), N = sample(5:15, n, TRUE))
db$y <- rbinom(n, db$N, plogis(0.3 + 0.5 * db$x))
fb <- brm(y | trials(N) ~ x, family = binomial(), data = db, chains = 1,
          iter = 600, warmup = 300, refresh = 0, seed = 7)
ndb <- data.frame(x = c(0, 0), N = c(1000, 2))
pp <- posterior_predict(fb, newdata = ndb)
cat("brms posterior_predict(newdata N = 1000, 2) col means",
    format(colMeans(pp), digits = 5), " max", apply(pp, 2, max), "\n")
try_it("brms predict(newdata without N)",
       predict(fb, newdata = data.frame(x = 0)))
dg <- data.frame(x = rnorm(n), s = runif(n, 0.2, 1))
dg$y <- 1 + dg$x + rnorm(n, 0, dg$s)
fg <- brm(y | se(s) ~ x, family = gaussian(), data = dg, chains = 1,
          iter = 600, warmup = 300, refresh = 0, seed = 7)
ndg <- data.frame(x = c(0, 0), s = c(0.05, 5))
pg <- posterior_predict(fg, newdata = ndg)
cat("brms se(s) posterior_predict(newdata s = 0.05, 5) col sd",
    format(apply(pg, 2, sd), digits = 5), "\n")
try_it("brms predict(newdata without s)",
       predict(fg, newdata = data.frame(x = 0)))
cat("done\n")
