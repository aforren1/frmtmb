# The +3e18 logLik of the consumer sweep (seed 167): what the stop is, and
# whether the objective is a function of its argument (the inner start).
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
set.seed(167)
d <- expand.grid(t = 1:6, g = factor(1:25))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
d$z <- rnorm(n)
d$f <- factor(sample(letters[1:3], n, TRUE))
e <- as.vector(replicate(25, arima.sim(list(ar = 0.5), 6, sd = 0.5)))
re <- rep(rnorm(25, 0, .4), each = 6)
d$y <- rpois(n, exp(0.5 + 0.4 * d$x + e + re))
d$bin <- rbinom(n, 1, plogis(0.3 * d$x + e))
d$pos <- rgamma(n, 4, 4 / exp(0.5 + 0.4 * d$x + e))
saveRDS(d, "dev/review-latres-out/gamma-d167.rds")
m <- msg_of(frm(pos ~ x + ar(t, g), family = Gamma(link = "log"), data = d))
f <- m$val
print(m$err)
cat(m$warn, sep = "\n")
print(f$opt$convergence); print(f$opt$message)
print(logLik(f))
print(f$opt$par)
print(capture.output(print(summary(f))))
print(AIC(f))
print(head(fitted(f)))
fx <- frm(pos ~ x, family = Gamma(link = "log"), data = d)
cat("no-latent logLik", as.numeric(logLik(fx)), "\n")
ft <- msg_of(frm(pos ~ x + ar(t, g, cov = TRUE), family = Gamma(link = "log"), data = d))
cat("cov=TRUE logLik", as.numeric(logLik(ft$val)), " code", ft$val$opt$convergence, "\n")
print(ft$val$opt$par)
# Is the objective reproducible: same point, two inner starts
par <- f$opt$par
env <- f$obj$env
r <- env$random
# restart from zero
env$last.par.best[r] <- 0
v0 <- f$obj$fn(par)
env$last.par.best[r] <- rnorm(length(r), 0, .3)
v1 <- f$obj$fn(par)
cat("fn at stop from two inner starts:", v0, v1, "\n")
# a sane point near the data's truth
par2 <- par
par2[] <- c(0.5, 0.4, log(0.5), qlogis(0.75), log(4))[seq_along(par)]
print(names(par))
for (k in 1:5) {
  env$last.par.best[r] <- rnorm(length(r), 0, .3)
  cat("fn at truth, random inner start", k, ":", f$obj$fn(par2), "\n")
}
env$last.par.best[r] <- 0
cat("fn at truth from zero:", f$obj$fn(par2), "\n")
# the same, from frm()'s own starts at the truth
fs <- msg_of(frm(pos ~ x + ar(t, g), family = Gamma(link = "log"), data = d,
                 start = list(theta = c(log(0.5), qlogis(0.75)))))
print(fs$err)
if (is.null(fs$err)) {
  cat("started at truth: logLik", as.numeric(logLik(fs$val)), " code", fs$val$opt$convergence, "\n")
  print(fs$val$opt$par)
  cat(fs$warn, sep = "\n")
}
