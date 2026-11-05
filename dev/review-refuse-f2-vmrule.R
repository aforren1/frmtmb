# Fix round 2, N5, von_mises branch of the ruling, on the data of the vm arm and
# the von_mises cells of table 1: the share of read rows whose mean direction is
# within 2 / sqrt(kappa) of the window's end. Also: update(), a fit saved without
# the field. Same seeds as the rule script (11600, 11207 for table 1's von_mises).
args <- "none"
txt <- readLines("dev/review-refuse-f2-rule.R")
i <- which(startsWith(txt, "mode <- args[1]"))
eval(parse(text = sub('args <- commandArgs(trailingOnly = TRUE)', 'args <- "none"', txt[seq_len(i - 1L)], fixed = TRUE)), envir = globalenv())
share <- function(a) {
  rows <- a$frame$mi_map$y$rows
  mu <- fitted(a, resp = "y", dpar = "mu")[rows, 1]
  k <- fitted(a, resp = "y", dpar = "kappa")[rows, 1]
  mean(pi - abs(mu) < 2 / sqrt(k))
}
set.seed(11600)
rvm <- function(mu, k) vapply(mu, function(m) { repeat { w <- runif(1, -pi, pi); if (runif(1) < exp(k * (cos(w) - 1))) break }; ((w + m + pi) %% (2 * pi)) - pi }, 0)
for (cfg in list(c(0, 3), c(3.05, 2), c(3.05, 3), c(3.14, 20))) for (mm in c(15, 60)) for (bz in c(0.2, 0.7)) for (k in 1:3) {
  n <- 200; x <- rnorm(n)
  y <- rvm(cfg[1] + 0.05 * x, cfg[2])
  d <- data.frame(x = x, y = y); d$z <- 0.5 + bz * y + 0.3 * x + rnorm(n, 0, 0.7)
  rows <- sample(n, mm); d$y[rows] <- NA
  a <- suppressWarnings(frm(bf(y | mi() ~ x) + von_mises() + bf(z ~ mi(y) + x) + gaussian(), data = d))
  pr("SHARE vm mean %.2f kappa %g read %d slope %.1f rep %d: share %.2f", cfg[1], cfg[2], mm, bz, k, share(a))
}
fam <- "von_mises"; sp <- speco(fam); set.seed(11200 + match(fam, fams))
for (cfg in list(c(200, 60, 0.7, 0.7), c(200, 60, 0.05, 1), c(200, 100, 0.7, 0.7), c(200, 15, 0.7, 0.7))) {
  n <- cfg[1]; mm <- cfg[2]; bz <- cfg[3]; sz <- cfg[4]
  x <- rnorm(n); y <- sp$rgen(0.3 + 0.5 * x)
  d <- data.frame(x = x, y = y); d$z <- 0.5 + bz * y + 0.3 * x + rnorm(n, 0, sz)
  d$y[sample(n, mm)] <- NA
  a <- suppressWarnings(frm(bf(y | mi() ~ x) + von_mises() + bf(z ~ mi(y) + x) + gaussian(), data = d))
  pr("SHARE t1 von_mises read %d slope %.2f: share %.2f", mm, bz, share(a))
}
# the field on update() and on a fit saved before it
set.seed(11901); n <- 150; x <- rnorm(n)
d <- data.frame(x = x, y = rgamma(n, 1.3, 1.3 / exp(0.3 + 0.5 * x)), w = rnorm(n)); d$z <- 0.5 + 0.7 * d$y + rnorm(n, 0, 0.7); d$y[seq(3, n, 3)] <- NA
a <- tr(frm(bf(y | mi() ~ x) + stats::Gamma(link = "log") + bf(z ~ mi(y) + x) + gaussian(), data = d))
pr("fit: warned %s; field %s", any(grepl("Laplace approximation of", attr(a, "w"))), paste(names(a$mi_laplace), collapse = ","))
u <- tr(update(a, . ~ . + w))
pr("update: %s; warned %s; field total %s against the first %s", msg(u, 100), any(grepl("Laplace approximation of", attr(u, "w"))), format(u$mi_laplace$total, digits = 4), format(a$mi_laplace$total, digits = 4))
s <- a; s$mi_laplace <- NULL
for (me in c("summary", "fitted", "logLik", "vcov")) pr("  field absent, %s: %s", me, msg(tr(do.call(me, list(s))), 80))
un <- tr(update(s, newdata = d)); pr("  field absent, update(newdata): %s", msg(un, 80))
print(system.time(frm(bf(y | mi() ~ x) + stats::Gamma(link = "log") + bf(z ~ mi(y) + x) + gaussian(), data = d, control = frmtmb_control(check_se = "ignore"))))
