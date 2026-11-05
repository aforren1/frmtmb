# N2, N3, N4, N6, N7, N8 of fix round 2, against hand deletion and brms. Seeds
# 12000 and up. Usage: Rscript dev/refuse-f2-chain.R
.libPaths(c("C:/Users/adf44/source/r/wt-refuse-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
pr <- function(...) cat(sprintf(...), "\n", sep = "")
tr <- function(e) tryCatch(suppressMessages(suppressWarnings(e)),
                           error = function(err) structure(list(
                             msg = conditionMessage(err)), class = "tr_err"))
is_err <- function(v) inherits(v, "tr_err")
Gam <- stats::Gamma(link = "log")
ll <- function(f) as.numeric(logLik(f))
gap <- function(a, b) {
  if (is_err(a)) return(paste("ERR", substr(a$msg, 1, 90)))
  if (is_err(b)) return(paste("REF ERR", substr(b$msg, 1, 90)))
  ea <- unlist(a$estimates[c("beta", "betad")])
  eb <- unlist(b$estimates[c("beta", "betad")])
  sprintf("logLik gap %.2e; max coef gap %.2e; latent %d vs %d", ll(a) - ll(b),
          max(abs(ea - eb)), length(a$estimates$miss), length(b$estimates$miss))
}

# the review's chain: Gamma, z missing on 6 of the 10
set.seed(12000)
n <- 120
d <- data.frame(x = rnorm(n)); d$y <- rgamma(n, 2, 2 / exp(0.2 + 0.4 * d$x))
d$z <- 0.5 + 0.6 * d$y + 0.3 * d$x + rnorm(n, 0, 0.7)
M <- seq(4, n, by = 12); d$y[M] <- NA
dz <- d; dz$z[M[1:6]] <- NA
a <- tr(frm(bf(y | mi() ~ x) + Gam + bf(z | mi() ~ mi(y) + x) + gaussian(),
            data = dz))
b <- tr(frm(bf(y | mi() ~ x) + Gam + bf(z ~ mi(y) + x) + gaussian(),
            data = dz[-M[1:6], ]))
pr("chain two, Gamma: %s", gap(a, b))

# the review's exgaussian: 40 of 60 rows
set.seed(9503)
n <- 200
x <- rnorm(n); eta <- 0.3 + 0.5 * x
d <- data.frame(x = x, y = eta + rnorm(n, 0, 0.5) + rexp(n, 1 / 2))
d$z <- 0.5 + 0.6 * d$y + rnorm(n, 0, 0.7)
rows <- sort(sample(n, 60))
d$y[rows] <- NA
d$z[rows[1:40]] <- NA
a <- tr(frm(bf(y | mi() ~ x) + exgaussian() + bf(z | mi() ~ mi(y)) + gaussian(),
            data = d))
b <- tr(frm(bf(y | mi() ~ x) + exgaussian() + bf(z ~ mi(y)) + gaussian(),
            data = d[-rows[1:40], ]))
pr("chain two, exgaussian: %s", gap(a, b))
if (!is_err(a)) {
  pr("  y dpars %s (review: -0.003, 0.484 held; -0.404, 0.689 deleted)",
     paste(round(a$estimates$betad[grep("^y_", names(a$estimates$betad))], 3), collapse = " "))
}

# a chain of three: x | mi() ~ w, y | mi() ~ mi(x), z | mi() ~ mi(y)
set.seed(12001)
n <- 150
d <- data.frame(w = rnorm(n)); d$a <- rgamma(n, 2, 2 / exp(0.1 + 0.3 * d$w))
d$b <- rgamma(n, 3, 3 / exp(0.1 + 0.3 * d$a))
d$c <- 0.2 + 0.5 * d$b + rnorm(n, 0, 0.5)
gone3 <- 1:6     # missing at all three: every row leaves
gone2 <- 7:12    # a and b missing, c observed: both stay
gone1 <- 13:18   # a missing, b observed: a stays
d3 <- d
d3$a[c(gone3, gone2, gone1)] <- NA
d3$b[c(gone3, gone2)] <- NA
d3$c[gone3] <- NA
f3 <- bf(a | mi() ~ w) + Gam + bf(b | mi() ~ mi(a)) + Gam + bf(c | mi() ~ mi(b)) +
  gaussian()
a <- tr(frm(f3, data = d3))
b <- tr(frm(bf(a | mi() ~ w) + Gam + bf(b | mi() ~ mi(a)) + Gam +
              bf(c ~ mi(b)) + gaussian(), data = d3[-gone3, ]))
pr("chain three: %s", gap(a, b))
if (!is_err(a)) pr("  latent a %d (12), b %d (6), c %d (0)",
                   length(a$frame$mi_map$a$rows), length(a$frame$mi_map$b$rows),
                   length(a$frame$mi_map$c$rows))

# a fork: y read by z1 and z2; a row leaves only when both readers have gone
set.seed(12002)
n <- 150
d <- data.frame(x = rnorm(n)); d$y <- rgamma(n, 2, 2 / exp(0.2 + 0.4 * d$x))
d$z1 <- 0.5 + 0.6 * d$y + rnorm(n, 0, 0.7)
d$z2 <- -0.2 + 0.4 * d$y + rnorm(n, 0, 0.5)
df <- d
both <- 1:5; one <- 6:10
df$y[c(both, one)] <- NA
df$z1[c(both, one)] <- NA
df$z2[both] <- NA
ff <- bf(y | mi() ~ x) + Gam + bf(z1 | mi() ~ mi(y)) + gaussian() +
  bf(z2 | mi() ~ mi(y)) + gaussian()
a <- tr(frm(ff, data = df))
b <- tr(frm(bf(y | mi() ~ x) + Gam + bf(z1 | mi() ~ mi(y)) + gaussian() +
              bf(z2 ~ mi(y)) + gaussian(), data = df[-both, ]))
pr("fork: %s", gap(a, b))
if (!is_err(a)) pr("  latent y %d (5), z1 %d (0), z2 %d (0)",
                   length(a$frame$mi_map$y$rows), length(a$frame$mi_map$z1$rows),
                   length(a$frame$mi_map$z2$rows))

# a reader weight of 0 takes its row's read away
set.seed(12003)
n <- 120
d <- data.frame(x = rnorm(n)); d$y <- rgamma(n, 2, 2 / exp(0.2 + 0.4 * d$x))
d$z <- 0.5 + 0.6 * d$y + rnorm(n, 0, 0.7)
M <- seq(4, n, by = 12); d$y[M] <- NA; d$wz <- ifelse(is.na(d$y), 0, 1)
d$wz2 <- d$wz; d$wz2[M[1:4]] <- 1
a <- tr(frm(bf(y | mi() ~ x) + Gam + bf(z | weights(wz) ~ mi(y)) + gaussian(), data = d))
b <- tr(frm(bf(y ~ x) + Gam + bf(z ~ y) + gaussian(), data = d[-M, ]))
if (!is_err(a) && !is_err(b)) {
  pr("reader weight 0, all: latent %d; logLik gap %.2e", length(a$estimates$miss), ll(a) - ll(b))
}
a <- tr(frm(bf(y | mi() ~ x) + Gam + bf(z | weights(wz2) ~ mi(y)) + gaussian(), data = d))
b <- tr(frm(bf(y | mi() ~ x) + Gam + bf(z ~ mi(y)) + gaussian(), data = d[-M[5:10], ]))
pr("reader weight 0, 6 of 10: %s", gap(a, b))

# fitted at a row both left: the reader takes the value's mean, not the placeholder
set.seed(12004)
n <- 120
d <- data.frame(x = rnorm(n)); d$y <- rgamma(n, 2, 2 / exp(0.2 + 0.8 * d$x))
d$z <- 0.5 + 0.6 * d$y + rnorm(n, 0, 0.7)
M <- seq(4, n, by = 12); d$y[M] <- NA; d$z[M] <- NA
a <- tr(frm(bf(y | mi() ~ x) + Gam + bf(z | mi() ~ mi(y)) + gaussian(), data = d))
if (!is_err(a)) {
  fy <- fitted(a, resp = "y")[M, 1]; fz <- fitted(a, resp = "z")[M, 1]
  b <- a$estimates$beta
  pr("fitted(z) at a row both left, minus b0 + b1 fitted(y): %.2e; residual NA %d of 10",
     max(abs(fz - b[["z_(Intercept)"]] - b[["z_miy"]] * fy)),
     sum(is.na(residuals(a, resp = "z")[M, 1])))
}

# N6: a reader through idx that skips every missing row
set.seed(11600)
n <- 120
x <- rnorm(n); y0 <- rgamma(n, 2, 2 / exp(0.2 + 0.4 * x))
d <- data.frame(x = x, y = y0, id = seq_len(n))
d$z <- 0.5 + 0.6 * y0 + rnorm(n, 0, 0.7)
M <- seq(4, n, by = 12); d$y[M] <- NA
d$ok <- TRUE; d$ok[M] <- FALSE; d$z2 <- ifelse(d$ok, d$z, 0)
a <- tr(frm(bf(y | mi() + index(id) ~ x) + Gam +
              bf(z2 | subset(ok) ~ mi(y, idx = id) + x) + gaussian(), data = d))
for (me in c("fitted", "predict", "residuals")) for (r in c("z2", "y")) {
  v <- tr(do.call(me, list(a, resp = r)))
  pr("N6 %s(resp = %s): %s", me, r, if (is_err(v)) v$msg else "ok")
}

# N4, N7, N8
set.seed(11500)
n <- 120; lev <- paste0("L", 1:8)
dm <- data.frame(g1 = sample(lev, n, TRUE), g2 = sample(lev, n, TRUE), x = rnorm(n))
re <- rnorm(8); res <- 0.7 * re + rnorm(8, 0, 0.2); names(re) <- names(res) <- lev
dm$y <- 1 + (re[dm$g1] + re[dm$g2]) / 2 +
  rnorm(n, 0, exp(-0.2 + (res[dm$g1] + res[dm$g2]) / 2))
dm$w1 <- runif(n); dm$w2 <- runif(n); dm$W <- cbind(dm$w1, dm$w2)
dm$h <- sample(lev, n, TRUE)
for (txt in c(
  "bf(y ~ x + (1 | mm(g1, g2, id = 'p')), sigma ~ (1 | mm(g1, g2, id = 'p', cor = FALSE)))",
  "bf(y ~ x + (1 | mm(g1, g2, id = 'p', cor = FALSE)), sigma ~ (1 | mm(g1, g2, id = 'p')))",
  "bf(y ~ x + (1 | p | mm(g1, g2)), sigma ~ (1 | p | mm(g1, g2, cor = FALSE)))",
  "bf(y ~ x + (1 | p | mm(g1, g2, cor = FALSE)), sigma ~ (1 | p | mm(g1, g2)))",
  "bf(y ~ x + (1 | mm(g1, g2, id = 'p', cor = FALSE)), sigma ~ (1 | mm(g1, g2, id = 'p', cor = FALSE)))",
  "bf(y ~ x + (1 | mm(g1, g2, id = 'p', weights = W)), sigma ~ (1 | mm(g1, g2, id = 'p', weights = cbind(w1, w2))))",
  "bf(y ~ x + (1 | mm(g1, g2, id = 'p', weights = cbind(w1, w2))), sigma ~ (1 | mm(g1, g2, id = 'p', weights = W)))",
  "bf(y ~ x + (1 | mm(g1, g2, id = 'p', weights = cbind(w1, w2))), sigma ~ (1 | mm(g1, g2, id = 'p', weights = cbind(w2, w1))))",
  "bf(y ~ x + (1 | gr(g1, id = 'p')), sigma ~ (1 | gr(h, id = 'p', cor = FALSE)))",
  "bf(y ~ x + (1 | gr(g1, id = 'p', cor = FALSE)), sigma ~ (1 | gr(h, id = 'p')))",
  "bf(y ~ x + (1 | gr(g1, id = 'p', cor = FALSE)), sigma ~ (1 | gr(h, id = 'p', cor = FALSE)))",
  "bf(y ~ x + (1 | p | gr(g1, cor = FALSE)), sigma ~ (1 | p | gr(h, cor = FALSE)))",
  "bf(y ~ x + (1 | p | gr(g1)), sigma ~ (1 | p | gr(h, cor = FALSE)))"
)) {
  a <- tr(frm(eval(parse(text = txt), envir = asNamespace("frmtmb")), data = dm))
  co <- tryCatch(brms::make_stancode(eval(parse(text = txt), envir = asNamespace("brms")),
                                     data = dm), error = function(e) e)
  pr("%s\n   frm: %s\n   brms: %s", txt,
     if (is_err(a)) paste("refuses:", substr(a$msg, 1, 130)) else
       sprintf("fits, logLik %.4f", ll(a)),
     if (inherits(co, "error")) paste("refuses:", substr(conditionMessage(co), 1, 90)) else
       paste("builds, cor_1", grepl("cor_1", co)))
}
a <- tr(frm(bf(y ~ x + (1 | p | mm(g1, g2, weights = W)),
               sigma ~ (1 | p | mm(g1, g2, weights = W))), data = dm))
b <- tr(frm(bf(y ~ x + (1 | mm(g1, g2, id = "p", weights = W)),
               sigma ~ (1 | mm(g1, g2, id = "p", weights = cbind(w1, w2)))), data = dm))
if (!is_err(a) && !is_err(b)) pr("N7: W twice minus W beside cbind: %.2e", ll(a) - ll(b))

# N8: importance on all unread, vcov(cluster =), mi_sampling("bad")
set.seed(12005)
n <- 240; g <- factor(rep(1:24, each = 10)); u <- rnorm(24, 0, 0.6)
d <- data.frame(g = g, x = rnorm(n))
d$y <- rgamma(n, 2, 2 / exp(0.2 + 0.4 * d$x + u[g])); M <- seq(4, n, by = 12); d$y[M] <- NA
set.seed(1); a <- tr(frm(bf(y | mi() ~ x + (1 | g)), family = Gam, data = d, importance = 200))
set.seed(1); b <- tr(frm(bf(y ~ x + (1 | g)), family = Gam, data = d[-M, ], importance = 200))
pr("importance, all unread: %s", if (is_err(a)) a$msg else gap(a, b))
f <- tr(frm(bf(y | mi() ~ x), family = Gam, data = d))
fc <- tr(frm(bf(y ~ x), family = Gam, data = d[-M, ]))
va <- tr(vcov(f, cluster = d$g))
vb <- tr(vcov(fc, cluster = d$g[-M]))
pr("vcov(cluster =), all unread: %s", if (is_err(va)) va$msg else
  sprintf("max gap to complete case %.2e", max(abs(va - vb))))
m <- tr(frmtmb::mi_sampling(f, "bad"))
pr("mi_sampling(fit, 'bad'): %s", if (is_err(m)) m$msg else "taken")
