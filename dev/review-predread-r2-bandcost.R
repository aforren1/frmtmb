# Short check of fix round 2, item 3: the draw's own quantile for a category, on expected
# counts of a multinomial and on a row of NA, and what it costs on many new rows.
# Seed 7972. REV_LIB picks the build.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7972)
n <- 300
d <- data.frame(lat = rnorm(n), sx = 0.5, tr = 20L)
d$xo <- d$lat + rnorm(n, 0, 0.5)
pr <- cbind(1, exp(0.3 + 0.6 * d$lat), exp(-0.2 + 0.9 * d$lat))
pr <- pr / rowSums(pr)
d$Y <- t(vapply(seq_len(n), function(i) as.integer(rmultinom(1L, 20L, pr[i, ])), integer(3L)))
colnames(d$Y) <- c("a", "b", "c")
fit <- quiet(frm(bf(Y | trials(tr) ~ me(xo, sx)), family = multinomial(), data = d))
nd <- data.frame(xo = c(0, 1, NA), sx = c(0.5, 1.2, 0.5), tr = 20L)
f <- suppressWarnings(fitted(fit, newdata = nd))
print(round(f, 3))
cf <- fixef(fit)[, "Estimate"]
cat(names(cf), "\n")
b <- matrix(cf, 2L, 2L, byrow = FALSE)
if (grepl("mu", names(cf)[1L])) {
  S <- 40000L
  for (i in 1:2) {
    l <- rnorm(S, nd$xo[i], nd$sx[i])
    e <- cbind(1, exp(cf[grep("mub_Intercept|mub_\\(Intercept\\)", names(cf))] + cf[grep("mub_me", names(cf))] * l),
               exp(cf[grep("muc_Intercept|muc_\\(Intercept\\)", names(cf))] + cf[grep("muc_me", names(cf))] * l))
    p <- 20 * e / rowSums(e)
    cat("row", i, "fixed coefficient draw q:", round(apply(p, 2L, quantile, c(0.025, 0.975)), 3), "\n")
  }
}
# cost on many new rows
o <- data.frame(lat = rnorm(n), sx = 0.4)
o$xo <- o$lat + rnorm(n, 0, 0.4)
o$ord <- as.integer(cut(0.4 + 0.7 * o$lat + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf), labels = FALSE))
fo <- quiet(frm(bf(ord ~ me(xo, sx)), family = cumulative(), data = o))
for (m in c(500L, 5000L, 20000L)) {
  big <- data.frame(xo = rnorm(m), sx = 0.4)
  t1 <- replicate(2L, system.time(fitted(fo, newdata = big))[["elapsed"]])
  cat(sprintf("fitted(cumulative me(), %d new rows): %.2f %.2f s\n", m, t1[1L], t1[2L]))
}
