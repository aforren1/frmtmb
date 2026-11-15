# Re-check of fix round 2, follow-ups: bitwise compare of part 4b, the
# data = globalenv() fallback, and the theta of R5.
source("C:/Users/adf44/source/r/frmtmb-wt-refits/dev/reviews/refits-review-lib.R")
out <- "C:/Users/adf44/source/r/frmtmb-wt-refits/dev/reviews"
a <- readRDS(file.path(out, "refits-short-r2-infl-before.rds"))
b <- readRDS(file.path(out, "refits-short-r2-infl-after.rds"))
for (k in names(a)) {
  cat(sprintf("4b %-10s identical(rellib-r21, lane): %s\n", k,
              identical(a[[k]], b[[k]])))
}
G <- 20
set.seed(1)
dd <- data.frame(g = factor(rep(seq_len(G), each = 6)),
                 cl = factor(rep(seq_len(G / 2), each = 12)),
                 x = rnorm(G * 6))
dd$y <- 1 + 0.5 * dd$x + rnorm(G * 6) + rnorm(G / 2)[dd$cl]
se_of <- function(v) round(sqrt(diag(v)), 7)
yv <- dd$y; xv <- dd$x; clv <- dd$cl
fenv <- frm(yv ~ xv, data = globalenv())
go("N2c' data = globalenv(): ~ clv vs explicit",
   list(se_of(vcov_cluster(fenv, ~ clv)), se_of(vcov_cluster(fenv, clv))))
# R5 again: the theta of each deletion, refit vs fresh fit
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(8)
nn <- 30
du <- data.frame(x = rnorm(nn), w = runif(nn, 0.5, 2), N = 10L,
                 g = factor(rep(1:6, each = 5)), t = rep(1:5, 6))
du$y <- rnorm(nn)
du$k <- rbinom(nn, 10, 0.4)
du$s <- rbinom(nn, 1, 0.5)
du <- du[sample(nn), ]
ynew <- rnorm(nn, 2)
f <- q(frm(y ~ x + (1 | g), data = du))
r <- q(refit(f, ynew))
d2 <- du; d2$y <- ynew
f2 <- q(frm(y ~ x + (1 | g), data = d2))
ia <- q(influence(r, groups = "g"))
ib <- q(influence(f2, groups = "g"))
go("R5 theta (log sd_g) per deletion: refit | fresh",
   cbind(refit = ia$theta[, 1], fresh = ib$theta[, 1]))
go("R5 sd_g = exp(theta): max |diff|",
   max(abs(exp(ia$theta[, 1]) - exp(ib$theta[, 1]))))
go("R5 full fits: theta refit vs fresh", c(r$estimates$theta, f2$estimates$theta))
