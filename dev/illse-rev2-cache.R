# Reviewer of lane illse, round 2, item 7: the Hessian caches on
# obj$env. refit() to a new response, update() to new data, and a refit
# to the same response (whose optimum is the cached point), each against
# a fresh frm() on the same data; whether any objective environment is
# shared; and a direct probe: the cached AD Hessian after the data
# behind an objective change.
#   REV_ARM=lane Rscript dev/illse-rev2-cache.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
mk <- function(seed, deg = 8) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(200, 1, 2))
  d$y <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  d
}
fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8
cmpf <- function(a, b) {
  Va <- suppressWarnings(vcov(a, full = TRUE)); Vb <- suppressWarnings(vcov(b, full = TRUE))
  c(par = max(abs(a$opt$par - b$opt$par) / sqrt(diag(Vb))),
    se = max(abs(sqrt(diag(Va)) / sqrt(diag(Vb)) - 1)))
}
keys <- function(f) c(he = !is.null(f$obj$env$.frmtmb_he), rh = !is.null(f$obj$env$.frmtmb_reparam_h),
                      oh = !is.null(f$obj$env$.frmtmb_outer_h))
d <- mk(1)
f <- suppressWarnings(frm(fo, data = d))
cat("caches on the first fit:", keys(f), "\n")
d2 <- d; set.seed(99); d2$y <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
r1 <- suppressWarnings(refit(f, d2$y))
g1 <- suppressWarnings(frm(fo, data = d2))
cat("refit(new response) vs fresh frm: ", sprintf("%s %.1e", names(cmpf(r1, g1)), cmpf(r1, g1)),
    "| env shared with original", identical(r1$obj$env, f$obj$env), "\n")
r2 <- suppressWarnings(refit(f, d$y))
cat("refit(same response) vs original: ", sprintf("%s %.1e", names(cmpf(r2, f)), cmpf(r2, f)),
    "| env shared", identical(r2$obj$env, f$obj$env), "| cached point = original optimum",
    identical(r2$obj$env$.frmtmb_he$p, unname(f$opt$par)), "\n")
u1 <- suppressWarnings(update(f, data = d2))
cat("update(data = d2) vs fresh frm:   ", sprintf("%s %.1e", names(cmpf(u1, g1)), cmpf(u1, g1)),
    "| env shared", identical(u1$obj$env, f$obj$env), "\n")
u2 <- suppressWarnings(update(f, . ~ . - x8))
g2 <- suppressWarnings(frm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7, data = d))
cat("update(. ~ . - x8) vs fresh frm:  ", sprintf("%s %.1e", names(cmpf(u2, g2)), cmpf(u2, g2)), "\n")
# with (1 | g): the outer-Hessian cache
dg <- d; dg$g <- factor(rep(1:20, length.out = 200)); set.seed(5)
dg$y <- dg$y + stats::rnorm(20, 0, 0.5)[dg$g]
fg <- suppressWarnings(frm(update(fo, . ~ . + (1 | g)), data = dg))
dg2 <- dg; set.seed(6); dg2$y <- dg$y + stats::rnorm(200, 0, 0.05)
rg <- suppressWarnings(refit(fg, dg2$y)); gg <- suppressWarnings(frm(update(fo, . ~ . + (1 | g)), data = dg2))
cat("(1|g) refit vs fresh frm:         ", sprintf("%s %.1e", names(cmpf(rg, gg)), cmpf(rg, gg)),
    "| env shared", identical(rg$obj$env, fg$obj$env), "\n")
# probe: does any frmtmb path change an objective's data in place? Mutate
# the data behind a fitted objective and ask the cache at the same point
p <- f$opt$par
H0 <- ns$obj_he_cached(f$obj, p)
ok <- tryCatch({
  f$obj$env$data$y <- f$obj$env$data$y + 1   # not something frmtmb does
  H1 <- ns$obj_he_cached(f$obj, p)
  identical(H0, H1)
}, error = function(e) paste("error:", conditionMessage(e)))
cat("after an in-place data change the cache serves the old Hessian:", ok,
    "(the key is the point only)\n")
