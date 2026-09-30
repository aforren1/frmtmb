# Item 1 probe, part 2: each expression term against the same model
# with the expression precomputed as a column, at fit time and on
# newdata; multivariate; index(); mi(sd); vreal(). Seed 11.
LIB <- Sys.getenv("FORMROBUST_LIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
set.seed(11)
n <- 60
d <- data.frame(x = rnorm(n), wt = runif(n, 0.5, 2), time = runif(n, 1, 3),
                n = rpois(n, 5) + 2L, s = runif(n, 0.2, 0.6),
                g = gl(6, 10), id = sample(1000, n))
d$yc <- rpois(n, exp(0.2 + 0.3 * d$x) * 2 * d$time)
d$yb <- rbinom(n, d$n + 1, plogis(0.3 * d$x))
d$y <- rnorm(n, 0.5 * d$x, 1)
d$y2 <- rnorm(n, -0.5 * d$x, 1)
d$c <- sample(c(0, 1), n, TRUE, prob = c(0.8, 0.2))
d$lb <- -3
# precomputed columns
d$w2 <- d$wt * 2; d$t2 <- d$time * 2; d$n1 <- d$n + 1; d$s2 <- d$s / 2
d$lbm <- d$lb - 1; d$id2 <- d$id * 3
d$sx <- d$x > 0; d$xm <- d$x; d$xm[c(3, 9)] <- NA
nd <- d[c(1:2, 4:8), ]
cmp <- function(lab, fa, fb, resp = NULL) {
  r <- tryCatch({
    A <- fa(); B <- fb()
    same_ll <- identical(as.numeric(logLik(A)), as.numeric(logLik(B)))
    pa <- fitted(A, newdata = nd, resp = resp)
    pb <- fitted(B, newdata = nd, resp = resp)
    sim <- function(F) if (is.null(resp)) {
      as.matrix(as.data.frame(simulate(F, nsim = 3, seed = 2, newdata = nd)))
    } else 0
    sa <- sim(A); sb <- sim(B)
    set.seed(3); qa <- predict(A, newdata = nd, resp = resp)
    set.seed(3); qb <- predict(B, newdata = nd, resp = resp)
    paste("ll identical:", same_ll, "| fitted identical:",
          identical(unname(as.matrix(pa)), unname(as.matrix(pb))),
          "| simulate identical:", identical(unname(sa), unname(sb)),
          "| predict identical:",
          identical(unname(as.matrix(qa)), unname(as.matrix(qb))))
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-22s %s\n", lab, r))
}
cmp("weights(wt * 2)",
    function() frm(bf(y | weights(wt * 2) ~ x), data = d),
    function() frm(bf(y | weights(w2) ~ x), data = d))
cmp("rate(time * 2)",
    function() frm(bf(yc | rate(time * 2) ~ x), data = d, family = poisson()),
    function() frm(bf(yc | rate(t2) ~ x), data = d, family = poisson()))
cmp("trials(n + 1)",
    function() frm(bf(yb | trials(n + 1) ~ x), data = d,
                   family = binomial()),
    function() frm(bf(yb | trials(n1) ~ x), data = d, family = binomial()))
cmp("se(s / 2)",
    function() frm(bf(y | se(s / 2) ~ x), data = d),
    function() frm(bf(y | se(s2) ~ x), data = d))
cmp("se(s / 2, sigma=TRUE)",
    function() frm(bf(y | se(s / 2, sigma = TRUE) ~ x), data = d),
    function() frm(bf(y | se(s2, sigma = TRUE) ~ x), data = d))
cmp("trunc(lb = lb - 1)",
    function() frm(bf(y | trunc(lb = lb - 1) ~ x), data = d),
    function() frm(bf(y | trunc(lb = lbm) ~ x), data = d))
cmp("mv weights(wt * 2)",
    function() frm(bf(y | weights(wt * 2) ~ x) + bf(y2 | se(s / 2) ~ x),
                   data = d),
    function() frm(bf(y | weights(w2) ~ x) + bf(y2 | se(s2) ~ x), data = d),
    resp = "y")
cmp("mv subset+rate",
    function() frm(bf(yc | rate(time * 2) + subset(x > 0) ~ x,
                      family = poisson()) + bf(y ~ x), data = d),
    function() {
      frm(bf(yc | rate(t2) + subset(sx) ~ x, family = poisson()) +
            bf(y ~ x), data = d)
    }, resp = "yc")
cmp("vreal(wt * 2) gaussian?",
    function() frm(bf(y | weights(wt * 2) + mi() ~ x), data = d),
    function() frm(bf(y | weights(w2) + mi() ~ x), data = d))
k <- 2
cmp("weights(wt * k) env k",
    function() frm(bf(y | weights(wt * k) ~ x), data = d),
    function() frm(bf(y | weights(w2) ~ x), data = d))
cmp("index(id * 3) + mi(idx)",
    function() {
      frm(bf(y | subset(x > -5) ~ mi(xm, idx = id2)) +
            bf(xm | mi() + index(id * 3) ~ 1), data = d)
    },
    function() {
      frm(bf(y | subset(x > -5) ~ mi(xm, idx = id2)) +
            bf(xm | mi() + index(id2) ~ 1), data = d)
    }, resp = "y")
