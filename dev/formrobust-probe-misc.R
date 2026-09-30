# Probes of paths the NEWS names: pp_check(newdata = ) on a -1/-2
# bernoulli fit, and core predict(re_formula = NA) with an NA response
# under cov = FALSE. Seeds as in the tests.
.libPaths(c("C:/Users/adf44/source/r/wt-formrobust-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(5)
n <- 120
d <- data.frame(x = rnorm(n))
yy <- rbinom(n, 1, plogis(-0.3 + 1.2 * d$x))
d$y01 <- yy
d$ym <- ifelse(yy == 1, -1, -2)
fm <- frm(bf(ym ~ x), data = d, family = bernoulli())
f0 <- frm(bf(y01 ~ x), data = d, family = bernoulli())
nd <- d[d$ym == -2, ][1:6, ]
pdf(NULL)
a <- tryCatch(pp_check(fm, newdata = nd, ndraws = 5, type = "bars"),
              error = function(e) conditionMessage(e))
b <- tryCatch(pp_check(f0, newdata = nd, ndraws = 5, type = "bars"),
              error = function(e) conditionMessage(e))
cat("pp_check class:", class(a)[1], class(b)[1], "\n")
if (is.list(a) && !is.null(a$data)) {
  cat("pp_check data identical:", identical(a$data, b$data), "\n")
}
yv <- frmtmb:::pp_check_newdata_y(fm, fm$spec$responses[[1]], nd)
cat("pp_check_newdata_y on -2 rows:", yv, "\n")
set.seed(31)
G <- 30; Tn <- 8
dd <- expand.grid(t = 1:Tn, g = factor(1:G))
dd$x <- rnorm(nrow(dd))
dd$u <- rep(rnorm(G, 0, 0.5), each = Tn)
dd$y <- 1 + 0.5 * dd$x + dd$u + rnorm(nrow(dd))
fa <- frm(bf(y ~ x + (1 | g) + ar(t, g)), data = dd)
nd2 <- dd[dd$g %in% c("1", "2"), ]
nd2$y[nd2$t >= 6] <- NA
set.seed(2)
p <- predict(fa, newdata = nd2, re_formula = NA, ndraws = 50)
cat("predict re_formula = NA, NA response: anyNA", anyNA(p), "\n")
