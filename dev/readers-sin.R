# The sin(a)^2 + c case of dev/review-readers-sc.R part D (data seed 7208):
# fitted() against the exact mean, and the rule sizes the net reads.
source("dev/readers-setup.R")
ns <- asNamespace("frmtmb")
set.seed(7208)
n <- 800
d <- data.frame(l1 = rnorm(n), l2 = rnorm(n), s1 = 0.2, s2 = 0.2)
d$x1 <- d$l1 + rnorm(n, 0, 0.2)
d$x2 <- d$l2 + rnorm(n, 0, 0.2)
d$y <- sin(0.3 + 0.8 * d$l1)^2 + (0.5 + 0.5 * d$l2) + rnorm(n, 0, 0.1)
fit <- frm(bf(y ~ sin(a)^2 + c, a ~ me(x1, s1), c ~ me(x2, s2), nl = TRUE),
           family = gaussian(), data = d)
lp <- function(nd, dp) as.numeric(frm_linpred(fit, newdata = nd, dpar = dp))
q <- data.frame(x1 = 0, x2 = 0, s1 = 1, s2 = 1)
ba <- lp(transform(q, x1 = 1), "a") - lp(q, "a")
nd <- data.frame(x1 = -lp(q, "a") / ba, x2 = 0, s1 = c(pi, 10 * pi) / abs(ba),
                 s2 = 0.5)
sizes <- integer()
trace("me_gh_moments", quote(sizes <<- c(sizes, K)), where = ns,
      print = FALSE)
f <- fitted(fit, newdata = nd)
untrace("me_gh_moments", where = ns)
exact <- (1 - exp(-2 * (c(pi, 10 * pi))^2)) / 2 + lp(nd, "c")
cat("a at the row", lp(nd, "a"), " fitted", f[, "Estimate"], " exact", exact,
    " relative", f[, "Estimate"] / exact - 1, "\n")
cat("rule sizes read:", sizes, "\n")
