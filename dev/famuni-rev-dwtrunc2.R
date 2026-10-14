# Reviewer debug: discrete_weibull() trunc(2, 9), per row.
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1014)
n <- 400
x <- rnorm(n)
q <- plogis(1.2 + 0.4 * x)
yt <- pmax(ceiling((log1p(-runif(n)) / log(q))^(1 / 1.1) - 1), 0)
d <- data.frame(x = x, y = yt)[yt >= 2 & yt <= 9, ]
f <- frm(y | trunc(lb = 2, ub = 9) ~ x, family = discrete_weibull(), data = d)
print(fixef(f))
print(summary(f)$spec_pars)
eta <- frm_linpred(f, dpar = "mu", type = "link")
mu <- plogis(eta)
k <- frm_linpred(f, dpar = "shape", type = "response")
fam <- family(f)
lp <- fam$lpdf(d$y, list(mu = mu, shape = k, .eta_mu = eta), list())
S <- function(v) mu^(v^k)
pm <- log(S(d$y) - S(d$y + 1))
cat("family lpdf vs hand pmf, max abs diff:", max(abs(lp - pm)), "\n")
lo <- fam$lcdf(rep(1, nrow(d)), list(mu = mu, shape = k, .eta_mu = eta), list())
hi <- fam$lcdf(rep(9, nrow(d)), list(mu = mu, shape = k, .eta_mu = eta), list())
cat("F(1) family vs hand:", max(abs(lo - (1 - S(2)))), " F(9):",
    max(abs(hi - (1 - S(10)))), "\n")
cat("sum lpdf", sum(lp), " sum log window", sum(log(hi - lo)),
    " sum lpdf - log window", sum(lp - log(hi - lo)), " logLik",
    as.numeric(logLik(f)), "\n")
cat("rows with log pmf > log window:", sum(lp > log(hi - lo)), "of", nrow(d), "\n")
i <- which.max(lp - log(hi - lo))
cat("worst row: y", d$y[i], "mu", mu[i], "k", k[i], "pmf", exp(lp[i]),
    "window", hi[i] - lo[i], "\n")
# the objective's own value at the optimum
cat("obj$fn at opt:", f$obj$fn(f$opt$par), "\n")
cat("n rows in fit:", nobs(f), " rows in d:", nrow(d), "\n")
