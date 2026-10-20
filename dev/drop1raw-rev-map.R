# Reviewer, item 5: what logLik() of a MAP fit is, in frmtmb and in
# glmmTMB (priors =), and what their drop1()/anova() do with it.
# Arg: lane. Log: dev/drop1raw-rev-out/map-lane.txt (seed 5, n = 200)
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(glmmTMB)})
cat("frmtmb from", dirname(find.package("frmtmb")), "| glmmTMB",
    format(packageVersion("glmmTMB")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(5)
n <- 200
d <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1, x3 = rnorm(n),
                g = factor(sample(letters[1:8], n, TRUE)))
d$z <- 1 + d$x + log(d$x2) + 0.1 * d$x3 + rnorm(n)

cat("\n== frmtmb: logLik(MAP fit) against the data log density ==\n")
fp <- q(frm(bf(z ~ x + log(x2) + x3), data = d,
            prior = set_prior("normal(0, 1)", class = "b")))
b <- fixef(fp)[, 1]
X <- model.matrix(~ x + log(x2) + x3, d)
mu <- drop(X %*% b)
# sigma at the MAP estimate: the ML sigma given b would differ, so read
# the fit's own value
sig <- unname(sigma(fp))[1]
lik <- sum(dnorm(d$z, mu, sig, log = TRUE))
lp <- sum(dnorm(b[-1], 0, 1, log = TRUE))
cat(sprintf("logLik(fit) %.10f | data log density %.10f | + log prior %.10f\n",
            as.numeric(logLik(fp)), lik, lik + lp))
cat(sprintf("logLik - (data + prior) = %.3e ; logLik - data = %.6f\n",
            as.numeric(logLik(fp)) - (lik + lp), as.numeric(logLik(fp)) - lik))
cat("names(estimates):", paste(names(fp$estimates), collapse = ","),
    "| sigma", sig, "\n")
dt <- q(drop1(fp, test = "Chisq"))
print(dt)
# the LRT for x3 by hand: refit without x3 under the same prior; then
# the data-only statistic at those MAP estimates, and the ML statistic
f3 <- q(frm(bf(z ~ x + log(x2)), data = d,
            prior = set_prior("normal(0, 1)", class = "b")))
b3 <- fixef(f3)[, 1]
mu3 <- drop(model.matrix(~ x + log(x2), d) %*% b3)
sig3 <- unname(sigma(f3))[1]
lik3 <- sum(dnorm(d$z, mu3, sig3, log = TRUE))
lp3 <- sum(dnorm(b3[-1], 0, 1, log = TRUE))
cat(sprintf("drop1 LRT(x3) %.6f = 2 (lp0 - lp3) %.6f + 2 (lik0 - lik3) %.6f\n",
            dt["x3", "LRT"], 2 * (lp - lp3), 2 * (lik - lik3)))
m0 <- lm(z ~ x + log(x2) + x3, d); m3 <- lm(z ~ x + log(x2), d)
cat(sprintf("ML LRT(x3) %.6f\n", 2 * (logLik(m0) - logLik(m3))))
# a bound-only prior adds no density: is it a MAP fit?
fbnd <- tryCatch(q(frm(bf(z ~ x + log(x2) + x3), data = d,
                       prior = set_prior("", class = "b", coef = "x",
                                         lb = 0))),
                 error = function(e) e)
if (inherits(fbnd, "error")) {
  cat("bound-only prior: ERROR", conditionMessage(fbnd), "\n")
} else {
  f0 <- q(frm(bf(z ~ x + log(x2) + x3), data = d))
  cat(sprintf("bound-only prior logLik %.10f | no prior %.10f\n",
              as.numeric(logLik(fbnd)), as.numeric(logLik(f0))))
  cat("summary marks MAP:", any(grepl("MAP", capture.output(print(fbnd)))),
      "| prior class:", class(fbnd$prior)[1], "\n")
}

cat("\n== glmmTMB priors = ==\n")
pr <- data.frame(prior = "normal(0,1)", class = "fixef", coef = "",
                 stringsAsFactors = FALSE)
gp <- glmmTMB(z ~ x + log(x2) + x3, data = d, priors = pr)
g0 <- glmmTMB(z ~ x + log(x2) + x3, data = d)
bg <- fixef(gp)$cond
siggp <- sigma(gp)
likg <- sum(dnorm(d$z, drop(X %*% bg), siggp, log = TRUE))
cat(sprintf("glmmTMB logLik(prior fit) %.10f | data log density %.10f | no-prior fit %.10f\n",
            as.numeric(logLik(gp)), likg, as.numeric(logLik(g0))))
cat(sprintf("glmmTMB: logLik - data = %.6f ; log prior of fixef = %.6f (with intercept %.6f)\n",
            as.numeric(logLik(gp)) - likg, sum(dnorm(bg[-1], 0, 1, log = TRUE)),
            sum(dnorm(bg, 0, 1, log = TRUE))))
r <- tryCatch(drop1(gp, test = "Chisq"), error = function(e) e,
              warning = function(w) w)
cat("glmmTMB drop1 on the prior fit:\n"); print(r)
gp3 <- glmmTMB(z ~ x + log(x2), data = d, priors = pr)
r2 <- tryCatch(anova(gp, gp3), error = function(e) e, warning = function(w) w)
cat("glmmTMB anova on prior fits:\n"); print(r2)
