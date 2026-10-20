# Reviewer: step() against stats::step() on lm() and glm(), and three
# drop1() cases of the reviewer's own against lm/glm/lmer.
# Arg: lane. Log: dev/drop1raw-rev-out/step-lane.txt (seed 21, n = 300)
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(21)
n <- 300
d <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1, x3 = rnorm(n),
                x4 = rnorm(n), f = factor(sample(c("a", "b", "c"), n, TRUE)),
                g = factor(sample(letters[1:10], n, TRUE)))
d$z <- 1 + 0.5 * d$x + log(d$x2) + 0.12 * d$x3 + rnorm(10)[d$g] + rnorm(n)
d$cnt <- rpois(n, exp(0.2 + 0.3 * d$x + 0.3 * log(d$x2) + 0.08 * d$x4))

path <- function(s) {
  a <- s$anova
  paste(sprintf("%s:%.6f", trimws(a$Step), a$AIC), collapse = " ")
}
cat("\n== step(): lm ==\n")
f4 <- z ~ x + log(x2) + splines::ns(x3, 2) + x4
sl <- step(lm(f4, d), trace = 0)
sf <- q(step(frm(bf(f4), data = d), trace = 0))
cat("lm    ", deparse1(formula(sl)), "|", path(sl), "\n")
cat("frmtmb", deparse1(formula(sf)), "|", path(sf), "\n")
# frmtmb's AIC is -2 logLik + 2 (p + 1); lm's extractAIC drops the
# constant n (log(2 pi) + 1) + 2
k <- n * (log(2 * pi) + 1) + 2
cat("AIC diffs minus the constant:",
    paste(sprintf("%.2e", sf$anova$AIC - sl$anova$AIC - k), collapse = " "),
    "\n")
cat("\n== step(): glm poisson ==\n")
g4 <- cnt ~ x + log(x2) + poly(x3, 2) + x4
sg <- step(glm(g4, poisson, d), trace = 0)
sfg <- q(step(frm(bf(g4), family = poisson(), data = d), trace = 0))
cat("glm   ", deparse1(formula(sg)), "|", path(sg), "\n")
cat("frmtmb", deparse1(formula(sfg)), "|", path(sfg), "\n")
cat("AIC diffs:", paste(sprintf("%.2e", sfg$anova$AIC - sg$anova$AIC),
                        collapse = " "), "\n")
cat("\n== step() with NA rows in x4 ==\n")
dn <- d; dn$x4[c(2, 50, 77)] <- NA
r <- tryCatch(step(lm(f4, dn), trace = 0), error = function(e) e)
cat("lm:", if (inherits(r, "error")) conditionMessage(r) else
  deparse1(formula(r)), "\n")
r2 <- tryCatch(q(step(frm(bf(f4), data = dn), trace = 0)),
               error = function(e) e)
cat("frmtmb:", if (inherits(r2, "error")) conditionMessage(r2) else
  paste(deparse1(formula(r2)), "|", path(r2), "| nobs", nobs(r2)), "\n")

cat("\n== reviewer drop1() cases ==\n")
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1))
dn3 <- d; dn3$x[c(4, 80, 150)] <- NA
# (a) spline x factor interaction, NA rows in another variable, lm
fa <- z ~ splines::ns(x2, 3) * f + x
ma <- lm(fa, dn3); ta <- drop1(ma)
fra <- q(drop1(frm(bf(fa), data = dn3), test = "Chisq"))
la <- nobs(ma) * log(ta$RSS[-1] / ta$RSS[1])
cat(sprintf("%-40s rows %s | df %s | LRT rel diff %.2e\n",
            "lm ns(x2, 3) * f + x, NA in x", identical(rownames(fra),
                                                     rownames(ta)),
            identical(as.numeric(fra$Df[-1]), as.numeric(ta$Df[-1])),
            rel(fra$LRT[-1], la)))
# (b) poly() x factor interaction, poisson, with NA rows
fb <- cnt ~ poly(x2, 2) * f + scale(x4)
mb <- glm(fb, poisson, dn3); tb <- drop1(mb, test = "Chisq")
frb <- q(drop1(frm(bf(fb), family = poisson(), data = dn3), test = "Chisq"))
cat(sprintf("%-40s rows %s | df %s | LRT rel diff %.2e\n",
            "glm poly(x2, 2) * f + scale(x4)",
            identical(rownames(frb), rownames(tb)),
            identical(as.numeric(frb$Df[-1]), as.numeric(tb$Df[-1])),
            rel(frb$LRT[-1], tb$LRT[-1])))
# (c) lmer with scale() x factor and ns(), no NA (lme4 refuses NA rows)
fc <- z ~ scale(x) * f + splines::ns(x2, 3) + (1 | g)
mc <- lme4::lmer(fc, d, REML = FALSE); tc <- drop1(mc, test = "Chisq")
frc <- q(drop1(frm(bf(fc), data = d), test = "Chisq"))
cat(sprintf("%-40s rows %s | df %s | LRT rel diff %.2e\n",
            "lmer scale(x) * f + ns(x2, 3) + (1|g)",
            identical(rownames(frc), rownames(tc)),
            identical(as.numeric(frc$Df[-1]), as.numeric(tc$npar[-1])),
            rel(frc$LRT[-1], tc$LRT[-1])))
print(frc); print(tc)
