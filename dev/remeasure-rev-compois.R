# Reviewer probe for R6: what mu means in brms 2.23.0's com_poisson and
# in frmtmb 0.70.0's compois, by FITTING both on the same data.
# Seed 20261013; n = 400 draws from a COM-Poisson with lambda = 9,
# nu = 2 (so lambda^(1/nu) = 3), by inversion of the exact pmf.
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb)})
set.seed(20261013)
lam <- 9; nu <- 2
k <- 0:200
lp <- k * log(lam) - nu * lfactorial(k)
pk <- exp(lp - max(lp)); pk <- pk / sum(pk)
cat(sprintf("true: mean %.6f, lambda^(1/nu) %.6f, nu %.3f\n",
            sum(k * pk), lam^(1 / nu), nu))
y <- sample(k, 400, replace = TRUE, prob = pk)
d <- data.frame(y = y)
cat(sprintf("sample mean %.6f\n", mean(y)))
## brms's Stan code for the family
sc <- brms::stancode(y ~ 1, data = d, family = brms::brmsfamily("com_poisson"))
fl <- strsplit(sc, "\n")[[1]]
i <- grep("com_poisson_lpmf|com_poisson_log_lpmf", fl)
cat("## brms Stan lines naming com_poisson:\n")
cat(fl[unique(pmin(length(fl), c(i, i + 1)))], sep = "\n")
j <- grep("real com_poisson_log_lpmf", fl)
if (length(j)) cat(fl[j:(j + 14)], sep = "\n")
## frmtmb
ff <- frm(bf(y ~ 1), family = compois(), data = d)
print(coef(summary(ff)))
cat(sprintf("frmtmb exp(Intercept) %.6f\n", exp(fixef(ff)[1, 1])))
## brms
bf <- brms::brm(y ~ 1, data = d, family = brms::brmsfamily("com_poisson"),
                chains = 2, iter = 1000, seed = 20261013, refresh = 0,
                silent = 2)
s <- brms::fixef(bf)
print(s)
print(summary(bf)$spec_pars)
cat(sprintf("brms exp(Intercept) %.6f\n", exp(s[1, 1])))
