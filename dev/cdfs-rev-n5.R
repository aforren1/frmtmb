# Reviewer, fix round 2: N5, censored counts at or above the trials
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(2)
d <- data.frame(x = rnorm(80), y = rbinom(80, 12, 0.5), tr = 12, cc = 0)
d$y2 <- d$y
try1 <- function(lab, dd, f = y | trials(tr) + cens(cc, y2) ~ x, fam = binomial()) {
  r <- tryCatch({ fit <- frm(f, family = fam, data = dd); sprintf("fits, logLik %.6f", as.numeric(logLik(fit))) },
                error = function(e) paste("REFUSED:", conditionMessage(e)))
  cat(sprintf("%-46s %s\n", lab, r))
}
a <- d; a$y[1:3] <- 12; a$cc[1:3] <- 1; try1("binomial right-censored at trials", a)
b <- d; b$y[1:3] <- 11; b$cc[1:3] <- 1; try1("binomial right-censored at trials - 1", b)
e <- d; e$y[1:3] <- 0; e$cc[1:3] <- -1; try1("binomial left-censored at 0", e)
g <- d; g$y[1:3] <- 12; g$y2[1:3] <- 13; g$cc[1:3] <- 2; try1("binomial interval starting at trials", g)
try1("beta_binomial right-censored at trials", a, fam = beta_binomial())
