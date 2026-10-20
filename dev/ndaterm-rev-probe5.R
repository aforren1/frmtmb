# Reviewer of lane ndaterm: subset() responses on newdata. Seed 7.
# Usage: Rscript dev/ndaterm-rev-probe5.R <arm>   (mine | r9)
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "mine"
base <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "mine") file.path(base, "wt-ndaterm-lib"),
            file.path(base, "rellib-r9"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, find.package("frmtmb"), "\n")
set.seed(7)
n <- 200
d <- data.frame(x = rnorm(n), N = sample(5:15, n, TRUE),
                sb = rep(c(TRUE, FALSE), n / 2))
d$y <- rbinom(n, d$N, plogis(0.3 + 0.5 * d$x))
d$gy <- 1 + d$x + rnorm(n)
f <- frm(mvbf(bf(y | trials(N) + subset(sb) ~ x), bf(gy ~ x)),
         family = list(binomial(), gaussian()), data = d)
nd <- data.frame(x = c(0, 0), N = c(1000, 2), sb = TRUE)
set.seed(1)
p <- tryCatch(predict(f, newdata = nd, resp = "y", ndraws = 400),
              error = function(e) paste("ERROR:", conditionMessage(e)))
print(p)
print(fitted(f, newdata = nd, resp = "y"))
p2 <- tryCatch(predict(f, newdata = nd[, c("x", "N")], resp = "y",
                       ndraws = 20),
               error = function(e) paste("ERROR:", conditionMessage(e)))
cat("newdata without sb:\n"); print(p2)
p3 <- tryCatch(predict(f, newdata = data.frame(x = 0), resp = "gy",
                       ndraws = 20),
               error = function(e) paste("ERROR:", conditionMessage(e)))
cat("resp = gy, newdata x only:\n"); print(p3)
cat("done\n")
