# Punch round 2, optional: refit() with an NA in newresp. Seed 7.
# FORMROBUST_LIB="" for the base arm.
LIB <- Sys.getenv("FORMROBUST_LIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
set.seed(7)
d <- data.frame(x = rnorm(50)); d$y <- 1 + d$x + rnorm(50)
f <- frm(bf(y ~ x), data = d)
r <- tryCatch(refit(f, replace(d$y, 1, NA)), error = function(e) e,
              warning = function(w) w)
if (inherits(r, "condition")) cat(class(r)[1], ":", conditionMessage(r), "\n") else
  cat("logLik", format(logLik(r)), "\n")
