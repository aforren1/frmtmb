# Punch round 1: reproduce fitted() on disc ~ 0 + gp(x, k = 6) and see
# which lp slot is empty.
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(11)
n <- 120
d <- data.frame(x = round(runif(n, 0, 5), 1), z = rnorm(n))
lat <- sin(d$x) + 0.5 * d$z + rlogis(n)
d$y <- factor(cut(lat, quantile(lat, c(0, 0.3, 0.6, 1)),
                  include.lowest = TRUE, labels = FALSE), ordered = TRUE)
fd <- suppressWarnings(frm(bf(y ~ z, disc ~ 0 + gp(x, k = 6)), data = d,
                           family = cumulative()))
for (nm in names(fd$frame$linpreds)) {
  lp <- fd$frame$linpreds[[nm]]
  cat(nm, "| par", format(lp$par), "| idx", length(lp$idx), "| ncol X",
      ncol(lp$X), "| class est", class(fd$estimates[[lp$par]])[1], "\n")
}
r <- tryCatch(fitted(fd), error = function(e) e)
cat("fitted():", if (inherits(r, "error")) conditionMessage(r) else "ok",
    "\n")
nd <- data.frame(x = c(5.6, 2.55), z = c(0, 1))
r2 <- tryCatch(fitted(fd, newdata = nd), error = function(e) e)
cat("fitted(newdata):", if (inherits(r2, "error")) conditionMessage(r2)
    else "ok", "\n")
if (!inherits(r2, "error")) print(r2)
