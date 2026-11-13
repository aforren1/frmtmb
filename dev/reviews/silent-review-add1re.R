# add1() of a new random-effect term: what refit does the row hold?
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r20",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-silent-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
set.seed(11)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), k = sample(1:4, n, TRUE),
                g = factor(sample(letters[1:8], n, TRUE)))
d$y <- 1 + 0.6 * d$k + 0.4 * d$z + sin(2 * d$x) + rnorm(8)[d$g] +
  rnorm(n, 0, 0.5)
f1 <- frm(y ~ x, data = d)
a <- tryCatch(add1(f1, ~ . + (1 | g), test = "Chisq"),
              error = function(e) paste("ERROR:", conditionMessage(e)))
print(a)
ll <- function(f) as.numeric(logLik(f))
h1 <- frm(y ~ x + (1 | g), data = d)
h2 <- frm(y ~ (x + 1 | g), data = d)
cat(sprintf("hand y ~ x + (1|g): LRT %.4f Df %d\n", 2 * (ll(h1) - ll(f1)),
            attr(logLik(h1), "df") - attr(logLik(f1), "df")))
cat(sprintf("hand y ~ (x + 1 | g): LRT %.4f Df %d\n", 2 * (ll(h2) - ll(f1)),
            attr(logLik(h2), "df") - attr(logLik(f1), "df")))
u <- update(f1, ~ . + 1 | g)
cat("update(f1, ~ . + 1 | g) formula:", deparse(u$bform$formula), "\n")
# string scope, as step() passes one
b <- tryCatch(add1(f1, "1 | g"), error = function(e)
  paste("ERROR:", conditionMessage(e)))
print(b)
