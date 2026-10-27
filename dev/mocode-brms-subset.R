## brms 2.23.0, mo(k) in a subset() response of a multivariate model:
## how new data are coded (standata on an empty fit; nothing compiles).
## Copied from the reviewer's dev/review-mocode-brms-subset.R and rerun.
##   Rscript dev/mocode-brms-subset.R
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] + rnorm(n)
d$y2 <- d$y + rnorm(n)
d$s1 <- d$k > 1L
d$s2 <- TRUE
f <- bf(y | subset(s1) ~ mo(k)) + bf(y2 | subset(s2) ~ mo(k)) +
  set_rescor(FALSE)
sd0 <- q(standata(f, data = d))
cat("fit: Jmo_y =", sd0$Jmo_y, " Xmo_y_1 head =", head(sd0$Xmo_y_1),
    " k on those rows =", head(d$k[d$s1]), "\n")
fit <- q(brm(f, data = d, empty = TRUE))
for (x in list(2:4, 1:4)) {
  nd <- data.frame(k = x, s1 = TRUE, s2 = TRUE, y = 0, y2 = 0)
  r <- tryCatch(q(standata(fit, newdata = nd)), error = function(e) e)
  cat("newdata k =", paste(x, collapse = " "), ": ")
  if (inherits(r, "error")) cat("ERROR:", conditionMessage(r), "\n") else
    cat("Jmo_y =", r$Jmo_y, " Xmo_y_1 =", r$Xmo_y_1, "; Jmo_y2 =", r$Jmo_y2,
        " Xmo_y2_1 =", r$Xmo_y2_1, "\n")
}

cat("\n== the same without subset(), and a univariate control\n")
f2 <- bf(y ~ mo(k)) + bf(y2 ~ mo(k)) + set_rescor(FALSE)
fit2 <- q(brm(f2, data = d, empty = TRUE))
r <- q(standata(fit2, newdata = data.frame(k = 2:4, y = 0, y2 = 0)))
cat("mv, no subset, newdata k = 2 3 4: Xmo_y_1 =", r$Xmo_y_1,
    "; Xmo_y2_1 =", r$Xmo_y2_1, "\n")
fit1 <- q(brm(y2 ~ mo(k), data = d, empty = TRUE))
r <- q(standata(fit1, newdata = data.frame(k = 2:4, y2 = 0)))
cat("univariate, newdata k = 2 3 4: Xmo_1 =", r$Xmo_1, "\n")
fit3 <- q(brm(bf(y | subset(s1) ~ mo(k)), data = d, empty = TRUE))
r <- tryCatch(q(standata(fit3, newdata = data.frame(k = 2:4, s1 = TRUE,
                                                     y = 0))),
              error = function(e) e)
cat("univariate subset, newdata k = 2 3 4:",
    if (inherits(r, "error")) conditionMessage(r) else r$Xmo_1, "\n")
