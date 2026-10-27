## brms 2.23.0's mo() coding, measured with standata() and stancode()
## (no compile). Run:
##   Rscript dev/mocode-brms-coding.R > dev/mocode-out/brms-coding.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
show <- function(lab, x, f = y ~ mo(x), nd = NULL, extra = NULL) {
  set.seed(1)
  d <- data.frame(x = x, y = rnorm(length(x)))
  if (!is.null(extra)) d <- cbind(d, extra)
  r <- tryCatch(q(standata(f, data = d)), error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-28s ERROR: %s\n", lab, conditionMessage(r)))
    return(invisible())
  }
  xm <- r$Xmo_1
  cat(sprintf("%-28s Jmo = %s; Xmo_1 = %s; attr min of x: %s\n", lab,
              paste(r$Jmo, collapse = ","),
              paste(xm, collapse = " "),
              format(attr(d$x, "min"))))
}
show("0..3", c(0, 1, 2, 3, 0, 1))
show("1..4", c(1, 2, 3, 4, 1, 2))
show("2..5", c(2, 3, 4, 5, 2, 3))
show("-2..1", c(-2, -1, 0, 1, -2, -1))
show("gaps 1,2,5", c(1, 2, 5, 1, 2, 5))
show("gaps 0,3", c(0, 3, 0, 3))
show("two values 1,2", c(1, 2, 1, 2))
show("two values 0,1", c(0, 1, 0, 1))
show("single value 3", c(3, 3, 3))
show("double 1.0..3.0", c(1, 2, 3, 1))
show("non-integer 1.5", c(1, 1.5, 2))
show("integer type 1L..4L", c(1L, 2L, 3L, 4L))
x0 <- c(1, 2, 3, 4); attr(x0, "min") <- 0
show("1..4 with attr min 0", x0)
show("ordered a<b<c<d", factor(c("a", "b", "c", "d", "a"), ordered = TRUE))
show("ordered unused level e",
     factor(c("a", "b", "c", "a"), levels = c("a", "b", "c", "e"),
            ordered = TRUE))
show("ordered, unused first", factor(c("b", "c", "d"),
     levels = c("a", "b", "c", "d"), ordered = TRUE))
show("unordered factor", factor(c("a", "b", "c")))
show("mo(x)*z, 1..4", c(1, 2, 3, 4, 1, 2), f = y ~ mo(x) * z,
     extra = data.frame(z = c(0.5, -1, 2, 0.1, 0.3, 1)))
show("mo(x + 1), 1..4", c(1, 2, 3, 4, 1, 2), f = y ~ mo(x + 1))
show("(mo(x) | g), 1..4", c(1, 2, 3, 4, 1, 2, 3, 4),
     f = y ~ mo(x) + (mo(x) | g),
     extra = data.frame(g = rep(1:2, 4)))

cat("\n== stancode for 1..4 (data and parameters blocks, mo lines)\n")
d <- data.frame(x = c(1, 2, 3, 4, 1, 2), y = rnorm(6))
sc <- q(stancode(y ~ mo(x), data = d))
l <- strsplit(sc, "\n")[[1]]
cat(grep("mo|simo|Jmo", l, value = TRUE), sep = "\n")

cat("\n== default_prior for 1..4\n")
print(q(default_prior(y ~ mo(x), data = d))[, c("prior", "class", "coef")])

cat("\n== newdata: standata(fit, newdata) on an empty (chains = 0) fit\n")
set.seed(2)
d <- data.frame(x = sample(1:4, 40, TRUE)); d$y <- d$x + rnorm(40)
fit <- q(brm(y ~ mo(x), data = d, chains = 0, empty = TRUE))
nd_show <- function(lab, x) {
  r <- tryCatch(q(standata(fit, newdata = data.frame(x = x, y = 0))),
                error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-26s ERROR: %s\n", lab, conditionMessage(r)))
  } else {
    cat(sprintf("%-26s Jmo = %s; Xmo_1 = %s\n", lab, r$Jmo,
                paste(r$Xmo_1, collapse = " ")))
  }
}
nd_show("newdata 2..4", 2:4)
nd_show("newdata 1..4", 1:4)
nd_show("newdata 3", 3)
nd_show("newdata 0 (below min)", c(0, 1))
nd_show("newdata 5 (above max)", c(4, 5))
nd_show("newdata 2.5 (non-int)", c(2.5))
cat("brms's stored data minimum: min(fit$data$x) =", min(fit$data$x), "\n")
cat("fit$data$x has attr min:", format(attr(fit$data$x, "min")), "\n")

cat("\n== where brms drops unused ordered levels\n")
d2 <- data.frame(x = factor(c("b", "c", "d", "b"), levels = c("a", "b", "c", "d"),
                            ordered = TRUE), y = rnorm(4))
fit2 <- q(brm(y ~ mo(x), data = d2, chains = 0, empty = TRUE))
cat("levels of fit$data$x:", levels(fit2$data$x), "\n")
r <- tryCatch(q(standata(fit2, newdata = data.frame(
  x = factor("a", levels = c("a", "b", "c", "d"), ordered = TRUE), y = 0))),
  error = function(e) e)
cat("newdata level a (unused at fit):",
    if (inherits(r, "error")) conditionMessage(r) else
      paste("Jmo", r$Jmo, "Xmo_1", r$Xmo_1), "\n")

cat("\n== attr(x, 'min') set by the user at the fit, then newdata\n")
d3 <- data.frame(x = c(1, 2, 3, 4, 1, 2), y = rnorm(6))
attr(d3$x, "min") <- 0
fit3 <- q(brm(y ~ mo(x), data = d3, chains = 0, empty = TRUE))
s3 <- q(standata(fit3))
cat("fit: Jmo", s3$Jmo, "Xmo_1", s3$Xmo_1, "\n")
s3n <- q(standata(fit3, newdata = data.frame(x = c(1, 2, 3, 4), y = 0)))
cat("newdata 1..4: Jmo", s3n$Jmo, "Xmo_1", s3n$Xmo_1, "\n")
