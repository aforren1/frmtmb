# Reviewer's short check of fix round 1: us_theta() and the changed
# us_sigma() paths. Lane library first.
.libPaths(c("C:/Users/adf44/source/r/wt-mvnlatent-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
show <- function(label, expr) {
  r <- tryCatch({
    v <- force(expr)
    paste("OK", paste(format(signif(unlist(v), 6)), collapse = " "))
  }, error = function(e) paste("ERROR:", conditionMessage(e)),
  warning = function(w) paste("WARNING:", conditionMessage(w)))
  cat(sprintf("%-40s %s\n", label, substr(gsub("\n", " ", r), 1, 170)))
}
cat("us_theta exported:", "us_theta" %in% getNamespaceExports("frmtmb"), "\n")

## refusals
show("non-square", us_theta(matrix(1, 2, 3)))
show("character", us_theta(matrix("1", 1, 1)))
show("vector", us_theta(c(1, 2)))
show("NA entry", us_theta(matrix(c(1, NA, NA, 1), 2)))
show("Inf entry", us_theta(matrix(c(Inf, 0, 0, 1), 2)))
show("asymmetric", us_theta(matrix(c(1, 0.5, 0.2, 1), 2)))
show("zero diagonal", us_theta(diag(c(1, 0))))
show("negative 1x1", us_theta(matrix(-1)))
show("pos diag, not PD", us_theta(matrix(c(1, 2, 2, 1), 2)))
show("PSD singular", us_theta(matrix(1, 2, 2)))
show("0 x 0", us_theta(matrix(numeric(0), 0, 0)))
show("zero matrix", us_theta(matrix(0, 2, 2)))
show("tiny scale symmetric", us_theta(diag(2) * 1e-300))
show("rounding asymmetry 1e-14", {
  S <- matrix(c(2, 0.3, 0.3 + 1e-14, 1), 2); us_theta(S)
})
show("integer matrix", us_theta(diag(2L)))
show("dimnames dropped", names(us_theta(matrix(c(1, .2, .2, 1), 2,
                              dimnames = list(c("a", "b"), c("a", "b"))))))

## round trips, d = 1..6, theta -> S -> theta and S -> theta -> S
set.seed(99)
for (d in 1:6) {
  np <- d + d * (d - 1) / 2
  th <- rnorm(np, 0, 0.7)
  S <- as.matrix(us_sigma(th, d))
  th2 <- us_theta(S)
  A <- crossprod(matrix(rnorm(d * (d + 3)), d + 3, d))
  A2 <- as.matrix(us_sigma(us_theta(A), d))
  cat(sprintf("round trip d=%d theta rel=%.2e  S rel=%.2e\n", d,
              max(abs(th2 - th)) / max(abs(th)),
              max(abs(A2 - A)) / max(abs(A))))
}
# a strongly correlated matrix
S <- matrix(0.999, 3, 3); diag(S) <- 1
cat(sprintf("round trip rho=0.999 S rel=%.2e\n",
            max(abs(as.matrix(us_sigma(us_theta(S), 3)) - S))))

## us_sigma changes
show("d=1 floored == sd2_floored", identical(
  as.matrix(us_sigma(-400, 1, floored = TRUE))[1, 1],
  frmtmb:::sd2_floored(-400)))
show("d = 2 + 1e-12 message", us_sigma(numeric(3), 2 + 1e-12))
show("d = 'a' message", us_sigma(numeric(3), "a"))
show("d = c(2,3) message", us_sigma(numeric(3), c(2, 3)))
# d = 4 order: theta entries fill L by column
th <- c(numeric(4), 1:6 / 10)
L <- diag(4)
L[2, 1] <- .1; L[3, 1] <- .2; L[4, 1] <- .3
L[3, 2] <- .4; L[4, 2] <- .5; L[4, 3] <- .6
L <- L / sqrt(rowSums(L^2))
show("d=4 column order", isTRUE(all.equal(as.matrix(us_sigma(th, 4)),
                                          L %*% t(L))))
