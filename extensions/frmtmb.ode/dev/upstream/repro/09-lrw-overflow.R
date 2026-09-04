## The ceiling is an integer overflow inside deSolve, not an RTMBode defect.
##
## deSolve sizes lsoda's real work array as
##   lrw = 22 + neq * max(16, neq + 9)
## in R *integer* arithmetic.  That exceeds .Machine$integer.max at
##   neq = 46337  (22 + 46337 * 46346 = 2,147,534,624 > 2^31 - 1)
## after which the length reaching the allocator is NA, and the call dies with
## a nonsense message ("cannot allocate memory block of size 134217728 Tb").
## Below the overflow the call succeeds but still allocates 8*neq^2 bytes.
##
## RTMBode reaches those sizes on its own: the augmented system at order k has
##   neq = sum(nstate * (nstate + nparms)^(0:k))
## and a Laplace outer gradient needs order 3.
##
## Usage: Rscript 09-lrw-overflow.R

f <- function(t, y, p) list(-y)
probe <- function(neq, method) {
    y <- rep(1, neq)
    tryCatch({ deSolve::ode(y, c(0, 1), f, NULL, method = method); "ok" },
             error = function(e) paste("ERROR:", conditionMessage(e)))
}
lrw <- function(neq) 22 + neq * max(16, neq + 9)

cat("integer.max =", .Machine$integer.max, "\n")
cat("predicted first overflowing neq:",
    min(which(sapply(1:60000, function(n) 22 + n * (n + 9)) > .Machine$integer.max)), "\n\n")

for (neq in c(40000, 46000, 46340, 46341, 46500, 84210)) {
    cat(sprintf("neq=%-6d lrw=%-12.4g (%5.1f GB)  lsoda: %-70s adams: %s\n",
                neq, suppressWarnings(lrw(neq)), suppressWarnings(lrw(neq)) * 8 / 1024^3,
                suppressWarnings(probe(neq, "lsoda")), probe(neq, "adams")))
}

cat("\naugmented system size, nstate = nparms = n, by differentiation order:\n")
aug <- function(k, n) sum(n * (2 * n)^(0:k))
cat(sprintf("%4s %10s %10s %10s %10s\n", "n", "order0", "order1", "order2", "order3"))
for (n in c(2, 4, 6, 8, 9, 10, 12, 16))
    cat(sprintf("%4d %10d %10d %10d %10d\n", n, aug(0, n), aug(1, n), aug(2, n), aug(3, n)))
cat("\norder 1 = a plain gradient; order 2 = a Laplace objective;",
    "order 3 = the gradient of a Laplace objective.\n")
