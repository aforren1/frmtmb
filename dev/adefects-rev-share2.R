insp1 <- function(x) utils::capture.output(.Internal(inspect(x)))[1L]
set.seed(1)
mf <- data.frame(a = rnorm(20000), b = rnorm(20000))
l1 <- list(x = mf, y = mf)
td <- tempdir(); p <- file.path(td, "l1b.rds")
saveRDS(l1, p, compress = FALSE)
r <- readRDS(p)
ax <- insp1(r$x); ay <- insp1(r$y)
cat("mem   x:", insp1(l1$x), "\n")
cat("mem   y:", insp1(l1$y), "\n")
cat("disk  x:", ax, "\n")
cat("disk  y:", ay, "\n")
# a mutation test, which does not depend on reading an address at all
r$x$a[1] <- -999
cat("after mutating r$x$a[1]: r$y$a[1] =", r$y$a[1], "\n")
l1$x$a[1] <- -999
cat("in-memory copy-on-modify: l1$y$a[1] =", l1$y$a[1], "\n")
