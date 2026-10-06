# Punch round 1: the a + b + log(c0) ridge (test-se-check.R): what the
# projection rule sees on sigma.
.libPaths(c("C:/Users/adf44/source/r/wt-nanse-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(45)
d <- data.frame(x = rnorm(200))
d$y <- 1 + 0.5 * d$x + log(5e-5) + rnorm(200, 0, 0.3)
f <- suppressWarnings(frm(bf(y ~ a + b + log(c0), a ~ 1 + x, b ~ 1 + x,
                             c0 ~ 1, nl = TRUE),
                          data = d, start = list(beta = c(0.5, 0.25, 0.5,
                                                          0.25, 5e-5))))
p <- f$opt$par
Hx <- f$obj$he(p)
E <- abs(Hx - t(Hx)) / 2
H <- (Hx + t(Hx)) / 2
D <- sqrt(abs(diag(H)))
S <- H / outer(D, D)
e <- eigen(S, TRUE)
cat("eig:", signif(e$values, 3), "\n")
cat("enorm:", sqrt(sum((E / outer(D, D))^2)), "\n")
fl <- which(e$values < 1e-9 * max(e$values))
print(round(e$vectors[, fl], 8))
cat("projection:", signif(sqrt(rowSums(e$vectors[, fl]^2)), 3), "\n")
