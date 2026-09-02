# Run the ?frmtmb-student-re examples exactly as R CMD check would.
# Run: Rscript dev/tre/check-examples.R
sink("dev/tre/check-examples.txt")
suppressMessages(pkgload::load_all(".", quiet = TRUE))

set.seed(1)
n <- 12
d <- data.frame(x = rnorm(20 * n), g = factor(rep(1:20, each = n)))
b <- rnorm(20)
b[20] <- b[20] + 6          # one outlying group
d$y <- 1 + 0.5 * d$x + b[d$g] + rnorm(20 * n)

fit_t <- frm(bf(y ~ x + (1 | gr(g, dist = "student"))),
             family = gaussian(), data = d)
fit_n <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)

print(VarCorr(fit_t))
print(VarCorr(fit_n))

print(frm(bf(y ~ x + (1 | gr(g, dist = "student", dist_nu = 3))),
          family = gaussian(), data = d))
sink()
cat("done\n")
