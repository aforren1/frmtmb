# test-nl-start.R's binomial Michaelis-Menten on a ridge to the linear
# limit (seed 22001): what each build says at the default stop.
#   CODE8_ARM=before|after Rscript dev/code8-mm22.R
source("dev/code8-common.R")
set.seed(22001)
n <- 300
d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3),
                w = round(rnorm(n), 3), xp = round(runif(n, 0.1, 5), 3),
                xq = round(runif(n, 0, 5), 3),
                xb = round(runif(n, 0, 3000), 1),
                g = factor(sample(1:15, n, TRUE)))
d$ug <- rnorm(15, 0, 0.4)[d$g]
d$nt <- 10L
d$y <- rbinom(n, 10, plogis(3 * d$xp / (1 + d$xp) - 2))
fm <- bf(y | trials(nt) ~ Vm * xp / (K + xp), Vm ~ 1, K ~ 1, nl = TRUE)
r <- rfit(fm, family = binomial(), data = d)
f <- r$f
cat("objective", format(f$opt$objective, digits = 10), "code",
    f$opt$convergence, f$opt$message, "\n par", format(f$opt$par), "\n")
cat(gsub(" || ", "\n", cond_txt(r), fixed = TRUE), "\n")
g <- drop(f$obj$gr(f$opt$par))
cat("grad", format(g), "\n")
cat("grad_verdict:\n"); str(fns$grad_verdict(f, g, f$control))
sdr <- quiet(fns$sdr_of(f))
cat("se_lost:", paste(names(sdr$se_lost), sdr$se_lost), "\n")
cat("se:", format(sqrt(diag(sdr$cov.fixed))), "\n")
cat("se_explained:", format(f$cache$se_explained), "\n")
