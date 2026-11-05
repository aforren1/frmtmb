# Fix round 1: why does importance on an uncentered ordinal land 3.4 nats below the
# Laplace fit on the same column, without a refusal? Data seed 2 of own2.R, m = 2000.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs(Sys.getenv("REV_ARM", "after"))
suppressPackageStartupMessages(library(frmtmb))
options(width = 170, digits = 7)
s <- 2L
m <- 2000
set.seed(s)
G <- 15L
g <- factor(rep(seq_len(G), each = 10))
z <- rnorm(150)
lat <- 0.7 * z + rnorm(G, 0, 0.7)[g] + rlogis(150)
y <- ordered(as.integer(cut(lat, c(-Inf, -1, 0.5, Inf))))
d <- data.frame(y, w = 5 * z + m, wc = 5 * z, g)
fit <- frm(bf(y ~ w + (1 | g)), family = cumulative(), data = d, importance = 300L, verbose = TRUE)
print(fit$opt$par)
print(fit$importance[c("rounds", "moves", "grads", "capped", "start_value", "ess_min", "mcse")])
print(as.numeric(logLik(fit)))
lap <- suppressWarnings(frm(bf(y ~ w + (1 | g)), family = cumulative(), data = d, verbose = TRUE))
print(lap$opt$par)
print(as.numeric(logLik(lap)))
cen <- suppressWarnings(frm(bf(y ~ wc + (1 | g)), family = cumulative(), data = d, importance = 300L))
print(cen$opt$par)
print(as.numeric(logLik(cen)))
print(names(fit$frame$par_template))
# the plain tape at the fit's estimates, Laplace objective at the two points
ns <- asNamespace("frmtmb")
fr <- fit$frame
nll <- ns$build_objective(fr)
tpl <- ns$make_start(fr, NULL, NULL)
l2 <- RTMB::MakeADFun(nll, tpl, random = "b", map = fr[["map"]], silent = TRUE)
cat("start of the outer parameters:", format(l2$par, digits = 5), "\n")
cat("Laplace objective at importance's estimate", l2$fn(fit$opt$par), "\n")
cat("Laplace objective at the Laplace fit's", l2$fn(lap$opt$par), "\n")
o <- nlminb(l2$par, l2$fn, l2$gr, control = list(iter.max = 1000, eval.max = 1000))
cat("plain nlminb from the start:", format(o$par, digits = 6), "obj", o$objective, "code", o$convergence, o$message, "\n")
print(l2$gr(o$par))
