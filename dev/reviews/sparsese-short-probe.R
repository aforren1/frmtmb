# Short check of lane sparsese: names and shapes the reference scripts read.
.libPaths(c("C:/Users/adf44/source/r/wt-sparsese-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
options(frmtmb.joint_root_max = 20)
q <- function(expr) suppressWarnings(suppressMessages(expr))
set.seed(1)
ng <- 30L; n <- 300L
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(rep(seq_len(ng), length.out = n)))
for (k in 1:5) d[[paste0("x", k)]] <- d$x^k
d$y <- sin(3 * d$x) + stats::rnorm(ng, 0, 0.5)[d$g] +
  stats::rnorm(ng, 0, 0.4)[d$g] * d$x + stats::rnorm(n, 0, 0.3)
fit <- q(frm(y ~ x1 + x2 + x3 + x4 + x5 + (1 + x1 | g), data = d, REML = TRUE))
str(fit$estimates)
print(fit$opt$par)
print(q(VarCorr(fit)))
str(unclass(q(VarCorr(fit))))
print(q(sigma(fit)))
s <- frmtmb:::sdr_of(fit)
print(rownames(s$jointPrecision)[1:12])
lb <- q(frm_lp_basis(fit, newdata = d[1:3, ]))
str(lb)
d$yc <- factor(cut(d$y + stats::rlogis(n), c(-Inf, 0, 0.8, Inf),
                   labels = c("lo", "mid", "hi")), ordered = TRUE)
fc <- q(frm(yc ~ x1 + x2 + x3 + x4 + x5 + (1 | g), data = d,
            family = cumulative(), REML = TRUE))
str(fc$estimates)
print(fc$opt$par)
print(q(summary(fc))$fixed)
print(q(VarCorr(fc)))
print(names(frmtmb:::fit_outer_vector(fc, frmtmb:::fit_draw_space(fc)$map)))
print(frmtmb:::fit_draw_space(fc)$map)
