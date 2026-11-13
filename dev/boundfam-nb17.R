# Lane boundfam fix round, m1: why the negbinomial at log shape 23.2
# (seed 17) is not marked by se_edge_limit().
.libPaths(c("C:/Users/adf44/source/r/wt-boundfam-lib", "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(17)
d <- data.frame(x = rnorm(200))
d$y <- rpois(200, exp(1 + 0.4 * d$x))
f <- suppressWarnings(frm(y ~ x, family = negbinomial(), data = d))
p <- f$opt$par
H <- f$obj$he(p)
cat("diag H:", diag(H), "\n")
cat("par_units:", f$par_units, "\n")
cat("se_explained:", f$cache$se_explained, "\n")
lim <- ns$se_edge_limit(f, H)
cat("se_edge_limit:", lim, "\n")
str(ns$se_limit_info(f))
an <- ns$cov_from_hessian(f, (H + t(H)) / 2, abs(H - t(H)) / 2, exact = TRUE)
print(an$lost)
s <- suppressWarnings(ns$sdr_of(f))
print(s$se_lost)
