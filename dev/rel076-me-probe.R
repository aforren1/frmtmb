# The marginaleffects test of test-mo-coding.R: how large is the finite
# difference error of marginaleffects' Jacobian, measured two ways, so
# the release fix can bound the error by a quantity the run measures
# rather than a fixed 1e-4. Reads the data helper from the test file.
#   Rscript dev/rel076-me-probe.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
src <- readLines("tests/testthat/test-mo-coding.R")
s <- grep("^mo_code_data <- function", src)
e <- s + which(src[(s + 1):length(src)] == "}")[1]
eval(parse(text = src[s:e]))
d <- mo_code_data()
fit <- frm(y ~ mo(k) + z, data = d)
m <- stats::lm(y ~ factor(k) + z, data = d)
nd <- data.frame(k = 1:4, z = 0)
ex <- stats::predict(m, newdata = nd, se.fit = TRUE)$se.fit *
  sqrt(m$df.residual / nrow(d))
P <- function(nm) marginaleffects::predictions(fit, newdata = nd,
                                               numderiv = nm)$std.error
f <- marginaleffects::predictions(fit, newdata = nd)$std.error
fc <- P("fdcenter")
fr <- P("richardson")
print(rbind(ex = ex, fwd = f, ctr = fc, rich = fr), digits = 12)
cat("fwd/ex-1", format(abs(f / ex - 1)), "\n")
cat("ctr/ex-1", format(abs(fc / ex - 1)), "\n")
cat("rich/ex-1", format(abs(fr / ex - 1)), "\n")
cat("fwd/ctr-1", format(abs(f / fc - 1)), "\n")
V <- marginaleffects::get_vcov(fit)
print(colnames(V))
sx <- grep("^zeta", colnames(V))
V0 <- V; V0[sx, ] <- 0; V0[, sx] <- 0
se0 <- marginaleffects::predictions(fit, newdata = nd, vcov = V0)$std.error
cat("se0", format(se0, digits = 10), "\n")
cat("se0/ex-1", format(abs(se0 / ex - 1)), "\n")
cat("ratio", max(abs(f / ex - 1)) / max(abs(se0 / ex - 1)), "\n")
cs <- frm_linpred(fit, newdata = nd, se.fit = TRUE)$se.fit
cat("core/ex-1", format(abs(cs / ex - 1)), "\n")
cat("fwd/core-1", format(abs(f / cs - 1)), "\n")
cat("ctr/core-1", format(abs(fc / cs - 1)), "\n")
sml <- sqrt(sum(stats::residuals(m)^2) / nrow(d))
cat("sigma/sigmaML-1", format(abs(stats::sigma(fit) / sml - 1)), "\n")
pl <- stats::predict(m, newdata = nd)
pf <- frm_linpred(fit, newdata = nd)
cat("pred rel", format(abs(pf - pl) / max(abs(pl))), "\n")
lb <- frmtmb:::frm_lp_basis(fit, newdata = nd, re_formula = NA)
jc <- frmtmb:::get_joint_cov(fit)
cat("root null:", is.null(jc$root), " Vrest null:", is.null(jc$Vrest), "\n")
A <- as.matrix(lb$A); Vb <- as.matrix(jc$V)[lb$coef_pos, lb$coef_pos]
zc <- which(jc$names[lb$coef_pos] == "zeta1")
sfull <- sqrt(rowSums((A %*% Vb) * A))
kp <- setdiff(seq_len(ncol(A)), zc)
slost <- sqrt(rowSums((A[, kp, drop = FALSE] %*% Vb[kp, kp]) * A[, kp, drop = FALSE]))
cat("sfull/core-1", format(abs(sfull / cs - 1)), "\n")
cat("slost", format(slost, digits = 8), "\n")
cat("slost/ex-1", format(abs(slost / ex - 1)), "\n")
