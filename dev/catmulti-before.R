# frmtmb's behavior on each item of lane catmulti, on whichever build the
# library path puts first. Arg 1: "base" (rellib-r9) or "lane".
# Output: dev/catmulti-out/before-<arm>.txt
arm <- commandArgs(TRUE)[1]
lib <- "C:/Users/adf44/source/r/wt-catmulti-lib"
base <- "C:/Users/adf44/source/r/rellib-r9"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "lane")) c(lib, base, user) else c(base, user))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
set.seed(1)
n <- 200
dd <- data.frame(x = rnorm(n), z = rnorm(n))
E <- cbind(0, 0.3 + 0.8 * dd$x, -0.4 - 0.5 * dd$x)
P <- exp(E) / rowSums(exp(E))
dd$f <- factor(apply(P, 1, function(p) sample(c("a", "b", "c"), 1, prob = p)))
Yc <- t(sapply(seq_len(n), function(i) rmultinom(1, 10, P[i, ])))
colnames(Yc) <- c("y1", "y2", "y3")
dd$Yc <- Yc
dd$N <- 10L
sec <- function(t) cat("\n========", t, "========\n")
show <- function(expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  print(r)
}
sec("1. categorical(refcat = 'b')")
show(args(categorical))
show(fixef(frm(bf(f ~ x), family = categorical(refcat = "b"), data = dd)))
sec("1. brms::categorical(refcat = 'b')")
show(fixef(frm(bf(f ~ x), family = brms::categorical(refcat = "b"),
               data = dd)))
sec("1. multinomial refcat")
show(fixef(frm(bf(Yc | trials(N) ~ x), family = brms::multinomial(refcat = "y2"),
               data = dd)))
sec("2. bare multinomial()")
show(fixef(frm(bf(Yc | trials(N) ~ x), family = multinomial(), data = dd)))
show(fixef(frm(bf(Yc | trials(N) ~ x), family = brms::multinomial(),
               data = dd)))
sec("2. multinomial(K = 3) names")
fm <- frm(bf(Yc | trials(N) ~ x), family = multinomial(K = 3), data = dd)
show(fixef(fm))
sec("3. fitted() on multinomial")
show(dim(fitted(fm)))
show(dimnames(fitted(fm)))
sec("4. posterior_linpred on categorical fit")
fc <- frm(bf(f ~ x), family = categorical(), data = dd)
show(dim(frmtmb::posterior_linpred(fc, ndraws = 5)))
show(dimnames(frmtmb::posterior_linpred(fc, ndraws = 5)))
show(dim(frmtmb::posterior_epred(fc, ndraws = 5)))
sec("5. mvbf dirichlet / logistic_normal")
G <- matrix(rgamma(n * 3, 3), n, 3); Y <- G / rowSums(G)
colnames(Y) <- c("a", "b", "c"); dd$Y <- Y
dd$w <- rnorm(n)
show(frm(mvbf(bf(Y ~ x, family = dirichlet()),
              bf(w ~ x, family = gaussian())), data = dd))
show(frm(mvbf(bf(Y ~ x, family = logistic_normal()),
              bf(w ~ x, family = gaussian())), data = dd))
