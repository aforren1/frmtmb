# Reviewer: negbinomial2 (lane nb2dm) on the cdfs CDF, through the
# nb2dm reviewer's trial merge (scratchpad nb2dm-rev-merge/tree, which
# holds this lane's R/families-cdf.R byte for byte), loaded with
# pkgload::load_all(); nothing installed.
.libPaths(c("C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
tree <- paste0("C:/Users/adf44/AppData/Local/Temp/3/claude/",
               "c--Users-adf44-source-r-frmtmb/",
               "7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad/",
               "nb2dm-rev-merge/tree")
suppressMessages(pkgload::load_all(tree, quiet = TRUE, export_all = FALSE))
options(frmtmb.notices = FALSE)
ns <- asNamespace("frmtmb")
cat("frmtmb from", getNamespaceInfo("frmtmb", "path"), "\n")
f2 <- ns$as_frmtmb_family(brmsfamily("negbinomial2"))
cat("slots: lcdf", !is.null(f2$lcdf), " lccdf", !is.null(f2$lccdf),
    " trunc_mean_fn", !is.null(f2$post$trunc_mean_fn), " aterms",
    paste(f2$accepts_aterms, collapse = ","), "\n")
cat("== lccdf / lcdf against stats::pnbinom(size = 1 / sigma)\n")
for (sg in c(1e-8, 1e-3, 0.4, 3, 100)) {
  for (mu in c(0.5, 7, 300)) {
    q <- c(0, 1, floor(mu), floor(3 * mu + 10), floor(20 * mu + 60))
    lS <- f2$lccdf(q, list(mu = mu, sigma = sg), list())
    Fv <- f2$lcdf(q, list(mu = mu, sigma = sg), list())
    rS <- stats::pnbinom(q, size = 1 / sg, mu = mu, lower.tail = FALSE, log.p = TRUE)
    rF <- stats::pnbinom(q, size = 1 / sg, mu = mu)
    cat(sprintf("sigma %-6g mu %-4g: max rel err lS %.2e  F %.2e\n", sg, mu,
                max(abs(lS - rS) / pmax(1, abs(rS))), max(abs(Fv - rF) / pmax(1e-300, rF))))
  }
}
cat("== censored and truncated fits: negbinomial2 vs negbinomial\n")
set.seed(4)
n <- 250
d <- data.frame(x = rnorm(n))
y <- rnbinom(n, mu = exp(1 + 0.4 * d$x), size = 1.7)
d$y <- pmin(y, 9)
d$cc <- as.numeric(y >= 9)
d$cc[sample(which(d$cc == 0), 40)] <- 2
d$y2 <- pmin(d$y + 2, 9)
d$cc[d$cc == 2 & d$y2 <= d$y] <- 0
a <- frm(y | cens(cc, y2) ~ x, family = negbinomial(), data = d)
b <- frm(y | cens(cc, y2) ~ x, family = brmsfamily("negbinomial2"), data = d)
cat(sprintf("cens: logLik nb %.9f nb2 %.9f | shape %.6f 1/sigma %.6f\n",
            logLik(a), logLik(b), exp(a$opt$par[3]), 1 / exp(b$opt$par[3])))
g <- b$obj$gr(b$opt$par + c(0.1, -0.1, 0.3))
fd <- numDeriv::grad(b$obj$fn, b$opt$par + c(0.1, -0.1, 0.3))
cat("nb2 gradient vs FD at a moved point, max rel err:",
    max(abs(g - fd) / pmax(1, abs(fd))), "\n")
d1 <- d[d$y >= 1 & d$cc == 0, ]
a <- frm(y | trunc(lb = 1) ~ x, family = negbinomial(), data = d1)
b <- frm(y | trunc(lb = 1) ~ x, family = brmsfamily("negbinomial2"), data = d1)
cat(sprintf("trunc(lb = 1): logLik nb %.9f nb2 %.9f | fitted max rel diff %.2e\n",
            logLik(a), logLik(b), max(abs(fitted(a)[, 1] / fitted(b)[, 1] - 1))))
cat("== Poisson limit, sigma -> 0: right-censored row at 12, mu 3\n")
for (ls in c(-5, -10, -20, -30)) {
  v <- f2$lccdf(11, list(mu = 3, sigma = exp(ls)), list())
  cat(sprintf("log sigma %g: lccdf(11) %.12f  poisson %.12f\n", ls, v,
              stats::ppois(11, 3, lower.tail = FALSE, log.p = TRUE)))
}
cat("== frm_compat rows on the merge\n")
for (at in c("cens()", "trunc()")) {
  for (fm in c("negbinomial", "negbinomial2")) {
    r <- frm_compat(at, fm)
    cat(sprintf("%-8s %-13s %s\n", at, fm, paste(r$status, collapse = ",")))
  }
}
