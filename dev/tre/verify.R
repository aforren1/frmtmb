# Direct checks of the registry entries against independent references,
# separated from the optimizer so a flat direction cannot be mistaken
# for a wrong density.
# Run: Rscript dev/tre/verify.R
sink("dev/tre/verify.txt")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressMessages(library(RTMB))

mk_blk <- function(cs, d, nlev, nu = NULL) {
  list(covstruct = cs, dim = d, n_levels = nlev, dist_nu = nu,
       cnms = paste0("c", seq_len(d)))
}
reg <- frmtmb:::covstruct_registry

cat("1. us_t / diag_t against mvtnorm's dmvt at matched theta.\n\n")
set.seed(1)
for (d in c(1, 2, 3)) {
  for (nu in c(2.5, 3, 5, 30)) {
    nlev <- 7
    th <- if (d == 1) 0.2 else
      c(seq(-0.3, 0.3, length.out = d), rep(0.35, d * (d - 1) / 2))
    b <- rnorm(d * nlev)
    blk <- mk_blk("us_t", d, nlev, nu)
    got <- reg$us_t$nll(b, th, blk)
    S <- reg$us_t$vcov(th, blk)
    B <- matrix(b, nrow = d)
    want <- sum(mvtnorm::dmvt(t(B), sigma = unname(S), df = nu,
                              log = TRUE))
    cat(sprintf("us_t   d %d nu %-5s  %14.9f  %14.9f  diff %.2e\n",
                d, nu, got, want, abs(got - want)))
    blk2 <- mk_blk("diag_t", d, nlev, nu)
    th2 <- seq(-0.3, 0.3, length.out = d)
    got2 <- reg$diag_t$nll(b, th2, blk2)
    S2 <- reg$diag_t$vcov(th2, blk2)
    want2 <- sum(mvtnorm::dmvt(t(B), sigma = unname(S2), df = nu,
                               log = TRUE))
    cat(sprintf("diag_t d %d nu %-5s  %14.9f  %14.9f  diff %.2e\n",
                d, nu, got2, want2, abs(got2 - want2)))
  }
}

cat("\n2. The gaussian limit of the DENSITY (no optimizer involved).\n\n")
cat(sprintf("%-8s %-4s %14s %14s\n", "nu", "d", "|us_t - us|",
            "|diag_t - diag|"))
for (nu in c(1e4, 1e6, 1e8, 1e10)) {
  for (d in c(1, 2, 3)) {
    nlev <- 20
    b <- rnorm(d * nlev)
    th <- if (d == 1) 0.2 else
      c(seq(-0.3, 0.3, length.out = d), rep(0.35, d * (d - 1) / 2))
    th2 <- seq(-0.3, 0.3, length.out = d)
    a1 <- reg$us_t$nll(b, th, mk_blk("us_t", d, nlev, nu))
    a0 <- reg$us$nll(b, th, mk_blk("us", d, nlev))
    c1 <- reg$diag_t$nll(b, th2, mk_blk("diag_t", d, nlev, nu))
    c0 <- reg$diag$nll(b, th2, mk_blk("diag", d, nlev))
    cat(sprintf("%-8.0e %-4d %14.3e %14.3e\n", nu, d, abs(a1 - a0),
                abs(c1 - c0)))
  }
}

cat("\n3. AD gradient of each entry against finite differences.\n\n")
for (cs in c("us_t", "diag_t")) {
  for (d in c(1, 2, 3)) {
    for (nu in c(2.5, 5, 30)) {
      nlev <- 6
      np <- reg[[cs]]$npar(d)
      th0 <- if (np == 1) 0.15 else
        c(seq(-0.2, 0.2, length.out = d),
          rep(0.3, np - d))[seq_len(np)]
      b0 <- rnorm(d * nlev)
      blk <- mk_blk(cs, d, nlev, nu)
      f <- function(p) {
        reg[[cs]]$nll(p[seq_len(d * nlev)],
                      p[d * nlev + seq_len(np)], blk)
      }
      p0 <- c(b0, th0)
      tp <- MakeTape(f, p0)
      cat(sprintf("%-7s d %d nu %-5s  max|AD - FD| %.2e\n", cs, d, nu,
                  max(abs(tp$jacobian(p0) -
                            numDeriv::grad(f, p0)))))
    }
  }
}

cat("\n4. draw_b(): does simulate() draw a multivariate t?\n\n")
set.seed(2)
G <- 30; n <- 6
dd <- data.frame(x = rnorm(G * n), g = factor(rep(1:G, each = n)))
dd$y <- 1 + 0.5 * dd$x + rt(G, 5)[dd$g] + rnorm(G * n)
f <- frm(bf(y ~ x + (1 | gr(g, dist = "student", dist_nu = 4))) +
           gaussian(), data = dd)
S <- unname(VarCorr(f)[[1]])
B <- replicate(20000, frmtmb:::draw_b(f))
q <- as.vector(B^2) / S[1, 1]
cat("KS of b^2/scale^2 against F(1, 4): ")
k <- suppressWarnings(ks.test(q, "pf", 1, 4))
cat(sprintf("D = %.5f, p = %.4f\n", k$statistic, k$p.value))
cat(sprintf("sample var %.4f, theory scale^2*nu/(nu-2) = %.4f\n",
            var(as.vector(B)), S[1, 1] * 4 / 2))

fg <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
Bg <- replicate(20000, frmtmb:::draw_b(fg))
kg <- suppressWarnings(ks.test(as.vector(Bg) /
                                 sqrt(unname(VarCorr(fg)[[1]])[1, 1]),
                               "pnorm"))
cat(sprintf("gaussian control: KS D = %.5f, p = %.4f\n", kg$statistic,
            kg$p.value))

cat("\n5. Correlated draw: is the mixing variable shared per level?\n\n")
dd$z <- rnorm(G * n)
dd$y2 <- 1 + dd$x + rt(G, 5)[dd$g] + rnorm(G * n)
f2 <- frm(bf(y2 ~ x + (x | gr(g, dist = "student", dist_nu = 5))) +
            gaussian(), data = dd)
S2 <- unname(VarCorr(f2)[[1]])
B2 <- t(replicate(20000, matrix(frmtmb:::draw_b(f2), nrow = 2)[, 1]))
qq <- rowSums((B2 %*% solve(S2)) * B2)
k2 <- suppressWarnings(ks.test(qq / 2, "pf", 2, 5))
cat(sprintf("KS of q/d against F(2, 5): D = %.5f, p = %.4f\n",
            k2$statistic, k2$p.value))
cat(sprintf("empirical cor %.4f, scale cor %.4f\n",
            cor(B2[, 1], B2[, 2]), stats::cov2cor(S2)[1, 2]))

cat("\n6. Unknown structure name around a bar (pre-existing?).\n")
gg <- function(e) tryCatch({ eval(e); "NO ERROR" },
                           error = function(x) conditionMessage(x))
cat("  foo(x | g) : ",
    gg(quote(frm(bf(y ~ foo(x | g)) + gaussian(), data = dd))), "\n")
cat("  us_t(x | g): ",
    gg(quote(frm(bf(y ~ us_t(x | g)) + gaussian(), data = dd))), "\n")

sink()
cat("done\n")
