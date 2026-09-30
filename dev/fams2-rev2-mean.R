# Reviewer, punch round 1, item 2 follow-up: RTMB::pbeta()'s third
# derivatives at x == a / (a + b) exactly, how wide the bad set is, and
# whether log_ibeta_half() and an xbeta fit reach it. Then the cost per
# row, timed properly.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
lib <- frmtmb:::log_ibeta_half
G3 <- MakeTape(function(p) log(RTMB::pbeta(p[1], p[2], p[3])), c(0.3, 200, 400))
J1 <- G3$jacfun(); H2 <- J1$jacfun(); T3 <- H2$jacfun()
cat("== 1. RTMB::pbeta at x = a / (a + b) (1 + delta) ==\n")
for (ab in list(c(300, 700), c(268.941421, 731.058579), c(150, 150), c(2e4, 8e4))) {
  a <- ab[1]; b <- ab[2]; mm <- a / (a + b)
  for (dl in c(0, 1e-16, -1e-16, 1e-15, 1e-13, 1e-10, 1e-6)) {
    p <- c(mm * (1 + dl), a, b)
    cat(sprintf("  a %g b %g delta %+.0e: x == mean %s | grad finite %s hess finite %s third finite %s",
                a, b, dl, p[1] == mm, all(is.finite(J1(p))), all(is.finite(H2(p))),
                all(is.finite(T3(p)))))
    if (!all(is.finite(T3(p)))) {
      t3 <- T3(p); cat("  non-finite entries:", which(!is.finite(t3)))
    }
    cat("\n")
  }
}
cat("\n== 2. log_ibeta_half at x = a / (a + b) exactly, inside the pbeta region ==\n")
F <- MakeTape(function(p) lib(p[1], exp(p[2]), exp(p[3])), c(0.1, 0, 0))
T3f <- F$jacfun()$jacfun()$jacfun()
Hf <- F$jacfun()$jacfun()
for (ab in list(c(300, 700), c(2000, 8000), c(2e4, 8e4), c(2e5, 8e5), c(40, 960),
                c(100, 100 / 0.45 - 100))) {
  a <- ab[1]; b <- ab[2]; x <- a / (a + b)
  if (x >= 0.5) next
  p <- c(x, log(a), log(b))
  cat(sprintf("  a %g b %g x %.10g (ab/s %.0f): value %.10g hessian finite %s third finite %s\n",
              a, b, x, a * b / (a + b), F(p), all(is.finite(Hf(p))), all(is.finite(T3f(p)))))
}
cat("\n== 3. an xbeta fit whose q meets mu exactly: kappa held at 0.25 so q = 1/6 ==\n")
suppressPackageStartupMessages(library(frmtmb))
fam <- frmtmb:::as_frmtmb_family(xbeta())
# the xbeta row density at y = 0 with mu = 1/6 exactly and phi = 3000
lp0 <- MakeTape(function(p) {
  sum(fam$lpdf(c(0, 0.3), list(mu = p[1], phi = p[2], kappa = p[3]), list()))
}, c(0.2, 100, 0.1))
T3x <- lp0$jacfun()$jacfun()$jacfun()
for (mu in c(1 / 6, 1 / 6 * (1 + 1e-12))) {
  p <- c(mu, 3000, 0.25)
  cat(sprintf("  mu %.17g, q %.17g: value %.10g third finite %s\n", mu, 0.25 / 1.5,
              lp0(p), all(is.finite(T3x(p)))))
}
cat("\n== 4. cost per row of a gradient sweep (500 rows, seed 5) ==\n")
sp <- paste0("C:/Users/adf44/AppData/Local/Temp/1/claude/",
             "c--Users-adf44-source-r-frmtmb/",
             "66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad")
ex1 <- parse(file.path(sp, "r1", "families.R"))
r1 <- new.env(parent = asNamespace("frmtmb"))
for (e in ex1) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      as.character(e[[2]]) %in% c("log_ibeta_half", "log_ibeta_cf")) eval(e, r1)
}
set.seed(5)
n <- 500
xr <- runif(n, 0.02, 0.3); ar <- exp(runif(n, 0, 6)); br <- exp(runif(n, 0, 6))
mk <- function(f) {
  obj <- MakeADFun(function(p) -sum(f(p$s * xr, p$u * ar, p$v * br)),
                   list(s = 1, u = 1, v = 1), silent = TRUE)
  obj
}
On <- mk(lib); O1 <- mk(r1$log_ibeta_half)
Ob <- mk(function(x, a, b) log(RTMB::pbeta(x, a, b)))
cat("objectives agree: new vs pbeta", On$fn(c(1, 1, 1)) - Ob$fn(c(1, 1, 1)), "\n")
tm <- function(o, reps = 200) {
  best <- Inf
  for (r in 1:7) {
    t0 <- proc.time()[["elapsed"]]
    for (k in seq_len(reps)) o$gr(c(1, 1, 1))
    best <- min(best, (proc.time()[["elapsed"]] - t0) / reps)
  }
  best
}
tt <- c(new = tm(On), round1 = tm(O1), pbeta = tm(Ob), new_again = tm(On))
print(signif(tt * 1e3, 3))
cat(sprintf("gradient per row: new %.1f us, round 1 %.1f us, pbeta %.1f us; new / round 1 %.2f, new / pbeta %.1f\n",
            tt[["new"]] / n * 1e6, tt[["round1"]] / n * 1e6, tt[["pbeta"]] / n * 1e6,
            tt[["new"]] / tt[["round1"]], tt[["new"]] / tt[["pbeta"]]))
