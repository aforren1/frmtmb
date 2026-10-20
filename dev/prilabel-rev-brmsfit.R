# Review of lane prilabel: two tiny brms 2.23.0 fits on cases the lane
# did not fit, and the lane's frmtmb on the same data, for the names a
# fit carries (variables(), ranef(), VarCorr(), ngrps()) and re_formula.
# Data: tests/testthat/test-group-labels.R lab_data() (seed 5) plus
# columns drawn with seed 12. brms: 1 chain, 200 iterations, seed 1;
# the draws are used for names only, except the re_formula contrasts.
.libPaths(c("C:/Users/adf44/source/r/wt-prilabel-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb)
})
lg <- file.path(wt, "dev/prilabel-rev-log")
sink(file.path(lg, "brmsfit.txt"), split = TRUE)
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "; brms",
    format(packageVersion("brms")), "\n")
lab_data <- function(n = 300, seed = 5) {
  set.seed(seed)
  g1 <- factor(sample(sprintf("p%d", 1:6), n, TRUE))
  g2 <- factor(sample(sprintf("q%d", 1:5), n, TRUE))
  m1 <- factor(sample(sprintf("m%02d", 1:10), n, TRUE))
  m2 <- factor(sample(sprintf("m%02d", 1:10), n, TRUE))
  g <- factor(sample(sprintf("g%02d", 1:12), n, TRUE))
  x <- rnorm(n)
  y <- 1 + 0.5 * x + rnorm(6, 0, 0.7)[as.integer(g1)] +
    rnorm(30, 0, 0.5)[(as.integer(g1) - 1L) * 5L + as.integer(g2)] +
    0.5 * (rnorm(10, 0, 0.5)[as.integer(m1)] +
             rnorm(10, 0, 0.5)[as.integer(m2)]) + rnorm(n, 0, 0.6)
  data.frame(y, x, g1, g2, m1, m2, g)
}
dd <- lab_data()
set.seed(12)
n <- nrow(dd)
dd$g3 <- factor(sample(c("r1", "r2"), n, TRUE))
dd$f <- factor(sample(c("a", "b", "c"), n, TRUE))
dd$z <- rnorm(n)
dd$yz <- dd$y + rnorm(n)
dd$w1 <- runif(n)
dd$w2 <- 1 - dd$w1
dd$g1b <- factor(ifelse(as.integer(dd$g1) <= 3L, "lo", "hi"))
# levels with a ':' (no '_' collision)
dd$k1 <- factor(sample(c("c:d", "e"), n, TRUE))
dd$k2 <- factor(sample(c("f", "h"), n, TRUE))
dd$y <- dd$y + c(0.6, -0.6)[as.integer(dd$k1)] * c(1, -1)[as.integer(dd$k2)]
try_ <- function(expr) tryCatch(expr, error = function(e) {
  paste("ERROR:", conditionMessage(e))
})
names_of <- function(fit) {
  v <- variables(fit)
  re <- ranef(fit)
  list(sdcor = grep("^(sd|cor)_", v, value = TRUE),
       ranef = names(re),
       levels = lapply(re, function(a) dimnames(a)[[1L]]),
       coefs = lapply(re, function(a) {
         dn <- dimnames(a)
         dn[[length(dn)]]
       }),
       vc = names(VarCorr(fit)),
       ngrps = names(ngrps(fit)))
}
cmp <- function(a, b, what) {
  same <- identical(a, b)
  cat(sprintf("%-8s %s\n", what, if (same) "SAME" else "DIFF"))
  if (!same) {
    cat(" brms:\n"); str(a)
    cat(" frmtmb:\n"); str(b)
  }
}
fA <- y ~ x + (1 + x | g1/g2/g3) + (1 + f || g) +
  (1 | mm(m1, m2, weights = cbind(w1, w2))) + (1 | k1/k2)
bA <- brm(brms::bf(fA), data = dd, chains = 1, iter = 200, warmup = 100,
          seed = 1, refresh = 0)
rA <- frm(frmtmb::bf(fA), data = dd)
cat("\n==== fit A:", deparse1(fA), "\n")
nb <- names_of(bA)
nf <- names_of(rA)
for (k in names(nb)) cmp(nb[[k]], nf[[k]], k)

cat("\n-- re_formula on fit A: brms posterior-mean fitted, 4 rows\n")
nd <- dd[1:4, ]
rfs <- list(NA, ~ (1 + x | g1:g2), ~ (1 + x | g2:g1), ~ (1 + x | g1),
            ~ (1 + x | g1/g2), ~ (1 | k1:k2))
for (rf in rfs) {
  pb <- try_(round(fitted(bA, newdata = nd, re_formula = rf)[, 1], 3))
  pf <- try_(round(fitted(rA, newdata = nd, re_formula = rf)[, 1], 3))
  cat(format(if (is.na(rf)[1]) "NA" else deparse1(rf), width = 22),
      "brms:", pb, " frmtmb:", pf, "\n")
}
cat("\n-- new level messages\n")
nd2 <- dd[1:2, ]
nd2$k2 <- factor("zz")
cat("brms:  ", try_(fitted(bA, newdata = nd2)), "\n")
cat("frmtmb:", try_(fitted(rA, newdata = nd2)), "\n")

fB <- brms::bf(y ~ x + (1 + x | p | g) + (1 + x | gr(g1, by = g1b))) +
  brms::bf(yz ~ z + (1 | p | g)) + brms::set_rescor(FALSE)
bB <- brm(fB, data = dd, chains = 1, iter = 200, warmup = 100, seed = 1,
          refresh = 0)
rB <- frm(frmtmb::bf(y ~ x + (1 + x | p | g) + (1 + x | gr(g1, by = g1b))) +
            frmtmb::bf(yz ~ z + (1 | p | g)) + frmtmb::set_rescor(FALSE),
          data = dd)
cat("\n==== fit B: mvbf with |p| across responses and gr(by = )\n")
nb <- names_of(bB)
nf <- names_of(rB)
for (k in names(nb)) cmp(nb[[k]], nf[[k]], k)
cat("\n-- per-coef sd prior on fit B, brms prior_summary against frmtmb\n")
pr_b <- brms::set_prior("normal(0, 0.5)", class = "sd", group = "g1",
                        coef = "x", resp = "y")
pr_f <- frmtmb::set_prior("normal(0, 0.5)", class = "sd", group = "g1",
                          coef = "x", resp = "y")
ps_b <- as.data.frame(brms::validate_prior(pr_b, fB, data = dd))
ps_f <- as.data.frame(frmtmb::validate_prior(
  pr_f, frmtmb::bf(y ~ x + (1 + x | p | g) + (1 + x | gr(g1, by = g1b))) +
    frmtmb::bf(yz ~ z + (1 | p | g)) + frmtmb::set_rescor(FALSE), data = dd))
keep <- c("prior", "class", "coef", "group", "resp", "source")
print(ps_b[ps_b$class == "sd", keep], row.names = FALSE)
print(ps_f[ps_f$class == "sd", keep], row.names = FALSE)
sink()
