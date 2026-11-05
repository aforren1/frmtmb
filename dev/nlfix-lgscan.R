# R3: the other discrete families with a precision that runs off on data with no
# overdispersion. Is the objective flat at the limit, as the true density is?
# Arg 1 is the library. A fit and a grid of the precision at the limit's estimate.
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(getNamespaceInfo("frmtmb", "path")), "\n")
ns <- asNamespace("frmtmb")
mk <- function(s, n = 300) {
  set.seed(s)
  d <- data.frame(x = round(rnorm(n), 3), ntr = 10L)
  d$yp <- rpois(n, exp(0.5 + 0.5 * d$x))
  d$yb <- rbinom(n, 10, plogis(0.3 + 0.5 * d$x))
  d
}
cases <- list(
  list("negbinomial", negbinomial(), bf(yp ~ x)),
  list("nbinom1", nbinom1(), bf(yp ~ x)),
  list("zi_nb", zero_inflated_negbinomial(), bf(yp ~ x, zi ~ 1)),
  list("hu_nb", hurdle_negbinomial(), bf(yp ~ x, hu ~ 1)),
  list("beta_binomial", beta_binomial(), bf(yb | trials(ntr) ~ x)),
  list("zibb", zero_inflated_beta_binomial(), bf(yb | trials(ntr) ~ x, zi ~ 1)))
for (cs in cases) {
  cat("\n---", cs[[1]], "\n")
  for (s in c(2, 4, 6, 9, 12)) {
    d <- mk(s)
    w <- character(0)
    f <- tryCatch(withCallingHandlers(
      suppressMessages(frm(cs[[3]], family = cs[[2]], data = d)),
      warning = function(c) {
        w <<- c(w, substr(conditionMessage(c), 1, 40))
        invokeRestart("muffleWarning")
      }), error = function(e) e)
    if (inherits(f, "error")) {
      cat(sprintf("seed %2d ERROR %s\n", s, substr(conditionMessage(f), 1, 90)))
      next
    }
    cat(sprintf("seed %2d ll %.6f code %d betad %s warn[%s]\n", s,
                as.numeric(logLik(f)), f$opt$convergence,
                paste(round(f$estimates$betad, 2), collapse = ","),
                paste(unique(w), collapse = " | ")))
    if (s != 9) next
    nll <- ns$build_objective(f$frame)
    tp <- f$estimates
    jd <- length(tp$betad)
    if (grepl("zi|hu", cs[[1]])) jd <- 1L
    for (lp in c(10, 15, 20, 25, 30, 35, 39, 45, 60, 100)) {
      t1 <- tp
      t1$betad[jd] <- if (identical(cs[[1]], "nbinom1")) -lp else lp
      v <- tryCatch(-as.numeric(nll(t1)), error = function(e) NA)
      cat(sprintf("  log precision %3d: %.8g\n", lp, v))
    }
  }
}
cat("\n--- multinomial data, dirichlet_multinomial(), the objective over log phi\n")
set.seed(9)
n <- 300
d <- data.frame(x = round(rnorm(n), 3))
d$y <- t(vapply(seq_len(n), function(i) {
  rmultinom(1, 10, c(1, exp(0.3 * d$x[i]), exp(-0.2 * d$x[i])))[, 1]
}, numeric(3)))
colnames(d$y) <- c("a", "b", "c")
d$ntr <- 10L
f <- tryCatch(suppressWarnings(suppressMessages(
  frm(bf(y | trials(ntr) ~ x), family = dirichlet_multinomial(), data = d))),
  error = function(e) e)
if (inherits(f, "error")) {
  cat("ERROR", conditionMessage(f), "\n")
} else {
  cat(sprintf("ll %.6f code %d betad %s\n", as.numeric(logLik(f)),
              f$opt$convergence, paste(round(f$estimates$betad, 2), collapse = ",")))
  nll <- ns$build_objective(f$frame)
  for (lp in c(10, 20, 30, 39, 45, 60, 100)) {
    t1 <- f$estimates
    t1$betad[length(t1$betad)] <- lp
    cat(sprintf("  log phi %3d: %.8f\n", lp, -as.numeric(nll(t1))))
  }
}
cat("\n--- continuous, for the record: Beta() and its inflated forms, tape against lbeta_ad\n")
for (lp in c(10, 20, 30, 39, 45)) {
  phi <- exp(lp)
  a <- 0.4 * phi
  b <- 0.6 * phi
  cat(sprintf("  log phi %2d: raw -lbeta %.6f, stable %.6f\n", lp,
              lgamma(a + b) - lgamma(a) - lgamma(b), -ns$lbeta_ad(a, b)))
}
cat("\n--- simulate: rbetabinom at large phi and at a mean at the edge\n")
d <- mk(9)
f <- suppressWarnings(suppressMessages(
  frm(bf(yb | trials(ntr) ~ x), family = beta_binomial(), data = d)))
sim <- tryCatch(simulate(f, nsim = 2, seed = 1), error = function(e) e)
if (inherits(sim, "error")) cat("simulate ERROR", conditionMessage(sim), "\n") else {
  print(summary(as.numeric(unlist(sim))))
}
fam <- f$spec$responses[[1]]$family
set.seed(1)
for (p in list(c(0.5, exp(39)), c(0.5, exp(60)), c(1, exp(5)), c(0, exp(5)))) {
  r <- tryCatch(fam$sim(list(mu = rep(p[1], 5), phi = rep(p[2], 5)),
                        list(trials = rep(10L, 5)), 5L), error = function(e) conditionMessage(e))
  cat(" mu", p[1], "phi", format(p[2], digits = 3), ":", format(r), "\n")
}
