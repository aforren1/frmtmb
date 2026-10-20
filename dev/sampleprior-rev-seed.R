# Reviewer of lane sampleprior: frm_sample(seed = ), serial and parallel
# chains, the caller's .Random.seed, and no seed.
lane <- !identical(Sys.getenv("REV_ARM"), "r9")
.libPaths(c(if (lane) "C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb.sample"), "\n")
set.seed(20261117)
d <- data.frame(x = rnorm(40), g = factor(rep(1:5, 8)))
d$y <- 1 + d$x + rnorm(5)[d$g] + rnorm(40)
fit <- frm(y ~ x + (1 | g), data = d)
run <- function(..., m = fit) {
  suppressMessages(suppressWarnings(frm_sample(
    m, chains = 3, iter = 400, refresh = 0, ...)))$draws
}
set.seed(1); a <- run(seed = 3)
runif(5)
b <- run(seed = 3)
cat("serial seed 3 twice identical:", identical(a, b), "\n")
c4 <- run(seed = 4)
cat("seed 4 differs from seed 3 (chains 1, 2, 3):",
    vapply(1:3, function(k) !identical(a[(k - 1) * 200 + 1:200, ],
                                       c4[(k - 1) * 200 + 1:200, ]), NA), "\n")
p1 <- run(seed = 3, cores = 3)
p2 <- run(seed = 3, cores = 3)
cat("parallel seed 3 twice identical:", identical(p1, p2),
    " parallel equals serial:", identical(p1, a),
    " max |parallel - serial|:", max(abs(p1 - a)), "\n")
# the caller's stream
set.seed(77); s0 <- .Random.seed
invisible(run(seed = 3))
cat("caller .Random.seed unchanged after frm_sample(seed = 3):",
    identical(s0, .Random.seed), "\n")
# what moves it: rstan itself (sampling() draws with sample.int)
set.seed(77)
fs <- asNamespace("frmtmb.sample")
if (exists("with_sample_seed", fs)) {
  invisible(fs$with_sample_seed(3, runif(3)))
  cat("with_sample_seed alone leaves it unchanged:",
      identical(s0, .Random.seed), "\n")
  rm(.Random.seed, envir = globalenv())
  invisible(fs$with_sample_seed(3, runif(3)))
  cat("no stream before, none after with_sample_seed:",
      !exists(".Random.seed", envir = globalenv()), "\n")
  set.seed(77)
  RNGkind("L'Ecuyer-CMRG"); set.seed(77); s1 <- .Random.seed
  invisible(fs$with_sample_seed(3, runif(3)))
  cat("L'Ecuyer kind restored:", identical(s1, .Random.seed),
      RNGkind()[1], "\n")
  RNGkind("Mersenne-Twister")
}
# no seed: set.seed() decides, as before
set.seed(9); n1 <- run()
set.seed(9); n2 <- run()
n3 <- run()
cat("no seed, set.seed(9) twice identical:", identical(n1, n2),
    " a third call differs:", !identical(n1, n3), "\n")
# formula route: random inits drawn by Stan from the seed
f1 <- run(seed = 3, m = y ~ x + (1 | g), data = d)
f2 <- run(seed = 3, m = y ~ x + (1 | g), data = d)
cat("formula route seed 3 twice identical:", identical(f1, f2), "\n")
