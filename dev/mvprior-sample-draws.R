# Lane wt-mvprior: the punch-round-2 reviewer's dev/mvprior-review2/r2-sample.R,
# rerun by the lane on its final build, saving under dev/mvprior-log.
# Short frm_sample() runs per arm on categorical, mixture, lca and a
# mixture with partial dpar priors. Saves prior_summary, draws column
# names and draws, for identical() across arms.
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-review2/r2-prelude.R")
suppressPackageStartupMessages({library(frmtmb.sample); library(frmtmb.latent)})
cat("StanHeaders", format(packageVersion("StanHeaders")), " rstan", format(packageVersion("rstan")), "\n")
source(file.path(R2_ROOT, "dev/mvprior-review2/r2-cases.R"))
d <- r2_data()
res <- list()
go <- function(nm, f, fam, dat, prior = NULL) {
  t0 <- proc.time()[["elapsed"]]
  msgs <- character(0)
  s <- tryCatch(withCallingHandlers(
    suppressWarnings(frm_sample(f, data = dat, family = fam, prior = prior,
                                chains = 1, iter = 300, warmup = 150, seed = 11,
                                refresh = 0)),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }),
    error = function(e) e)
  if (inherits(s, "error")) { cat(nm, "ERROR", conditionMessage(s), "\n"); res[[nm]] <<- conditionMessage(s); return() }
  dr <- as.matrix(s)
  ps <- tryCatch(as.data.frame(prior_summary(s)), error = function(e) conditionMessage(e))
  cat("\n==", nm, " ok in", round(proc.time()[["elapsed"]] - t0, 1), "s; draws", nrow(dr), "x", ncol(dr), "\n")
  cat("announce:\n", paste(msgs, collapse = ""), "\n")
  cat("columns:", paste(head(colnames(dr), 14), collapse = " "), "\n")
  if (is.data.frame(ps)) print(ps[, intersect(c("prior", "class", "coef", "group", "resp", "dpar", "source"), names(ps))], row.names = FALSE)
  res[[nm]] <<- list(cols = colnames(dr), draws = dr, ps = ps, msgs = msgs)
}
go("cat4re", bf(cat4 ~ x + (1 | g)), categorical(), d)
go("mix2", bf(ym ~ x), mixture(gaussian(), gaussian()), d)
go("mix2_dpar", bf(ym ~ x), mixture(gaussian(), gaussian()), d,
   set_prior("normal(0, 1)", class = "b", dpar = "mu1") +
     set_prior("normal(-2, 2)", class = "Intercept", dpar = "mu1") +
     set_prior("normal(3, 2)", class = "Intercept", dpar = "mu2") +
     set_prior("student_t(3, 0, 2.5)", class = "sigma1"))
go("mix2_partial", bf(ym ~ x), mixture(gaussian(), gaussian()), d,
   set_prior("normal(-2, 2)", class = "Intercept", dpar = "mu1"))
set.seed(3)
n <- 200; cl <- rbinom(n, 1, 0.4) + 1
pr <- rbind(c(0.85, 0.80, 0.75, 0.90), c(0.15, 0.20, 0.25, 0.10))
Y <- matrix(0L, n, 4); for (j in 1:4) Y[, j] <- 1L + rbinom(n, 1, pr[cl, j])
dd <- data.frame(x = rnorm(n)); dd$Y <- Y
go("lca", bf(Y ~ x), lca(K = 2), dd)
saveRDS(res, file.path(R2_ROOT, "dev/mvprior-log", paste0("sample-draws-", r2_arm, ".rds")))
