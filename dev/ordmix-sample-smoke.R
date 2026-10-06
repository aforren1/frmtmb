# frmtmb.sample on ordinal mixtures: draw names, posterior_epred against
# the theta-weighted sum of brms's R-side category probabilities at the
# stored columns, log_lik, posterior_predict. Data seed 20261005,
# sampler seed 3. Usage: Rscript dev/ordmix-sample-smoke.R [lane|base]
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "lane"
libs <- c("C:/Users/adf44/source/r/rellib-r5",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ordmix-lib", libs)
.libPaths(libs)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-ordmix/dev/stan-cache")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb:", find.package("frmtmb"), " sample:",
    find.package("frmtmb.sample"), "\n")
set.seed(20261005)
n <- 300
x <- rnorm(n)
z <- rnorm(n)
g <- factor(sample(c("a", "b"), n, TRUE))
cls <- rbinom(n, 1, 0.4)
lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
y <- as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf)))
d <- data.frame(y, x, z, g)
tryf <- function(label, expr) {
  cat("\n=====", label, "\n")
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  if (inherits(r, "error")) cat("ERROR:", conditionMessage(r), "\n") else
    print(r)
  invisible(r)
}
# brms's R-side category probabilities of one component
comp_p <- function(fam, eta, thres, disc = 1) {
  dens <- get(paste0("d", fam), asNamespace("brms"))
  n <- length(eta)
  dens(seq_len(length(thres) + 1L), eta = eta,
       thres = matrix(thres, n, length(thres), byrow = TRUE), disc = disc,
       link = "logit")
}
# brms's default student_t on the thresholds, written by dpar as brms
# keys the rows, and a weak prior on the slopes: flat thresholds leave
# the sampler an improper direction
pr_none <- set_prior("student_t(3, 0, 2.5)", class = "Intercept",
                     dpar = "mu1") +
  set_prior("student_t(3, 0, 2.5)", class = "Intercept", dpar = "mu2") +
  set_prior("normal(0, 2)", class = "b", dpar = "mu1") +
  set_prior("normal(0, 2)", class = "b", dpar = "mu2")
pr_mu <- set_prior("student_t(3, 0, 2.5)", class = "Intercept") +
  set_prior("normal(0, 2)", class = "b", dpar = "mu1") +
  set_prior("normal(0, 2)", class = "b", dpar = "mu2")
check <- function(fit, fams, label) {
  cat("\n############", label, "\n")
  ds <- tryf("frm_sample", frm_sample(fit, chains = 1, iter = 400,
                                      warmup = 200, seed = 3,
                                      refresh = 0))
  if (inherits(ds, "error")) return(invisible())
  tryf("variables", variables(ds))
  ep <- tryf("posterior_epred dim", {
    e <- posterior_epred(ds)
    dim(e)
  })
  M <- posterior::as_draws_matrix(ds)
  e <- posterior_epred(ds)
  # the theta-weighted sum at the first and last draw, from the stored
  # columns: thresholds b_mu<k>_Intercept[j], slopes b_mu<k>_x, theta
  for (s in c(1L, nrow(M))) {
    P <- 0
    for (k in seq_along(fams)) {
      th <- as.numeric(M[s, grep(paste0("^b_mu", k, "_Intercept\\["),
                                 colnames(M))])
      eta <- as.numeric(M[s, paste0("b_mu", k, "_x")]) * d$x
      P <- P + as.numeric(M[s, paste0("theta", k)]) *
        comp_p(fams[k], eta, th)
    }
    cat(sprintf("draw %d: max |epred - theta-weighted brms sum| = %.3g\n",
                s, max(abs(e[s, , ] - P))))
  }
  tryf("log_lik dim", dim(log_lik(ds)))
  tryf("posterior_predict table", table(posterior_predict(ds)[1, ]))
  tryf("fixef", fixef(ds))
  tryf("summary", summary(ds))
}
f1 <- frm(bf(y ~ x), family = mixture(cumulative(), sratio()), data = d,
          prior = pr_none)
check(f1, c("cumulative", "sratio"), "order none, cumulative + sratio")
f2 <- frm(bf(y ~ x), family = mixture(cumulative(), acat(), order = "mu"),
          data = d, prior = pr_mu)
check(f2, c("cumulative", "acat"), "order mu, cumulative + acat")
f3 <- frm(bf(y ~ x), family = mixture(cumulative(),
                                      cratio(threshold = "sum_to_zero"),
                                      order = "mu"), data = d,
          prior = pr_mu)
check(f3, c("cumulative", "cratio"), "order mu, flexible + sum_to_zero")
