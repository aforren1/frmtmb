# Item 6: every consumer of an importance fit, against a reference where one
# exists. Prints each result, or the error it raised.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs()
suppressPackageStartupMessages(library(frmtmb))
try_print <- function(lab, expr) {
  w <- character(0)
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "ERR")),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  cat("\n---", lab, "---\n")
  if (inherits(v, "ERR")) cat("ERROR:", substr(unclass(v), 1, 400), "\n") else print(v)
  if (length(w)) cat("WARN:", substr(paste(w, collapse = " | "), 1, 400), "\n")
  invisible(v)
}
set.seed(21)
G <- 20
g <- factor(rep(seq_len(G), each = 6))
x <- rnorm(120)
z <- rnorm(120)
u <- rnorm(G, 0, 0.9)
y <- rbinom(120, 1, plogis(0.2 + 0.6 * x + u[g]))
d <- data.frame(y, x, z, g)
f <- bf(y ~ x + (1 | g))
fi <- frm(f, family = bernoulli(), data = d, importance = 500L)
fq <- frm(f, family = bernoulli(), data = d, quadrature = TRUE)
fl <- frm(f, family = bernoulli(), data = d)

try_print("logLik imp / quad / lap", c(logLik(fi), logLik(fq), logLik(fl)))
try_print("mcse", fi$importance$mcse)
try_print("vcov imp / quad ratio", diag(vcov(fi)) / diag(vcov(fq)))
try_print("se theta", {
  s <- frmtmb:::sdr_of(fi)
  sq <- frmtmb:::sdr_of(fq)
  sqrt(diag(s$cov.fixed)) / sqrt(diag(sq$cov.fixed))
})
try_print("confint wald", confint(fi, method = "wald"))
ci_i <- try_print("confint profile theta imp", confint(fi, parm = "theta_1", method = "profile"))
ci_q <- try_print("confint profile theta quad", confint(fq, parm = "theta_1", method = "profile"))
try_print("confint profile x imp", confint(fi, parm = "x", method = "profile"))
try_print("confint profile x quad", confint(fq, parm = "x", method = "profile"))
try_print("predict response, newdata", head(predict(fi, newdata = d[1:3, ], type = "response")))
try_print("predict re_formula NA", head(predict(fi, re_formula = NA)))
try_print("fitted", head(fitted(fi)))
try_print("ranef", head(ranef(fi)$g))
try_print("simulate", dim(simulate(fi, nsim = 2)))
try_print("summary", summary(fi))
try_print("diagnose", diagnose(fi, quiet = TRUE)[c("max_grad")])
fi2 <- try_print("update drop x", update(fi, . ~ . - x))
try_print("update keeps importance", fi2$importance[c("draws", "seed", "rounds", "ess_min")])
try_print("update keeps seed", {
  a <- update(fi, control = frmtmb_control(importance_seed = 7L))
  c(a$importance$seed, logLik(a))
})
try_print("anova imp, imp", anova(fi2, fi))
try_print("anova imp, lap", anova(fi2, fl))
try_print("drop1", drop1(fi))
try_print("refit on own data, no data arg", {
  dd <- d
  fx <- frm(f, family = bernoulli(), data = dd, importance = 500L)
  rm(dd)
  u <- update(fx, importance = 200L)
  c(u$importance$draws, as.numeric(logLik(u)))
})
try_print("same seed same answer", {
  a <- frm(f, family = bernoulli(), data = d, importance = 500L)
  identical(a$opt$par, fi$opt$par)
})
try_print("session RNG undisturbed", {
  set.seed(1); r1 <- runif(1)
  set.seed(1); invisible(frm(f, family = bernoulli(), data = d, importance = 16L)); r2 <- runif(1)
  r1 == r2
})
try_print("frm_bootstrap", {
  if (exists("frm_bootstrap")) {
    b <- frm_bootstrap(fi, nsim = 2, seed = 1)
    class(b)
  }
})
try_print("frm_allfit", {
  if (exists("frm_allfit")) {
    a <- frm_allfit(fi)
    class(a)
  }
})
try_print("influence", {
  if (exists("frm_influence")) dim(frm_influence(fi))
})
try_print("cooks", head(stats::cooks.distance(fi)))
try_print("vcov_cluster refuses", vcov_cluster(fi))
try_print("marginal curve", {
  if (exists("frm_curve")) head(frm_curve(fi, "x"))
})

# saved fit, in a child process with nothing else loaded
tf <- tempfile(fileext = ".rds")
saveRDS(fi, tf)
ch <- tempfile(fileext = ".R")
writeLines(c(
  sprintf('source("%s")', "C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R"),
  'set_libs()',
  'suppressPackageStartupMessages(library(frmtmb))',
  sprintf('fi <- readRDS("%s")', gsub("\\\\", "/", tf)),
  'cat("logLik", logLik(fi), "\\n")',
  'print(try(summary(fi)$fixed))',
  'print(try(head(predict(fi, newdata = fi$data[1:2, ], type = "response"))))',
  'print(try(confint(fi, method = "wald")))',
  'print(try(confint(fi, parm = "theta_1", method = "profile")))',
  'print(try(vcov(fi)))',
  'print(try(update(fi, . ~ . - x)$importance$draws))',
  'print(try(frmtmb:::imp_ess_at(frmtmb:::imp_frozen_proposal(fi), fi$opt$par)$mcse))',
  'print(fi$importance$mcse)'), ch)
cat("\n--- saved fit, child process ---\n")
cat(system2(file.path(R.home("bin"), "Rscript.exe"), shQuote(ch), stdout = TRUE,
            stderr = TRUE), sep = "\n")

# prior: is it in the objective at all?
try_print("prior on b: objective vs no prior", {
  fp <- frm(f, family = bernoulli(), data = d, importance = 500L,
            prior = set_prior("normal(0, 0.1)", class = "b"))
  fp2 <- frm(f, family = bernoulli(), data = d,
             prior = set_prior("normal(0, 0.1)", class = "b"))
  rbind(imp_prior = fp$opt$par, lap_prior = fp2$opt$par, imp = fi$opt$par)
})
try_print("prior on sd", {
  fp <- frm(f, family = bernoulli(), data = d, importance = 500L,
            prior = set_prior("normal(0, 0.1)", class = "sd"))
  fp2 <- frm(f, family = bernoulli(), data = d,
             prior = set_prior("normal(0, 0.1)", class = "sd"))
  rbind(imp_prior = fp$opt$par, lap_prior = fp2$opt$par, imp = fi$opt$par)
})
try_print("bounds on b", {
  fb <- frm(f, family = bernoulli(), data = d, importance = 200L,
            prior = prior(normal(0, 1000), class = "b", lb = 0.8))
  fb$opt$par
})
# autoscale
try_print("autoscale", {
  dd <- d
  dd$xb <- dd$x * 1e6
  a <- frm(bf(y ~ xb + (1 | g)), family = bernoulli(), data = dd,
           importance = 200L, control = frmtmb_control(autoscale = TRUE))
  b <- frm(bf(y ~ x + (1 | g)), family = bernoulli(), data = dd,
           importance = 200L)
  c(logLik(a), logLik(b), a$opt$par[2] * 1e6, b$opt$par[2], a$importance$rounds)
})
# frmtmb.learn and frmtmb.spatial
try_print("learn loads", requireNamespace("frmtmb.learn", quietly = TRUE))
try_print("spatial loads", requireNamespace("frmtmb.spatial", quietly = TRUE))
try_print("spatial car refuses", {
  library(frmtmb.spatial)
  A <- matrix(0, 8, 8)
  A[cbind(1:7, 2:8)] <- 1
  A <- A + t(A)
  dimnames(A) <- list(as.character(1:8), as.character(1:8))
  ds <- data.frame(y = rpois(80, 2), x = rnorm(80),
                   loc = factor(rep(as.character(1:8), 10)))
  frm(bf(y ~ x + car(1 | loc, M = A)), family = poisson(), data = ds,
      data2 = list(A = A), importance = 100L)
})
try_print("spatial sar refuses", {
  frm(bf(y ~ x + (1 | loc) + sar(M)), family = gaussian(),
      data = data.frame(y = rnorm(80), x = rnorm(80),
                        loc = factor(rep(1:8, 10))),
      data2 = list(M = diag(80)), importance = 100L)
})
