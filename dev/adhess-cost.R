# Lane adhess: cost of the covariance, before (base) and after (lane).
# Per model: the fit (frm()) and the first vcov(), which builds the
# sdreport, each the minimum over `reps` repeats in this process, plus
# one profiled repeat (Rprof, 10 ms) whose share in the Hessian
# builders is printed. Seeds fixed per model.
#   ADHESS_ARM=lane|base Rscript dev/adhess-cost.R <model> [reps]
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-common.R")
a <- commandArgs(TRUE)
model <- a[1]
reps <- if (length(a) > 1) as.integer(a[2]) else 3L
mk <- list(
  glm_pois_p3_n1e3 = function() {
    set.seed(1); n <- 1e3
    d <- data.frame(x1 = rnorm(n), x2 = rnorm(n, 500, 100))
    d$y <- rpois(n, exp(0.3 + 0.3 * d$x1 + 0.002 * (d$x2 - 500)))
    list(f = bf(y ~ x1 + x2), fam = poisson(), d = d)
  },
  glm_pois_p10_n1e5 = function() {
    set.seed(2); n <- 1e5
    X <- matrix(rnorm(n * 9), n)
    d <- as.data.frame(X)
    d$y <- rpois(n, exp(0.3 + drop(X %*% rep(0.1, 9))))
    list(f = bf(as.formula(paste("y ~", paste(names(d)[1:9],
                                              collapse = " + ")))),
         fam = poisson(), d = d)
  },
  glm_bern_p50_n1e5 = function() {
    set.seed(3); n <- 1e5
    X <- matrix(rnorm(n * 49), n)
    d <- as.data.frame(X)
    d$y <- rbinom(n, 1, plogis(0.2 + drop(X %*% rep(0.03, 49))))
    list(f = bf(as.formula(paste("y ~", paste(names(d)[1:49],
                                              collapse = " + ")))),
         fam = bernoulli(), d = d)
  },
  lm_p10_n1e5 = function() {
    set.seed(4); n <- 1e5
    X <- matrix(rnorm(n * 9), n)
    d <- as.data.frame(X)
    d$y <- 1 + drop(X %*% rep(0.1, 9)) + rnorm(n)
    list(f = bf(as.formula(paste("y ~", paste(names(d)[1:9],
                                              collapse = " + ")))),
         fam = gaussian(), d = d)
  },
  glmm_pois_n1e3 = function() {
    set.seed(5); n <- 1e3
    d <- data.frame(x = rnorm(n, 500, 100), g = factor(rep(1:50, length.out = n)))
    d$y <- rpois(n, exp(0.3 + 0.004 * (d$x - 500) + rnorm(50, 0, 0.3)[d$g]))
    list(f = bf(y ~ x + (1 | g)), fam = poisson(), d = d)
  },
  glmm_pois_n1e5 = function() {
    set.seed(6); n <- 1e5
    d <- data.frame(x1 = rnorm(n), x2 = rnorm(n, 500, 100),
                    g = factor(rep(1:500, length.out = n)))
    d$y <- rpois(n, exp(0.3 + 0.2 * d$x1 + 0.002 * (d$x2 - 500) +
                          rnorm(500, 0, 0.3)[d$g]))
    list(f = bf(y ~ x1 + x2 + (1 | g)), fam = poisson(), d = d)
  }
)
cs <- mk[[model]]()
one <- function() {
  t0 <- proc.time()[["elapsed"]]
  f <- suppressWarnings(frm(cs$f, family = cs$fam, data = cs$d))
  t1 <- proc.time()[["elapsed"]]
  V <- suppressWarnings(vcov(f, full = TRUE))
  t2 <- proc.time()[["elapsed"]]
  list(fit = t1 - t0, vcov = t2 - t1, f = f, V = V)
}
tt <- t(vapply(seq_len(reps), function(i) {
  r <- one()
  c(r$fit, r$vcov)
}, numeric(2)))
pf <- file.path(tempdir(), "adhess-prof.out")
Rprof(pf, interval = 0.01)
r <- one()
Rprof(NULL)
sp <- summaryRprof(pf)$by.total
keep <- grep("fd_hessian|obj_he_cached|outer_hessian_ad|fit_outer_hessian|sdreport|optimHess|se_fit_analysis|outer_fd_steps|\"frm\"",
             rownames(sp))
cat(sprintf("%s %s | fit min %.2f s (%s) | vcov min %.2f s (%s) | np %d\n",
            model, arm, min(tt[, 1]), paste(sprintf("%.2f", tt[, 1]),
                                             collapse = ","),
            min(tt[, 2]), paste(sprintf("%.2f", tt[, 2]), collapse = ","),
            length(r$f$opt$par)))
print(sp[keep, c("total.time", "total.pct")])
cat("SE:", format(sqrt(diag(r$V))[seq_len(min(4, NROW(r$V)))],
                  digits = 8), "\n")
