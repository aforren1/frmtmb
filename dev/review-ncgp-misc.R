# Items 3, 5: mutants in process (gp_fd_step and gp_field_dtheta), saved fits and the
# retape, readers the first pass spelled wrong, emtrends against the kernel's derivative.
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
patch <- function(nm, f) { unlockBinding(nm, ns); assign(nm, f, envir = ns); lockBinding(nm, ns) }
# ---- mutants: where does gp_fd_step() or the analytic d f / d theta change a number?
orig_step <- ns$gp_fd_step; orig_dth <- ns$gp_field_dtheta
fd_dth <- function(z, theta, blk) {
  out <- matrix(0, length(z), length(theta))
  bp <- function(th) as.vector(backsolve(ns$gp_factor_ordered(th, blk), matrix(0, length(z), 1), upper.tri = FALSE))
  L0 <- ns$gp_factor_ordered(theta, blk)
  for (k in seq_along(theta)) {
    h <- 1e-5
    tp <- theta; tp[k] <- tp[k] + h; tm <- theta; tm[k] <- tm[k] - h
    # d (sd L z) / d theta / sd, as the analytic column holds it in the factor order
    out[, k] <- as.vector((ns$gp_factor_ordered(tp, blk) - ns$gp_factor_ordered(tm, blk)) %*% z) / (2 * h)
  }
  # the analytic form's first column is z itself; solve L u = dL z for the same meaning
  u <- backsolve(L0, out, upper.tri = FALSE)
  u[, 1L] <- z
  u
}
cases <- list(
  low_n60 = list(n = 60, s = 1e-3), low_n200 = list(n = 200, s = 1e-3), vlow_n60 = list(n = 60, s = 1e-4), ord_n60 = list(n = 60, s = 0.2))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  set.seed(800 + cs$n)
  x <- sort(runif(cs$n, 0, 8)); z <- rnorm(cs$n)
  d <- data.frame(x = x, z = z, y = sin(2 * x) + cs$s * rnorm(cs$n), yl = exp(0.3 * sin(2 * x) + cs$s * rnorm(cs$n)))
  nd <- data.frame(x = c(0.5, 2.1, 4.33, 7.4), z = 0)
  mods <- list(gauss = list(y ~ gp(x)), lnorm = list(yl ~ gp(x), family = lognormal()), sig = list(bf(y ~ 1, sigma ~ gp(x))))
  for (mn in names(mods)) {
    m <- mods[[mn]]
    fit <- tryCatch(suppressWarnings(do.call(frm, c(m[1], list(data = d), m[-1]))), error = function(e) e)
    if (inherits(fit, "error")) { cat("RMU", nm, mn, "ERR\n"); next }
    rd <- function(f) {
      f$cache <- list2env(as.list(f$cache))
      if (is.environment(f$cache)) rm(list = intersect(ls(f$cache, all.names = TRUE), c("gp_jacobians", "Vjoint")), envir = f$cache)
      list(a = tryCatch(suppressWarnings(fitted(f, newdata = nd)), error = function(e) NA), b = tryCatch(suppressWarnings(fitted(f)), error = function(e) NA),
           c = tryCatch(suppressWarnings(unlist(ranef(f))), error = function(e) NA),
           d = tryCatch(suppressWarnings(as.data.frame(frm_curve_deriv(fit = f, var = "x", newdata = nd))$.se), error = function(e) NA))
    }
    r0 <- rd(fit)
    patch("gp_fd_step", function(g, fit) g); r1 <- rd(fit); patch("gp_fd_step", orig_step)
    patch("gp_field_dtheta", fd_dth); r2 <- rd(fit); patch("gp_field_dtheta", orig_dth)
    df <- function(a, b) { a <- suppressWarnings(as.numeric(unlist(a))); b <- suppressWarnings(as.numeric(unlist(b))); if (length(a) != length(b) || all(is.na(a))) return(NA); ok <- is.finite(a) & is.finite(b); max(abs(a[ok] - b[ok]) / pmax(abs(b[ok]), 1e-8)) }
    cat(sprintf("RMU %s %s conv %d | nofdstep: new %.1e fit %.1e ran %.1e deriv %.1e | fd_dtheta: new %.1e fit %.1e ran %.1e deriv %.1e\n", nm, mn, fit$opt$convergence,
                df(r1$a, r0$a), df(r1$b, r0$b), df(r1$c, r0$c), df(r1$d, r0$d), df(r2$a, r0$a), df(r2$b, r0$b), df(r2$c, r0$c), df(r2$d, r0$d)))
    flush.console()
  }
}
# ---- saved fits
set.seed(5); x <- sort(runif(60, 0, 6)); d <- data.frame(x = x, y = sin(x) + rnorm(60, sd = 0.2), z = rnorm(60), g = factor(rep(1:6, 10)))
fit <- suppressWarnings(frm(y ~ z + gp(x) + (1 | g), d))
tf <- tempfile(fileext = ".rds"); saveRDS(fit, tf); f2 <- readRDS(tf)
tcmp <- function(lab, a, b) cat(sprintf("RSV %-10s %s\n", lab, tryCatch(sprintf("%.1e", max(abs(as.numeric(unlist(a)) - as.numeric(unlist(b))), na.rm = TRUE)), error = function(e) conditionMessage(e))))
tcmp("fitted", suppressWarnings(fitted(fit)), tryCatch(suppressWarnings(fitted(f2)), error = function(e) conditionMessage(e)))
tcmp("vcov", as.matrix(vcov(fit)), tryCatch(as.matrix(vcov(f2)), error = function(e) conditionMessage(e)))
tcmp("confint", suppressWarnings(confint(fit)), tryCatch(suppressWarnings(confint(f2)), error = function(e) conditionMessage(e)))
tcmp("predict", suppressWarnings(fitted(fit, newdata = d[1:5, ])), tryCatch(suppressWarnings(fitted(f2, newdata = d[1:5, ])), error = function(e) conditionMessage(e)))
cat("compress after readers:", TMB::config(DLL = "RTMB")$tmbad.sparse_hessian_compress, " matmul ", RTMB::TapeConfig()$matmul, "\n")
# ---- readers spelled right this time
print(args(frm_curve))
rr <- function(lab, e) cat(sprintf("RRD %-30s %s\n", lab, tryCatch({ v <- suppressWarnings(suppressMessages(e)); paste(class(v)[1], paste(utils::head(signif(suppressWarnings(as.numeric(unlist(v))), 5), 6), collapse = " ")) }, error = function(err) paste("ERR", substr(conditionMessage(err), 1, 150)))))
fg <- suppressWarnings(frm(y ~ z + gp(x), d))
rr("frm_curve", frm_curve(fg, focal = "x", newdata = d[1:4, ]))
rr("conditional_effects", conditional_effects(fg, effects = "x"))
fg0 <- suppressWarnings(frm(y ~ z, d))
rr("anova 2 fits", anova(fg0, fg))
rr("anova ML/REML", anova(fg0, suppressWarnings(frm(y ~ gp(x), d))))
rr("diagnose", diagnose(fg))
rr("importance=", frm(y ~ gp(x), d, importance = 10))
rr("importance= (1|g) poisson", frm(cnt ~ gp(x) + (1 | g), transform(d, cnt = rpois(60, 3)), family = poisson(), importance = 10))
rr("quadrature", frm(y ~ gp(x), d, quadrature = TRUE))
rr("vcov_cluster", vcov_cluster(fg, cluster = d$g))
rr("cluster_scores", cluster_scores(fg, cluster = d$g))
rr("frm_simulate", frm_simulate(y ~ gp(x), data = d[1:5, ], par = list(beta = 0, betad = 0, theta = c(0, 0)), nsim = 1, seed = 1))
rr("ranef", ranef(fg))
rr("kriging far from data", fitted(fg, newdata = data.frame(x = 100, z = 0)))
if (requireNamespace("frmtmb.spatial", quietly = TRUE)) {
  library(frmtmb.spatial)
  cat("spatial lib:", dirname(system.file(package = "frmtmb.spatial")), "\n")
  dsp <- data.frame(s1 = runif(40), s2 = runif(40)); dsp$y <- sin(3 * dsp$s1) + rnorm(40, sd = 0.2); dsp$reg <- factor(rep(1:10, 4))
  rr("spatial gp alongside car", frm(y ~ gp(s1, s2) + (1 | reg), dsp))
}
# ---- emtrends against the kernel's derivative at the same x
set.seed(8150)
n <- 90; dd <- data.frame(x = runif(n, 0, 6)); dd$y <- sin(dd$x) + rnorm(n, sd = 0.3)
fe <- frm(y ~ gp(x), dd)
xb <- mean(dd$x)
em <- as.data.frame(emmeans::emtrends(fe, ~ 1, var = "x"))
de <- as.data.frame(frm_curve_deriv(fe, var = "x", newdata = data.frame(x = xb)))
cat(sprintf("REM emtrends %.5f se %.5f | frm_curve_deriv %.5f se %.5f at x %.4f\n", em$x.trend, em$SE, de$.estimate, de$.se, xb))
cat("DONE\n")
