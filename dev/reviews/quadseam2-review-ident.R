# Review of lane quadseam2, attacks 1 and 2: fits with no declared
# stationary point stalled, and quadrature fits not stalled, must be
# bitwise the same as 0.82.0. Prints logLik and the estimate vector in
# hex, SEs, warnings, times, and the escape records.
#   Rscript quadseam2-review-ident.R <libpath|base>
a <- commandArgs(TRUE)
if (a[1] == "lib") a[1] <- "C:/Users/adf44/source/r/wt-quadseam2-lib"
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.latent)})
cat("lib:", dirname(find.package("frmtmb")), "\n")
`%||%` <- function(a, b) if (is.null(a)) b else a
hx <- function(v) paste(sprintf("%a", v), collapse = ",")
run <- function(lbl, expr) {
  w <- character()
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(withCallingHandlers(expr,
    warning = function(x) {
      w <<- c(w, class(x)[1]); invokeRestart("muffleWarning")
    }, message = function(x) {
      w <<- c(w, paste0("m:", class(x)[1])); invokeRestart("muffleMessage")
    }), error = function(e) e)
  t1 <- proc.time()[["elapsed"]] - t0
  if (inherits(f, "error")) {
    cat(sprintf("FIT %s ERROR %s\n", lbl,
                substr(gsub("\\s+", " ", conditionMessage(f)), 1, 150)))
    return(invisible())
  }
  t0 <- proc.time()[["elapsed"]]
  se <- tryCatch(withCallingHandlers(sqrt(diag(vcov(f))),
    warning = function(x) {
      w <<- c(w, paste0("v:", class(x)[1])); invokeRestart("muffleWarning")
    }, message = function(x) {
      w <<- c(w, paste0("vm:", class(x)[1])); invokeRestart("muffleMessage")
    }), error = function(e) NA)
  t2 <- proc.time()[["elapsed"]] - t0
  cat(sprintf("FIT %s ll %s par %s\n", lbl, hx(as.numeric(logLik(f))),
              hx(f$opt$par)))
  cat(sprintf("SE  %s %s\n", lbl, hx(se)))
  cat(sprintf("AUX %s ll %.9f conv %d fit %.2f s vcov %.2f s esc %s sdesc %s warn %s\n",
              lbl, as.numeric(logLik(f)), f$opt$convergence, t1, t2,
              paste(f$opt$stationary_escape %||% "-", collapse = "/"),
              paste(f$opt$sd_escape %||% "-", collapse = "/"),
              paste(unique(w), collapse = ",")))
}
sn_draw <- function(n, mu, sigma, alpha) {
  dl <- alpha / sqrt(1 + alpha^2)
  mu + sigma * (dl * abs(rnorm(n)) + sqrt(1 - dl^2) * rnorm(n))
}
hmm_data <- function(seed, alpha = 5, gauss = FALSE) {
  set.seed(seed)
  G <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE)
  do.call(rbind, lapply(1:20, function(id) {
    s <- integer(40); s[1] <- sample.int(2, 1)
    for (t in 2:40) s[t] <- sample.int(2, 1, prob = G[s[t - 1], ])
    y <- if (gauss) ifelse(s == 1, rnorm(40, -2), rnorm(40, 2)) else
      ifelse(s == 1, sn_draw(40, -3, 1, alpha), sn_draw(40, 3, 1, -alpha))
    data.frame(id = id, t = 1:40, y = y)
  }))
}
set.seed(1)
n <- 200
d <- data.frame(x = rnorm(n), g = factor(rep(1:20, length.out = n)))
u <- rnorm(20, 0, 0.6)[d$g]
d$y <- 1 + 0.5 * d$x + u + rnorm(n)
d$yp <- rpois(n, exp(0.2 + 0.4 * d$x + u))
d$yb <- rbinom(n, 1, plogis(0.2 + 0.6 * d$x + u))
d$ys <- 1 + 0.5 * d$x + u + sn_draw(n, 0, 1, 4)
# sd just above the escape's region (0.01 on the link scale): true 0.05
u2 <- rnorm(20, 0, 0.05)[d$g]
d$yp2 <- rpois(n, exp(0.5 + 0.4 * d$x + u2))
d$yb0 <- rbinom(n, 1, plogis(0.2 + 0.6 * d$x))
set.seed(2)
dm <- data.frame(y = c(rnorm(150, -2), rnorm(150, 2)))
dsm <- data.frame(y = c(sn_draw(200, -3, 1, 5), sn_draw(200, 3, 1, -5)))
set.seed(3)
dn <- data.frame(ga = factor(rep(1:15, each = 12)),
                 gb = factor(rep(1:45, each = 4)), x = rnorm(180))
dn$y <- rbinom(180, 1, plogis(0.3 * dn$x + rnorm(15, 0, 0.8)[dn$ga] +
                                rnorm(45, 0, 0.5)[dn$gb]))
dn$tt <- factor(rep(1:4, 45))
run("gauss_lm", frm(y ~ x, data = d))
run("gauss_glmm", frm(y ~ x + (1 | g), data = d))
run("pois_glmm", frm(yp ~ x + (1 | g), family = poisson(), data = d))
run("bern_glmm", frm(yb ~ x + (1 | g), family = bernoulli(), data = d))
run("bern_quad", frm(yb ~ x + (1 | g), family = bernoulli(), data = d,
                     quadrature = TRUE))
run("pois_quad", frm(yp ~ x + (1 | g), family = poisson(), data = d,
                     quadrature = TRUE))
run("pois_quad_sd05", frm(yp2 ~ x + (1 | g), family = poisson(), data = d,
                          quadrature = TRUE))
run("bern_quad_sd0", frm(yb0 ~ x + (1 | g), family = bernoulli(), data = d,
                         quadrature = TRUE))
run("sn_glmm", frm(ys ~ x + (1 | g), family = skew_normal(), data = d))
run("sn_quad", frm(ys ~ x + (1 | g), family = skew_normal(), data = d,
                   quadrature = TRUE))
run("mix_gauss", frm(y ~ 1, family = mixture(gaussian(), gaussian()),
                     data = dm))
run("mix_sn", frm(y ~ 1, family = mixture(skew_normal(), skew_normal()),
                  data = dsm))
run("nested_quad", frm(bf(y ~ 1 + x + (1 | ga/gb)) + bernoulli(), data = dn,
                       quadrature = TRUE))
run("crossed_quad", frm(bf(y ~ 1 + x + (1 | ga) + (1 | tt)) + bernoulli(),
                        data = dn, quadrature = TRUE))
hd <- hmm_data(1, gauss = TRUE)
run("hmm_gauss", frm(bf(y ~ 1), family = hmm(K = 2, gaussian(), time = t,
                                             group = id), data = hd))
hs <- hmm_data(1)
run("hmm_sn", frm(bf(y ~ 1), family = hmm(K = 2, skew_normal(), time = t,
                                          group = id), data = hs))
# symmetric data: both alphas truly 0, the escape's legitimate case
h0 <- hmm_data(4, alpha = 0)
run("hmm_sn_sym", frm(bf(y ~ 1), family = hmm(K = 2, skew_normal(), time = t,
                                              group = id), data = h0))
set.seed(5)
ds0 <- data.frame(y = c(rnorm(200, -3), rnorm(200, 3)))
run("mix_sn_sym", frm(y ~ 1, family = mixture(skew_normal(), skew_normal()),
                      data = ds0))
cat("END\n")
