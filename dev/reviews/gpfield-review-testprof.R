# Review of lane gpfield: profile the log-likelihood over the length scale
# of the test fits that now allow frmtmb_gp_basis, which the lane did not
# profile (test-marginaleffects-re.R k = 15, test-gp-by.R k = 8 by f and
# by w, test-gp-multidim.R k = 4 iso). Each warned block: the fit's
# log-likelihood against the best with that block's length scale fixed at
# 0.05 to 4 spacings, every other parameter re-optimized.
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs("lane")
ns <- asNamespace("frmtmb")
prof <- function(tag, fit, warn) {
  ll <- as.numeric(logLik(fit))
  obj <- fit$obj; par <- fit$opt$par
  tpos_all <- which(names(par) == "theta")
  for (bk in fit$frame$re_blocks) {
    if (bk$covstruct != "hsgp") next
    th <- fit$estimates$theta[bk$theta_idx]
    white <- ns$gp_basis_white(bk, th)
    sp <- min(bk$gp_spacing)
    tl <- tpos_all[bk$theta_idx[-1]]
    pv <- c()
    for (m in c(0.05, 0.25, 0.5, 1, 2, 4)) {
      p0 <- par; p0[tl] <- log(m * sp)
      free <- setdiff(seq_along(par), tl)
      fn <- function(v) { p <- p0; p[free] <- v; obj$fn(p) }
      gr <- function(v) { p <- p0; p[free] <- v; obj$gr(p)[free] }
      best <- Inf
      for (st in 1:2) {
        v0 <- par[free]
        if (st == 2) {
          sp_i <- which(free == tpos_all[bk$theta_idx[1]])
          v0[sp_i] <- th[1] + 0.5 * length(tl) * (mean(th[-1]) - log(m * sp))
        }
        o <- tryCatch(suppressWarnings(stats::nlminb(v0, fn, gr)),
                      error = function(e) NULL)
        if (!is.null(o) && is.finite(o$objective)) best <- min(best, o$objective)
      }
      pv <- c(pv, -best - ll)
    }
    cat(sprintf("TPROF %-22s %-28s white %d ls/spacing %s | prof - ll at .05 .25 .5 1 2 4 spacings: %s | warn [%s]\n",
                tag, bk$group_name %||% bk$label, white,
                paste(signif(exp(th[-1]) / bk$gp_spacing, 3), collapse = ","),
                paste(sprintf("%.3f", pv), collapse = " "),
                paste(unique(sub(":.*", "", warn)), collapse = ";")))
  }
}
set.seed(7411)
n <- 200
d <- data.frame(x = stats::rnorm(n))
d$y <- sin(2 * d$x) + stats::rnorm(n, 0, 0.4)
q <- gf_fit(y ~ gp(x, k = 15), data = d)
prof("marginaleffects k15", q$fit, q$warn)
set.seed(5)
n <- 60
d <- data.frame(x = round(stats::runif(n, 0, 6), 1),
                f = factor(rep(c("a", "b", "c"), length.out = n)),
                w = stats::runif(n, 0.5, 2))
d$y <- 0.5 + ifelse(d$f == "a", sin(d$x),
                    ifelse(d$f == "b", cos(d$x), 0.2 * d$x)) +
  stats::rnorm(n, 0, 0.3)
for (form in c("y ~ gp(x, by = f, k = 8)", "y ~ gp(x, by = w, k = 8)")) {
  q <- gf_fit(bf(stats::as.formula(form)), data = d)
  prof(form, q$fit, q$warn)
}
