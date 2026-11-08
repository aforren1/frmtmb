# Convergence of the nested rule as the nodes rise, against references
# written from the model alone; the inner sd near zero against the
# one-block Gauss-Kronrod route; gaussian against its closed form.
#   Rscript dev/nestquad-nodes.R <case: bern|pois|gaus|zero>
# Seeds in each case. Every value is the nested objective re-anchored
# at the stated parameters.
a <- commandArgs(TRUE)
wt <- "C:/Users/adf44/source/r/frmtmb-wt-nestquad"
.libPaths(c("C:/Users/adf44/source/r/wt-nestquad-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
source(file.path(wt, "dev/nestquad-ref.R"))

# the nested log-likelihood at outer parameters `par`, rule anchored there
nest_at <- function(fit, par, nodes) {
  fr <- fit$frame
  nll <- ns$build_objective(fr)
  tpl <- fr$par_template
  lap <- RTMB::MakeADFun(nll, tpl, random = "b", map = fr$map,
                         silent = TRUE)
  plan <- ns$quad_check_blocks(fr)
  anc <- ns$quad_nest_anchor(lap, plan, par)
  pl <- ns$imp_par_list(tpl, par, fr$map)
  qo <- ns$quad_nest_objective(fr, plan, anc, nodes, pl)
  pl <- ns$imp_par_list(tpl, par, fr$map)
  qo$prime(pl)
  -qo$fn(pl)
}
lap_at <- function(fit, par) {
  fr <- fit$frame
  lap <- RTMB::MakeADFun(ns$build_objective(fr), fr$par_template,
                         random = "b", map = fr$map, silent = TRUE)
  -lap$fn(par)
}
sd_of <- function(fit, par) {
  th <- ns$imp_par_list(fit$frame$par_template, par, fit$frame$map)$theta
  s <- vapply(fit$frame$re_blocks, function(bk) {
    sqrt(ns$cs_entry(bk$covstruct)$vcov(th[bk$theta_idx], bk)[1, 1])
  }, 1)
  names(s) <- vapply(fit$frame$re_blocks, `[[`, "", "group_name")
  s
}
ks <- c(1, 2, 3, 5, 7, 9, 11, 15, 21, 31)
report <- function(tag, fit, par, ref, extra = "") {
  la <- lap_at(fit, par)
  cat(sprintf("REF %s ref %.10f laplace %.10f err_laplace %.3e %s\n", tag,
              ref, la, abs(la - ref), extra))
  for (k in ks) {
    v <- nest_at(fit, par, c(k, k))
    cat(sprintf("NODES %s k %2d loglik %.10f err %.3e ratio_to_laplace %.3e\n",
                tag, k, v, abs(v - ref), abs(v - ref) / abs(la - ref)))
  }
}

case <- a[1]
if (case == "bern") {
  set.seed(20261020)
  ng <- 12; nt <- 6; n <- ng * nt
  d <- data.frame(ga = factor(rep(seq_len(ng), each = nt)),
                  gb = factor(rep(c("i", "ii"), length.out = n)),
                  x = rnorm(n), w = runif(n, 0.5, 2))
  ab <- interaction(d$ga, d$gb, drop = TRUE)
  d$eta <- 0.2 + 0.5 * d$x + rnorm(ng, 0, 0.7)[d$ga] +
    rnorm(nlevels(ab), 0, 1)[ab]
  d$y <- rbinom(n, 1, plogis(d$eta))
  fit <- frm(bf(y | weights(w) ~ 1 + x + (1 | ga/gb)), family = bernoulli(),
             data = d, quadrature = TRUE)
  par <- fit$opt$par
  s <- sd_of(fit, par)
  b <- unname(fit$estimates$beta)
  e0 <- b[1] + b[2] * d$x
  r80 <- ref_nested_gh("bernoulli", d$y, e0, d$ga, ab, s[["ga"]],
                       s[["ga:gb"]], w = d$w, k = 80)
  r120 <- ref_nested_gh("bernoulli", d$y, e0, d$ga, ab, s[["ga"]],
                        s[["ga:gb"]], w = d$w, k = 120)
  cat(sprintf("REFCONV bern gh80 %.12f gh120 %.12f diff %.3e\n", r80, r120,
              abs(r80 - r120)))
  report("bern", fit, par, r120,
         sprintf("fit_loglik %.10f", as.numeric(logLik(fit))))
} else if (case == "pois") {
  set.seed(20261107)
  no <- 8; npo <- 3; rpi <- 8; n <- no * npo * rpi
  d <- data.frame(ga = factor(rep(seq_len(no), each = npo * rpi)),
                  gb = factor(rep(rep(seq_len(npo), each = rpi), no)),
                  x = rnorm(n))
  ab <- interaction(d$ga, d$gb, drop = TRUE)
  d$y <- rpois(n, exp(0.5 + 0.3 * d$x + rnorm(no, 0, 0.8)[d$ga] +
                        rnorm(nlevels(ab), 0, 0.5)[ab]))
  fit <- frm(bf(y ~ 1 + x + (1 | ga/gb)), family = poisson(), data = d,
             quadrature = TRUE)
  par <- fit$opt$par
  s <- sd_of(fit, par)
  b <- unname(fit$estimates$beta)
  e0 <- b[1] + b[2] * d$x
  r1 <- ref_nested_integrate("poisson", d$y, e0, d$ga, ab, s[["ga"]],
                             s[["ga:gb"]], rel.tol = 1e-9)
  r2 <- ref_nested_integrate("poisson", d$y, e0, d$ga, ab, s[["ga"]],
                             s[["ga:gb"]], rel.tol = 1e-11)
  cat(sprintf("REFCONV pois int1e-9 %.12f int1e-11 %.12f diff %.3e\n", r1,
              r2, abs(r1 - r2)))
  report("pois", fit, par, r2,
         sprintf("fit_loglik %.10f", as.numeric(logLik(fit))))
} else if (case == "gaus") {
  set.seed(20261108)
  no <- 10; npo <- 3; rpi <- 4; n <- no * npo * rpi
  d <- data.frame(ga = factor(rep(seq_len(no), each = npo * rpi)),
                  gb = factor(rep(rep(seq_len(npo), each = rpi), no)),
                  x = rnorm(n))
  ab <- interaction(d$ga, d$gb, drop = TRUE)
  d$y <- 1 + 0.4 * d$x + rnorm(no, 0, 0.9)[d$ga] +
    rnorm(nlevels(ab), 0, 0.6)[ab] + rnorm(n, 0, 0.5)
  fit <- frm(bf(y ~ 1 + x + (1 | ga/gb)), family = gaussian(), data = d,
             quadrature = TRUE)
  par <- fit$opt$par
  s <- sd_of(fit, par)
  b <- unname(fit$estimates$beta)
  sig <- exp(ns$imp_par_list(fit$frame$par_template, par,
                             fit$frame$map)$betad[1])
  Za <- outer(as.integer(d$ga), seq_len(no), "==") * 1
  Zb <- outer(as.integer(ab), seq_len(nlevels(ab)), "==") * 1
  Vm <- s[["ga"]]^2 * tcrossprod(Za) + s[["ga:gb"]]^2 * tcrossprod(Zb) +
    diag(sig^2, n)
  ch <- chol(Vm)
  rr <- backsolve(ch, d$y - b[1] - b[2] * d$x, transpose = TRUE)
  ref <- -0.5 * sum(rr^2) - sum(log(diag(ch))) - n / 2 * log(2 * pi)
  report("gaus", fit, par, ref,
         sprintf("sigma %.6f fit_loglik %.10f", sig,
                 as.numeric(logLik(fit))))
} else if (case == "zero") {
  set.seed(20261020)
  ng <- 12; nt <- 6; n <- ng * nt
  d <- data.frame(ga = factor(rep(seq_len(ng), each = nt)),
                  gb = factor(rep(c("i", "ii"), length.out = n)),
                  x = rnorm(n))
  ab <- interaction(d$ga, d$gb, drop = TRUE)
  d$y <- rbinom(n, 1, plogis(0.2 + 0.5 * d$x + rnorm(ng, 0, 0.9)[d$ga]))
  fs <- frm(bf(y ~ 1 + x + (1 | ga)), family = bernoulli(), data = d,
            quadrature = TRUE)
  gk <- as.numeric(logLik(fs))
  b <- unname(fs$estimates$beta)
  s1 <- sd_of(fs, fs$opt$par)[["ga"]]
  ref1 <- 0
  for (aa in levels(d$ga)) {
    r <- d$ga == aa
    ref1 <- ref1 + ref_log_int(function(u) {
      stats::dnorm(u, 0, s1, log = TRUE) +
        sum(ref_lpdf("bernoulli", d$y[r], b[1] + b[2] * d$x[r] + u))
    })
  }
  fn <- frm(bf(y ~ 1 + x + (1 | ga/gb)), family = bernoulli(), data = d,
            quadrature = TRUE)
  th_a <- fs$opt$par[names(fs$opt$par) == "theta"]
  for (sd_in in c(1e-1, 1e-2, 1e-3, 1e-4)) {
    p <- fn$opt$par
    p[names(p) == "beta"] <- fs$opt$par[names(fs$opt$par) == "beta"]
    th <- which(names(p) == "theta")
    # theta is log sd for a dimension-1 block (checked below)
    p[th] <- c(th_a, log(sd_in))
    s <- sd_of(fn, p)
    v <- nest_at(fn, p, c(11, 11))
    cat(sprintf(paste("ZERO sd_inner %.0e sds %.6f %.6f nested %.10f",
                      "gk_one_block %.10f ref_one_block %.10f",
                      "nested-gk %.3e gk-ref %.3e\n"),
                sd_in, s[[1]], s[[2]], v, gk, ref1, v - gk, gk - ref1))
  }
}
if (case == "bern" && length(a) > 1) {
  for (kk in list(c(31, 11), c(11, 31), c(31, 5), c(5, 31), c(11, 15),
                  c(11, 21), c(15, 21), c(9, 15), c(9, 21), c(7, 21))) {
    v <- nest_at(fit, par, kk)
    cat(sprintf("ASYM bern ko %2d ki %2d loglik %.10f err %.3e\n", kk[1],
                kk[2], v, abs(v - r120)))
  }
}
