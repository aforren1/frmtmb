# L1 of the exact gp() on well-conditioned, UNSORTED designs, through
# the translator, on each spelling of gp() brms has.
#
#   Rscript dev/gpgap-identity.R <index>     # one spec; no index lists them
#
# Each spec is fitted once and read against brms's program three ways:
#   old    the translator before the fix (positions in frmtmb's order)
#   brms   brms's own program
#   frmnug brms's program with frmtmb's nugget in its one nugget line
# Seed 20261102. The translator is the worktree's helper-brms.R. Output is
# one GPID line per spec.
args <- commandArgs(trailingOnly = TRUE)
source("dev/gpgap-lib.R")

gpgap_data <- function() {
  set.seed(20261102)
  n <- 80
  # integer positions at a kernel distance the sine below keeps open
  d <- data.frame(x = sample(0:12, n, TRUE), z = sample(0:12, n, TRUE),
                  f = factor(sample(c("a", "b"), n, TRUE)),
                  w = runif(n, 0.5, 1.5))
  e <- rnorm(n, 0, 0.3)
  d$y <- 2 * sin(0.7 * d$x) + e
  d$y2 <- 2 * sin(0.7 * d$z) + rnorm(n, 0, 0.3)
  d$ybz <- 2 * sin(0.7 * d$x) + 1.5 * cos(0.7 * d$z) + e
  d$ybf <- ifelse(d$f == "a", 2 * sin(0.7 * d$x), 1.5 * cos(0.7 * d$x)) + e
  d$yw <- d$w * 2 * sin(0.7 * d$x) + e
  d$yl <- exp(0.5 * d$y / 2)
  d$yc <- rpois(n, exp(0.5 + sin(0.7 * d$x)))
  d$yb <- (2 + sin(0.7 * d$z)) * exp(0.3 * d$w) + e
  d$ys <- 2 * sin(0.7 * d$x) + rt(n, 5) * 0.3
  d$ysg <- rnorm(n, 0, exp(-1 + sin(0.7 * d$x)))
  d$xu <- round(d$x + runif(n, -0.04, 0.04), 3)
  d
}
sp <- list(
  list("gaussian gp(x)", "bf(y ~ gp(x))", "gaussian()"),
  list("lognormal gp(x)", "bf(yl ~ gp(x))", "lognormal()"),
  list("student gp(x)", "bf(ys ~ gp(x))", "student()"),
  list("poisson gp(x)", "bf(yc ~ gp(x))", "poisson()"),
  list("gp(x, z) iso", "bf(ybz ~ gp(x, z))", "gaussian()"),
  list("gp(x, z, iso = FALSE)", "bf(ybz ~ gp(x, z, iso = FALSE))",
       "gaussian()"),
  list("gp(x, by = f)", "bf(ybf ~ gp(x, by = f))", "gaussian()"),
  list("gp(x, by = f, cmc = FALSE)", "bf(ybf ~ gp(x, by = f, cmc = FALSE))",
       "gaussian()"),
  list("gp(x, by = w)", "bf(yw ~ gp(x, by = w))", "gaussian()"),
  list("gp(xu, gr = FALSE)", "bf(y ~ gp(xu, gr = FALSE))", "gaussian()"),
  list("gp(x, scale = FALSE)", "bf(y ~ gp(x, scale = FALSE))", "gaussian()"),
  list("sigma ~ gp(x)", "bf(ysg ~ 1, sigma ~ gp(x))", "gaussian()"),
  list("gp(x) + gp(z)", "bf(ybz ~ gp(x) + gp(z))", "gaussian()"),
  list("nl b1 ~ gp(z)",
       "bf(yb ~ b1 * exp(b2 * w), b1 ~ gp(z), b2 ~ 1, nl = TRUE)",
       "gaussian()"),
  list("mvbf gp(x), gp(z)",
       "mvbf(bf(y ~ gp(x)), bf(y2 ~ gp(z)), rescor = FALSE)", NA)
)
if (!length(args)) {
  for (i in seq_along(sp)) cat(i, sp[[i]][[1]], "\n")
  quit("no")
}
s <- sp[[as.integer(args[1])]]
E$s <- s
E$dd <- gpgap_data()
local(envir = E, {
  bform_f <- eq_eval(s[[2]], "frmtmb")
  a <- list(bform_f, data = dd)
  if (!is.na(s[[3]])) a$family <- eq_eval(s[[3]], "frmtmb")
  fit <- suppressWarnings(do.call(frmtmb::frm, a))
  bform <- eq_eval(s[[2]], "brms")
  fam <- if (is.na(s[[3]])) NULL else eq_eval(s[[3]], "brms")
  prior <- brms_flat_prior(bform, data = dd, family = fam)
  code <- brms::make_stancode(bform, data = dd, family = fam, prior = prior)
  sdat <- brms_standata(bform, data = dd, family = fam, prior = prior)
  rtab <- brms_ranef_table(bform, dd, fam, prior)
  sf0 <- gpgap_sf(code, sdat)
  sf1 <- gpgap_sf(gpgap_frm_nugget(code), sdat)
  ours0 <- -fit$obj$env$f(fit$obj$env$last.par.best)
  read <- function(p, sf) {
    rstan::log_prob(sf, rstan::unconstrain_pars(sf, p),
                    adjust_transform = FALSE) - (ours0 + attr(p, "logJ"))
  }
  pars <- stan_pars_from_fit(fit, sdat, code, rtab)
  new_order <- brms_gp_order
  # the mutant: frmtmb's order, which is the translator before the fix
  brms_gp_order <- function(sdat, zname) {
    seq_len(NROW(sdat[[sub("^zgp_", "Xgp_", zname)]]))
  }
  environment(stan_pars_from_fit) <- environment()
  p_old <- stan_pars_from_fit(fit, sdat, code, rtab)
  c_old <- read(p_old, sf0)
  c_brms <- read(pars, sf0)
  c_nug <- read(pars, sf1)
  u1 <- rstan::unconstrain_pars(sf1, pars)
  g1 <- rstan::grad_log_prob(sf1, u1, adjust_transform = FALSE)
  gmax <- max(abs(g1[brms_inner_index(sf1, pars)]))
  gp <- Filter(function(b) identical(b$covstruct, "gp"), fit$frame$re_blocks)
  conds <- vapply(gp, function(b) {
    v <- eigen(brms_block_cov(fit, b), symmetric = TRUE,
               only.values = TRUE)$values
    max(v) / min(v)
  }, 0)
  cat(sprintf(paste0("GPID %-27s old %12.6f  brms %12.5e  frmnug %12.5e  ",
                     "grad %.2e  ours %.3f  cond %s conv %d\n"),
              s[[1]], c_old, c_brms, c_nug, gmax, ours0 + attr(pars, "logJ"),
              paste(sprintf("%.2g", conds), collapse = ","),
              fit$opt$convergence))
})
