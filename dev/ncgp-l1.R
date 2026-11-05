# L1 against brms's own program, unedited, on the spellings of an exact gp().
# Through the translator stan_pars_from_fit(), comparing log_prob with frm()'s
# joint density at the estimate and at two stepped points, and the inner
# gradient at the estimate, as dev/equiv-findings.md 2.5 states them. Under
# "attr" the same fits are read on brms's program with its nugget line
# rewritten to frmtmb's relative one: ATTRIBUTION of the constant, never a
# verdict. Seed 20261102. One NCGPL line per spelling.
#   Rscript dev/ncgp-l1.R <index> [attr]
args <- commandArgs(TRUE)
LIB <- "C:/Users/adf44/source/r/wt-ncgp-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
})
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("helper-warnings.R", "helper-brms.R")) {
  sys.source(file.path("tests/testthat", f), envir = E)
}
set.seed(20261102)
n <- 60L
x <- round(stats::runif(n, 0, 6), 3)
z <- round(stats::runif(n, 0, 4), 2)
f <- factor(rep(c("a", "b", "c"), length.out = n))
w <- stats::runif(n, 0.5, 2)
g <- factor(rep(seq_len(6L), length.out = n))
e <- stats::rnorm(n)
y <- sin(2 * x) + 0.3 * e
d <- data.frame(x = x, z = z, f = f, w = w, g = g, e = e, y = y,
                y2 = cos(z) + 0.3 * stats::rnorm(n),
                yl = exp(0.5 * sin(2 * x) + 0.3 * e),
                cnt = stats::rpois(n, exp(1 + sin(2 * x))),
                bin = stats::rbinom(n, 1, plogis(sin(2 * x))),
                xr = rep(round(stats::runif(15L, 0, 6), 3), length.out = n),
                se = 0.3)
d$yr <- sin(2 * d$xr) + 0.3 * d$e
d$ynl <- (2 + sin(z)) * exp(0.3 * w) + 0.2 * e
d$ysg <- stats::rnorm(n, 0, exp(-1 + sin(x)))
d$ymx <- ifelse(as.integer(g) %% 2L == 0L, sin(2 * x), 3) + 0.3 * e
d$yt <- sin(2 * x) + 0.3 * stats::rt(n, 5)
spell <- list(
  list("gp(x)", "bf(y ~ gp(x))"),
  list("gp(x, scale = FALSE)", "bf(y ~ gp(x, scale = FALSE))"),
  list("gp(x, z)", "bf(y ~ gp(x, z))"),
  list("gp(x, z, iso = FALSE)", "bf(y ~ gp(x, z, iso = FALSE))"),
  list("gp(x, by = f)", "bf(y ~ gp(x, by = f))"),
  list("gp(x, by = f, cmc = FALSE)", "bf(y ~ gp(x, by = f, cmc = FALSE))"),
  list("gp(x, by = w)", "bf(y ~ gp(x, by = w))"),
  list("gp(xr, gr = FALSE)", "bf(yr ~ gp(xr, gr = FALSE))"),
  list("gp(xr) replicated", "bf(yr ~ gp(xr))"),
  list("gp(x) + gp(z)", "bf(y ~ gp(x) + gp(z))"),
  list("sigma ~ gp(x)", "bf(ysg ~ 1, sigma ~ gp(x))"),
  list("nl b1 ~ gp(z)", "bf(ynl ~ b1 * exp(b2 * w), b1 ~ gp(z), b2 ~ 1, nl = TRUE)"),
  list("mvbf gp(x), gp(z)", "mvbf(bf(y ~ gp(x)), bf(y2 ~ gp(z)), rescor = FALSE)"),
  list("gp(x) + (1 | g)", "bf(y ~ gp(x) + (1 | g))"),
  list("se() gp(x)", "bf(y | se(se) ~ gp(x))"),
  list("student gp(x)", "bf(yt ~ gp(x))", "student()"),
  list("lognormal gp(x)", "bf(yl ~ gp(x))", "lognormal()"),
  list("poisson gp(x)", "bf(cnt ~ gp(x))", "poisson()"),
  list("bernoulli gp(x)", "bf(bin ~ gp(x))", "bernoulli()"),
  list("mixture mu1 ~ gp(x)", "bf(ymx ~ 1, mu1 ~ gp(x), mu2 ~ 1)", "mixture(gaussian, gaussian)")
)

# the fit at an outer point: its conditional modes come from the tape there
fit_at <- function(fit, p) {
  o <- fit$obj
  o$fn(p)
  pl <- o$env$parList(p)
  fit$estimates[names(pl)] <- pl[names(pl)]
  fit$cache <- new.env(parent = emptyenv())
  fit
}
read1 <- function(sf, fit, sdat, code, rtab) {
  pars <- E$stan_pars_from_fit(fit, sdat, code, rtab)
  up <- rstan::unconstrain_pars(sf, pars)
  lp <- rstan::log_prob(sf, up, adjust_transform = FALSE, gradient = FALSE)
  e <- fit$obj$env
  # the joint at the point fit_at() or the caller left in the tape, with the map
  # of any other block; an exact gp() adds nothing to logJ
  ours <- -as.numeric(e$f(e$last.par)) + attr(pars, "logJ")
  list(const = as.numeric(lp) - ours, ours = ours, sf = sf, up = up, pars = pars)
}

run1 <- function(i, attr_mode) {
  sp <- spell[[i]]
  fam <- if (length(sp) > 2) sp[[3]] else "gaussian()"
  mkf <- function(ns) eval(parse(text = sp[[2]]), new.env(parent = ns))
  famf <- function(ns) eval(parse(text = fam), new.env(parent = ns))
  fns <- asNamespace("frmtmb")
  bns <- asNamespace("brms")
  t0 <- Sys.time()
  fit <- suppressWarnings(frm(mkf(fns), data = d, family = famf(fns)))
  bform <- mkf(bns)
  bfam <- famf(bns)
  prior <- E$brms_flat_prior(bform, data = d, family = bfam)
  code <- brms::make_stancode(bform, data = d, family = bfam, prior = prior)
  if (attr_mode) {
    stopifnot(grepl("cov[n, n] += 1e-12;", code, fixed = TRUE))
    code <- gsub("cov[n, n] += 1e-12;", "cov[n, n] += 1e-12 * square(sdgp);",
                 code, fixed = TRUE)
  }
  sdat <- E$brms_standata(bform, data = d, family = bfam, prior = prior)
  rtab <- E$brms_ranef_table(bform, d, bfam, prior)
  sf <- suppressMessages(rstan::sampling(E$brms_stan_model(code), data = sdat,
                                         chains = 0))
  p0 <- fit$opt$par
  fit$obj$fn(p0)
  r0 <- read1(sf, fit, sdat, code, rtab)
  grad <- rstan::grad_log_prob(sf, r0$up, adjust_transform = FALSE)
  idx <- E$brms_inner_index(sf, r0$pars)
  gi <- max(abs(grad[idx]))
  cs <- vapply(c(0.05, -0.05), function(s) {
    read1(sf, fit_at(fit, p0 + s), sdat, code, rtab)$const
  }, 0)
  tol <- 1e-6 * max(1, abs(r0$ours))
  bks <- Filter(function(b) identical(b$covstruct, "gp"), fit$frame$re_blocks)
  sds <- vapply(bks, function(b) exp(fit$estimates$theta[b$theta_idx[1L]]), 0)
  cat(sprintf(paste0("NCGPL %2d %-26s %s ours %.4f const %.6g steps %.6g %.6g ",
                     "tol %.2e over %.3g grad %.2e conv %d sdgp %s %.0fs\n"),
              i, sp[[1]], if (attr_mode) "attr" else "brms", r0$ours, r0$const,
              cs[1], cs[2], tol, max(abs(c(r0$const, cs))) / tol, gi,
              fit$opt$convergence,
              paste(format(sds, digits = 3), collapse = ","),
              as.numeric(Sys.time() - t0, units = "secs")))
}
if (!length(args)) {
  for (i in seq_along(spell)) cat(i, spell[[i]][[1]], "\n")
  quit("no")
}
i <- as.integer(args[1])
tryCatch(run1(i, length(args) > 1L), error = function(e) {
  cat(sprintf("NCGPL %2d %-26s ERROR %s\n", i, spell[[i]][[1]],
              conditionMessage(e)))
})
