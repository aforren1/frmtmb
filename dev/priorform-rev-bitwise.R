# Reviewer, lane wt-priorform: does every formula the base fitted still
# fit BITWISE the same on the lane build? Part A fits a spread of
# formulas through bf(), lf(), nlf(), mvbf() and cs(); part B runs every
# Rd example of core and the extensions (set.seed(1) before each) and
# keeps every frmtmb_fit it leaves behind. Output is an RDS of
# logLik, optimizer vector, objective and sdreport covariance per fit.
#   Rscript dev/priorform-rev-bitwise.R ref|lane [A|B]
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]; part <- if (length(args) > 1) args[2] else "A"
.libPaths(c(if (mode == "lane") "C:/Users/adf44/source/r/priorform-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-priorform"
Sys.setenv(FRMTMB_STAN_CACHE = file.path(root, "dev", "stan-cache"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")

digest_fit <- function(fit) {
  sdr <- tryCatch(frmtmb:::sdr_of(fit), error = function(e) NULL)
  list(ll = tryCatch(as.numeric(logLik(fit)), error = function(e) NA),
       par = fit$opt$par, obj = fit$opt$objective,
       cov = sdr$cov.fixed,
       est = fit$estimates)
}

out <- list()
if (part == "A") {
  set.seed(20260916)
  ng <- 15; per <- 12; n <- ng * per
  lev <- paste0("i", 1:ng)
  A <- diag(ng); for (i in seq(1, ng - 1, 2)) A[i, i + 1] <- A[i + 1, i] <- 0.5
  dimnames(A) <- list(lev, lev)
  d <- data.frame(g = factor(rep(lev, each = per), levels = lev),
                  x = runif(n, 0, 3), z = rnorm(n))
  d$id <- d$g; d$id2 <- d$g
  d$g1 <- factor(sample(lev, n, TRUE), levels = lev)
  d$g2 <- factor(sample(lev, n, TRUE), levels = lev)
  d$time <- factor(rep(1:4, length.out = n))
  u <- rnorm(ng, 0, .5)
  d$y <- 1 + 0.5 * d$x + u[d$g] + rnorm(n, 0, exp(-0.3 + 0.3 * d$z))
  d$ynl <- 3 * exp(-0.6 * d$x) + u[d$g] * 0.3 + rnorm(n, 0, 0.2)
  d$y1 <- d$y + rnorm(n); d$y2 <- 0.5 * d$z + u[d$g] + rnorm(n)
  lat <- 0.8 * d$x - 0.6 * d$z + u[d$g] + rlogis(n)
  d$yo <- cut(lat, c(-Inf, 0, 1.2, 2.4, Inf), labels = FALSE)
  d$cnt <- rpois(n, exp(0.2 + 0.3 * d$x)) * rbinom(n, 1, 0.8)
  d$yb <- plogis(0.3 * d$x - 0.5 + u[d$g] * 0.5 + rnorm(n, 0, .3))
  d$yb <- pmin(pmax(d$yb, 0.01), 0.99)
  d2 <- list(A = A)
  models <- list(
    lmm = quote(frm(bf(y ~ x + (x | g)) + gaussian(), data = d)),
    lmm_family_arg = quote(frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)),
    dist_sigma = quote(frm(bf(y ~ x + (1 | g), sigma ~ z) + gaussian(), data = d)),
    dist_lf = quote(frm(bf(y ~ x + (1 | g)) + lf(sigma ~ z) + gaussian(), data = d)),
    dist_const = quote(frm(bf(y ~ x, sigma = 1) + gaussian(), data = d)),
    nl_explicit_false = quote(frm(bf(y ~ x, nl = FALSE) + gaussian(), data = d)),
    nl = quote(frm(bf(ynl ~ a * exp(-b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE) +
                     gaussian(), data = d, start = list(beta = c(3, .6)))),
    nl_shared = quote(frm(bf(ynl ~ a * exp(-b * x), a + b ~ 1, nl = TRUE) +
                            gaussian(), data = d, start = list(beta = c(3, .6)))),
    nlf = quote(frm(bf(ynl ~ a, nl = TRUE) + nlf(a ~ c0 * exp(-c1 * x)) +
                      lf(c0 ~ 1, c1 ~ 1) + gaussian(), data = d,
                    start = list(beta = c(3, .6)))),
    mvbf = quote(frm(mvbf(bf(y1 ~ x + (1 | p | g)), bf(y2 ~ z + (1 | p | g))) ,
                     data = d)),
    mvbind = quote(frm(bf(mvbind(y1, y2) ~ x + (1 | g)), data = d)),
    sratio_cs = quote(frm(bf(yo ~ x + cs(z)) + sratio(), data = d)),
    acat_cs_two = quote(frm(bf(yo ~ cs(x) + cs(z)) + acat(), data = d)),
    cumul_re = quote(frm(bf(yo ~ x + z + (1 | g)) + cumulative(), data = d)),
    beta_phi_re = quote(frm(bf(yb ~ x + (1 | g), phi ~ (1 | g)) + Beta(), data = d)),
    beta_phi_re_prior = quote(frm(bf(yb ~ x + (1 | g), phi ~ (1 | g)) + Beta(),
                                  data = d,
                                  prior = set_prior("normal(0, 1)", class = "sd") +
                                    set_prior("normal(0, 2)", class = "b"))),
    beta_phi_re_prior_group = quote(frm(bf(yb ~ x + (1 | g), phi ~ (1 | g)) + Beta(),
                                  data = d,
                                  prior = set_prior("normal(0, 0.1)", class = "sd",
                                                    group = "g"))),
    zip = quote(frm(bf(cnt ~ x, zi ~ z) + zero_inflated_poisson(), data = d)),
    smooth = quote(frm(bf(y ~ s(x, k = 5) + (1 | g)) + gaussian(), data = d)),
    gp = quote(frm(bf(y ~ gp(x)) + gaussian(), data = d)),
    cs_cov = quote(frm(bf(y ~ cs(0 + time | g)) + gaussian(), data = d)),
    ar1_cov = quote(frm(bf(y ~ ar1(0 + time | g)) + gaussian(), data = d)),
    int_0slope = quote(frm(bf(y ~ x + (1 | g) + (0 + z | g)) + gaussian(), data = d)),
    dblbar = quote(frm(bf(y ~ x + (1 + z || g)) + gaussian(), data = d)),
    grcov_copy = quote(frm(bf(y ~ x + (1 | gr(id, cov = A)) + (1 | id2)) + gaussian(),
                           data = d, data2 = d2)),
    mm = quote(frm(bf(y ~ x + (1 | mm(g1, g2))) + gaussian(), data = d)),
    theta_prior = quote(frm(bf(y ~ x + (1 | g)) + gaussian(), data = d,
                            prior = set_prior("normal(-1, 0.5)", class = "theta"))),
    sigma_natural_prior = quote(frm(bf(y ~ x) + gaussian(), data = d,
                            prior = set_prior("student_t(3, 0, 1)", class = "sigma") +
                              set_prior("normal(0, 1)", coef = "x")))
  )
  for (nm in names(models)) {
    r <- tryCatch(suppressWarnings(suppressMessages(eval(models[[nm]]))),
                  error = function(e) conditionMessage(e))
    out[[nm]] <- if (is.character(r)) r else digest_fit(r)
    cat(nm, if (is.character(r)) paste("ERROR", r) else "ok", "\n")
  }
} else {
  pkgs <- c(frmtmb = ".", frmtmb.coupling = "extensions/frmtmb.coupling",
            frmtmb.eam = "extensions/frmtmb.eam",
            frmtmb.latent = "extensions/frmtmb.latent",
            frmtmb.learn = "extensions/frmtmb.learn",
            frmtmb.ode = "extensions/frmtmb.ode",
            frmtmb.spline = "extensions/frmtmb.spline")
  for (p in names(pkgs)) {
    suppressMessages(library(p, character.only = TRUE))
    rds <- list.files(file.path(root, pkgs[[p]], "man"), "[.]Rd$",
                      full.names = TRUE)
    for (rd in rds) {
      ex <- tempfile(fileext = ".R")
      tools::Rd2ex(rd, ex, commentDontrun = TRUE, commentDonttest = FALSE)
      if (!file.exists(ex)) next
      env <- new.env(parent = globalenv())
      set.seed(1)
      t0 <- Sys.time()
      err <- tryCatch({
        suppressWarnings(suppressMessages(utils::capture.output(
          sys.source(ex, envir = env))))
        ""
      }, error = function(e) conditionMessage(e))
      fits <- Filter(function(o) inherits(o, "frmtmb_fit") &&
                       !inherits(o, "frmtmb_unfitted"),
                     mget(ls(env), envir = env))
      key <- paste0(p, "::", basename(rd))
      out[[key]] <- list(err = substr(err, 1, 200),
                         fits = lapply(fits, digest_fit))
      cat(key, length(fits), "fits",
          sprintf("%.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs"))),
          if (nzchar(err)) paste("ERR:", substr(err, 1, 80)), "\n")
    }
  }
}
saveRDS(out, file.path(root, "dev", sprintf("priorform-rev-bitwise-%s-%s.rds",
                                            part, mode)))
cat("DONE\n")
