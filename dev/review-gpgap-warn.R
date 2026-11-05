# Review of lane gpgap, item 3: when does frmtmb_gp_nugget fire, when does it
# stay silent, and can gp_nugget_message() itself break a fit? Seed 5501.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-warn.R warn.txt
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(5501)
n <- 60
x <- round(runif(n, 0, 6), 3)
e <- rnorm(n)
d <- data.frame(x = x, z = round(runif(n, 0, 6), 3),
                g = factor(rep(1:6, length.out = n)),
                f = factor(rep(c("a", "b"), length.out = n)),
                e = e, w = runif(n, 0.5, 1.5))

run <- function(tag, expr) {
  warned <- character()
  msgs <- character()
  t0 <- proc.time()[["elapsed"]]
  fit <- withCallingHandlers(
    tryCatch(expr, error = function(er) {
      cat(sprintf("WARN %-34s ERROR: %s\n", tag,
                  substr(conditionMessage(er), 1, 200)))
      NULL
    }),
    warning = function(w) {
      warned <<- c(warned, class(w)[1])
      if (inherits(w, "frmtmb_gp_nugget")) msgs <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  if (is.null(fit)) return(invisible(NULL))
  sg <- tryCatch(paste(signif(range(as.numeric(sigma(fit))), 4),
                       collapse = ".."), error = function(e) "n/a")
  bks <- Filter(function(k) k$covstruct %in% c("gp", "hsgp"), fit$frame$re_blocks)
  sds <- vapply(bks, function(b) signif(exp(fit$estimates$theta[b$theta_idx][[1]]), 4), 0)
  cat(sprintf("WARN %-34s nug %-5s sigma %s sdgp %s ll %.3f conv %d others: %s\n",
              tag, "frmtmb_gp_nugget" %in% warned, sg,
              paste(sds, collapse = ","), -fit$opt$objective, fit$opt$convergence,
              paste(setdiff(unique(warned), "frmtmb_gp_nugget"), collapse = "/")))
  if (length(msgs)) cat("     msg:", substr(msgs, 1, 130), "\n")
  invisible(fit)
}
low <- function(s = 1e-3) sin(2 * d$x) + s * d$e

# 1. families with a sigma
d$y <- low()
fg <- run("gaussian low", frm(bf(y ~ gp(x)), family = gaussian(), data = d))
run("gaussian gp(x, k = 15) low", frm(bf(y ~ gp(x, k = 15)), data = d))
run("gaussian gp(x, k = 20) low", frm(bf(y ~ gp(x, k = 20)), data = d))
run("gaussian gp(x, k = 30) low", frm(bf(y ~ gp(x, k = 30)), data = d))
run("gaussian gp(x, k = 60) low", frm(bf(y ~ gp(x, k = 60)), data = d))
run("student low", frm(bf(y ~ gp(x)), family = student(), data = d))
d$yp <- exp(d$y)
run("lognormal low", frm(bf(yp ~ gp(x)), family = lognormal(), data = d))
run("shifted_lognormal low",
    frm(bf(yp ~ gp(x)), family = shifted_lognormal(), data = d))
run("skew_normal low", frm(bf(y ~ gp(x)), family = skew_normal(), data = d))
run("exgaussian low", frm(bf(y ~ gp(x)), family = exgaussian(), data = d))
run("asym_laplace low", frm(bf(y ~ gp(x)), family = asym_laplace(), data = d))
d$yh <- d$yp; d$yh[1:5] <- 0
run("hurdle_lognormal low", frm(bf(yh ~ gp(x)), family = hurdle_lognormal(), data = d))
d$yg <- d$yp
run("Gamma (no sigma)", frm(bf(yg ~ gp(x)), family = Gamma(link = "log"), data = d))

# 2. forms of sigma
run("sigma ~ z", frm(bf(y ~ gp(x), sigma ~ z), data = d))
d$yh1 <- sin(2 * d$x) + ifelse(seq_len(n) == 7, 1e-4, 0.3) * d$e
d$one <- factor(ifelse(seq_len(n) == 7, "a", "b"))
run("sigma low on 1 row only", frm(bf(yh1 ~ gp(x), sigma ~ one), data = d))
d$wsm <- 1e-3
run("sigma ~ offset(log(w)), w = 1e-3",
    frm(bf(y ~ gp(x), sigma ~ 0 + offset(log(wsm))), data = d))
run("sigma ~ 1 + offset(log(w))",
    frm(bf(y ~ gp(x), sigma ~ 1 + offset(log(wsm))), data = d))
run("sigma = 0.0005 constant", frm(bf(y ~ gp(x), sigma = 0.0005), data = d))
run("sigma = 0.3 constant", frm(bf(y ~ gp(x), sigma = 0.3), data = d))
run("sigma ~ (1|g)", frm(bf(y ~ gp(x), sigma ~ (1 | g)), data = d))
run("sigma ~ s(z)", frm(bf(y ~ gp(x), sigma ~ s(z, k = 5)), data = d))
run("sigma ~ gp(z)", frm(bf(y ~ gp(x), sigma ~ gp(z)), data = d))
d$yy <- 0.1 * d$z + d$y
run("nl body b1 ~ gp(x)", frm(bf(yy ~ b1 + b2 * z, b1 ~ gp(x), b2 ~ 1, nl = TRUE),
                              data = d))
run("by = f", frm(bf(y ~ gp(x, by = f)), data = d))
run("gp(x, z) 2d", frm(bf(y ~ gp(x, z)), data = d))
d$xr <- round(d$x + 1e-4 * seq_len(n), 4)
run("gr = FALSE", frm(bf(y ~ gp(xr, gr = FALSE)), data = d))
run("REML", frm(bf(y ~ gp(x)), data = d, REML = TRUE))
run("gp plus (1|g)", frm(bf(y ~ gp(x) + (1 | g)), data = d))
d$y2 <- 0.5 * d$z + rnorm(n, 0, 0.3)
run("mvbf rescor = FALSE",
    frm(mvbf(bf(y ~ gp(x)), bf(y2 ~ z), rescor = FALSE), data = d))
run("mvbf rescor = TRUE",
    frm(mvbf(bf(y ~ gp(x)), bf(y2 ~ z), rescor = TRUE), data = d))
run("mixture mu1 ~ gp(x)",
    frm(bf(y ~ 1, mu1 ~ gp(x), mu2 ~ 1), family = mixture(gaussian(), gaussian()),
        data = d))

# 3. se(): known noise below the floor, against the closed form at each nugget
d$se <- 3e-4
d$ys <- sin(2 * d$x) + 3e-4 * d$e
fs <- run("se(3e-4) low", frm(bf(ys | se(se) ~ gp(x)), data = d))
run("se(3e-4, sigma = TRUE) low", frm(bf(ys | se(se, sigma = TRUE) ~ gp(x)), data = d))
marg <- function(nugget, y, sdv) {
  u <- sort(unique(d$x)); D2 <- outer(u, u, "-")^2; J <- match(d$x, u)
  nll <- function(p) {
    s2 <- exp(2 * p[2])
    K <- s2 * exp(-D2 / (2 * exp(2 * p[3]))); diag(K) <- diag(K) + nugget(s2)
    V <- K[J, J]; diag(V) <- diag(V) + sdv^2
    R <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(R)) return(1e10)
    r <- backsolve(R, y - p[1], transpose = TRUE)
    0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * n * log(2 * pi)
  }
  best <- NULL
  for (st in list(c(0, 0, 0), c(0, 0.7, 0.3), c(0, 1.5, 0.5), c(0, 0.3, -0.3))) {
    o <- nlminb(st, nll, control = list(iter.max = 600, eval.max = 900))
    if (is.null(best) || o$objective < best$objective) best <- o
  }
  best
}
a <- marg(function(s2) 1e-6 * s2, d$ys, 3e-4)
b <- marg(function(s2) 1e-12, d$ys, 3e-4)
cat(sprintf(paste0("SE closed form, known se 3e-4: frm nugget sdgp %.4f rho %.4f ",
                   "ll %.3f | brms nugget sdgp %.4f rho %.4f ll %.3f\n"),
            exp(a$par[2]), exp(a$par[3]), -a$objective,
            exp(b$par[2]), exp(b$par[3]), -b$objective))

# 4. class, update, suppress, refits
if (!is.null(fg)) {
  got <- tryCatch({
    frm(bf(y ~ gp(x)), family = gaussian(), data = d)
    "none"
  }, frmtmb_gp_nugget = function(w) paste(class(w), collapse = ","))
  cat("WARN class seen by tryCatch:", got, "\n")
  run("update(fit)", update(fg))
  run("update(fit, sigma data)", update(fg, newdata = d[1:50, ]))
  run("control check_se = ignore",
      frm(bf(y ~ gp(x)), data = d, control = frmtmb_control(check_se = "ignore")))
  run("se = FALSE", frm(bf(y ~ gp(x)), data = d, se = FALSE))
  run("start = from fit", frm(bf(y ~ gp(x)), data = d, start = fg$estimates))
  w0 <- character()
  withCallingHandlers(
    drop1(fg, trace = FALSE),
    warning = function(w) {
      w0 <<- c(w0, class(w)[1])
      invokeRestart("muffleWarning")
    }, error = function(e) w0 <<- c(w0, "ERR"))
}

# 5. non-gaussian, high information: no sigma, same nugget
d$cnt <- rpois(n, exp(12 + sin(2 * d$x)))
pe <- run("poisson e^12 exact", frm(bf(cnt ~ gp(x)), family = poisson(), data = d))
pk <- run("poisson e^12 k = 40", frm(bf(cnt ~ gp(x, k = 40)), family = poisson(), data = d))
d$tr <- 1e7L
d$suc <- rbinom(n, d$tr, plogis(sin(2 * d$x)))
be <- run("binomial 1e7 exact", frm(bf(suc | trials(tr) ~ gp(x)), family = binomial(), data = d))
bk <- run("binomial 1e7 k = 40", frm(bf(suc | trials(tr) ~ gp(x, k = 40)), family = binomial(), data = d))
if (!is.null(pe) && !is.null(pk)) {
  cat(sprintf("POIS ll exact %.3f, k = 40 %.3f\n", -pe$opt$objective, -pk$opt$objective))
}

# 6. ordinary noise, a long trend: the field is huge against sigma
d$lin <- 5 * d$x + 0.3 * d$e
d$qd <- 1.5 * d$x^2 + 0.3 * d$e
d$expo <- exp(1.2 * d$x) + 0.3 * d$e
for (v in c("lin", "qd", "expo")) {
  fm <- as.formula(paste(v, "~ gp(x)"))
  fx <- run(paste("trend", v), frm(do.call(bf, list(fm)), data = d))
  a <- marg2 <- NULL
  u <- sort(unique(d$x)); D2 <- outer(u, u, "-")^2; J <- match(d$x, u)
  y <- d[[v]]
  fitm <- function(nugget) {
    nll <- function(p) {
      s2 <- exp(2 * p[2])
      K <- s2 * exp(-D2 / (2 * exp(2 * p[3]))); diag(K) <- diag(K) + nugget(s2)
      V <- K[J, J]; diag(V) <- diag(V) + exp(2 * p[4])
      R <- tryCatch(chol(V), error = function(e) NULL)
      if (is.null(R)) return(1e10)
      r <- backsolve(R, y - p[1], transpose = TRUE)
      0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * n * log(2 * pi)
    }
    best <- NULL
    for (st in list(c(mean(y), log(sd(y)), 0, -1),
                    c(mean(y), log(sd(y)) + 2, 1.5, -1),
                    c(mean(y), log(sd(y)) + 4, 2.5, -1))) {
      o <- nlminb(st, nll, control = list(iter.max = 800, eval.max = 1200))
      if (is.null(best) || o$objective < best$objective) best <- o
    }
    best
  }
  af <- fitm(function(s2) 1e-6 * s2)
  bb <- fitm(function(s2) 1e-12)
  cat(sprintf(paste0("TREND %s closed form: frm nugget sigma %.4f sdgp %.4g rho %.4g ll %.3f | ",
                     "brms nugget sigma %.4f sdgp %.4g rho %.4g ll %.3f\n"),
              v, exp(af$par[4]), exp(af$par[2]), exp(af$par[3]), -af$objective,
              exp(bb$par[4]), exp(bb$par[2]), exp(bb$par[3]), -bb$objective))
}
