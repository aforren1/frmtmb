# Item 5, 6, 7 and the density fix's accuracy. Arg picks a block; each in its own process.
a <- commandArgs(trailingOnly = TRUE)
blk <- a[1]
source("dev/review-nlfix-common.R")
`%||%` <- function(a, b) if (is.null(a)) b else a
suppressPackageStartupMessages(library(testthat))
ns <- asNamespace("frmtmb")
E <- new.env(parent = ns); E$`%||%` <- function(a, b) if (is.null(a)) b else a
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) sys.source(f, envir = E)
mkd <- function() {
  set.seed(9202); n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
  d$y <- ordered(cut(1.5 * exp(0.6 * d$x) + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
  set.seed(1); d$yh <- ordered(ifelse(runif(n) < 0.2, 0L, as.integer(d$y)), levels = 0:4)
  set.seed(9205)
  d$yg <- 1.5 * exp(0.6 * d$x) + rnorm(n, 0, 0.5)
  d$yc <- rpois(n, 1.5 * exp(0.6 * d$x))
  d
}
if (blk %in% c("pois_frm", "pois_opt", "pois_opt_guard")) {
  d <- mkd()
  fb <- brms::bf(yc ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE)
  fam <- brms::brmsfamily("poisson", link = "identity")
  if (blk == "pois_frm") {
    r <- tt(frm(bf(yc ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE), family = poisson("identity"), data = d))
    cat("frm default:", if (iserr(r)) substr(r$msg, 1, 130) else as.numeric(logLik(r)), "\n")
    for (s in list(c(.3, .3), c(-.3, .3), c(-.3, -.3), c(.3, -.3), c(1, 1))) {
      r <- tt(frm(bf(yc ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE), family = poisson("identity"), data = d, start = list(beta = c("b1_(Intercept)" = s[1], "b2_(Intercept)" = s[2]))))
      cat("frm start", s, ":", if (iserr(r)) substr(r$msg, 1, 100) else as.numeric(logLik(r)), "\n")
    }
    cat("frm arm finished, exit 0\n")
  } else {
    prior <- E$brms_flat_prior(fb, data = d, family = fam)
    code <- brms::make_stancode(fb, data = d, family = fam, prior = prior)
    sdat <- brms::make_standata(fb, data = d, family = fam, prior = prior)
    m <- E$brms_stan_model(code)
    cat(grep("target \\+=", strsplit(as.character(code), "\n")[[1]], value = TRUE)[1], "\n"); flush.console()
    for (i in 1:6) {
      cat("optimizing seed", 100 + i, "\n"); flush.console()
      r <- tryCatch(rstan::optimizing(m, data = sdat, seed = 100 + i, iter = 5000, as_vector = FALSE), error = function(e) e)
      cat(if (inherits(r, "error")) paste("ERROR", substr(conditionMessage(r), 1, 150)) else r$value, "\n"); flush.console()
    }
    cat("brms arm finished, exit 0\n")
  }
}
if (blk == "fitcheck") {
  d <- mkd()
  r <- tt(frm(bf(y ~ 0 + x) + nlf(disc ~ exp(a * z)) + lf(a ~ 1), family = cumulative(), data = d))
  cat("disc body, linear mu:", if (iserr(r)) r$msg else as.numeric(logLik(r)), "\n W:", substr(attr(r, "w"), 1, 220), "\n")
  r <- tt(frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE) + nlf(disc ~ exp(a * z)) + lf(a ~ 1), family = cumulative(), data = d))
  cat("disc body, mu body:", if (iserr(r)) substr(r$msg, 1, 200) else as.numeric(logLik(r)), "\n W:", substr(attr(r, "w"), 1, 220), "\n")
  r <- tt(frm(bf(y ~ 1) + nlf(mu ~ b1 * x) + nlf(disc ~ exp(a * z)) + lf(b1 ~ 1) + lf(a ~ 1), family = cumulative(), data = d))
  cat("nlf mu, nlf disc:", if (iserr(r)) substr(r$msg, 1, 200) else as.numeric(logLik(r)), "\n W:", substr(attr(r, "w"), 1, 220), "\n")
}
if (blk == "bbgrid") {
  dat <- function(s) { set.seed(s); n <- 300; d <- data.frame(x = round(rnorm(n), 3), ntr = 10L); d$y <- rbinom(n, 10, plogis(0.3 + 0.5 * d$x)); d }
  d <- dat(9)
  beta <- coef(glm(cbind(y, ntr - y) ~ x, binomial, d)); p <- plogis(beta[1] + beta[2] * d$x)
  ex <- function(phi) sum(vapply(seq_len(nrow(d)), function(i) {
    a <- p[i] * phi; b <- (1 - p[i]) * phi; y <- d$y[i]; n <- d$ntr[i]
    lchoose(n, y) + (if (y > 0) sum(log(a + 0:(y - 1))) else 0) + (if (n - y > 0) sum(log(b + 0:(n - y - 1))) else 0) - sum(log(phi + 0:(n - 1)))
  }, 0))
  sd <- ns$lgamma_shift_diff
  pa <- function(phi) sum(lgamma(d$ntr + 1) - lgamma(d$y + 1) - lgamma(d$ntr - d$y + 1) + sd(p * phi, d$y) + sd((1 - p) * phi, d$ntr - d$y) - sd(rep(phi, nrow(d)), d$ntr))
  for (lp in c(-3, 0, 2, 10, 20, 30, 39, 60, 100)) cat(sprintf("log phi %4d exact %.8f patched %.8f diff %.2e\n", lp, ex(exp(lp)), pa(exp(lp)), pa(exp(lp)) - ex(exp(lp))))
  print(ns$lgamma_shift_diff)
  cat("tape gradient in log phi at 39:\n")
  tp <- RTMB::MakeTape(function(lp) { "c" <- RTMB::ADoverload("c"); sum(lgamma(d$ntr + 1) - lgamma(d$y + 1) - lgamma(d$ntr - d$y + 1) + sd(p * exp(lp), d$y) + sd((1 - p) * exp(lp), d$ntr - d$y) - sd(rep(1, nrow(d)) * exp(lp), d$ntr)) }, 39)
  cat(" value", tp(39), " jac", as.numeric(tp$jacobian(39)), "\n")
}
if (blk == "ordse") {
  d <- mkd(); set.seed(66)
  f <- bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)
  r <- suppressWarnings(suppressMessages(frm(f, family = cumulative(), data = d)))
  ce <- conditional_effects(r, effects = "x", categorical = TRUE, resolution = 5)[[1]]
  print(names(ce)); print(head(ce[, intersect(c("x", "effect2__", "estimate__", "se__", "lower__", "upper__"), names(ce))], 3))
  cat("\ncategory probabilities: finite differences of the whole parameter vector with vcov\n")
  nd <- data.frame(x = sort(unique(ce$x)), z = 0)
  est <- r$estimates
  flat <- r$opt$par
  V <- as.matrix(vcov(r))
  print(dimnames(V)[[1]])
  probs <- function(b1, b2, th) {
    eta <- b1 * exp(b2 * nd$x)
    cu <- cbind(0, sapply(th, function(t) plogis(t - eta)), 1)
    t(apply(cu, 1, diff))
  }
  th <- summary(r)$fixed[grep("Intercept\\[", rownames(summary(r)$fixed)), 1]
  bb <- fixef(r)[, 1]
  nm <- rownames(V)
  theta <- c(bb[grep("b1", names(bb))], bb[grep("b2", names(bb))], th)
  names(theta) <- c("b1", "b2", names(th))
  # V is on the reported scale? compare both ways; summary's Est.Error on thresholds is the check
  Jg <- function(cat) {
    h <- 1e-5
    do.call(cbind, lapply(seq_along(theta), function(j) {
      tp <- theta; tm <- theta; tp[j] <- tp[j] + h; tm[j] <- tm[j] - h
      (probs(tp[1], tp[2], tp[-(1:2)])[, cat] - probs(tm[1], tm[2], tm[-(1:2)])[, cat]) / (2 * h)
    }))
  }
  ord <- match(c(grep("b1", nm, value = TRUE), grep("b2", nm, value = TRUE), grep("Intercept\\[", nm, value = TRUE)), nm)
  if (!length(grep("Intercept\\[", nm))) ord <- seq_len(nrow(V))
  Vo <- V[ord, ord]
  cat("V dim", dim(V), " order", nm[ord], "\n")
  cat("SE of thresholds vcov:", round(sqrt(diag(Vo))[3:5], 5), " summary:", round(summary(r)$fixed[1:3, 2], 5), "\n")
  cats <- sort(unique(ce$effect2__))
  for (k in seq_along(cats)) {
    J <- Jg(k)
    se <- sqrt(rowSums((J %*% Vo) * J))
    cc <- ce[ce$effect2__ == cats[k], ]
    mine_p <- probs(theta[1], theta[2], theta[-(1:2)])[, k]
    cat(sprintf("cat %s: |estimate - hand| %.2e | se_hand %s | se_ce %s | max ratio %.4f\n", cats[k], max(abs(cc$estimate__ - mine_p)), paste(round(se, 5), collapse = ","),
                paste(round(cc$se__ %||% ((cc$upper__ - cc$lower__) / (2 * qnorm(.975))), 5), collapse = ","),
                max(abs((cc$se__ %||% ((cc$upper__ - cc$lower__) / (2 * qnorm(.975)))) / se - 1))))
  }
  cat("\nthe linear spelling\n")
  fl <- frm(y ~ 0 + x, family = cumulative(), data = d); fn <- suppressWarnings(frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), family = cumulative(), data = d))
  a <- conditional_effects(fn, effects = "x", categorical = TRUE, resolution = 20)[[1]]; b <- conditional_effects(fl, effects = "x", categorical = TRUE, resolution = 20)[[1]]
  cat("categorical: est max diff", max(abs(a$estimate__ - b$estimate__)), " upper", max(abs(a$upper__ - b$upper__)), " lower", max(abs(a$lower__ - b$lower__)), " largest half-width", max(b$upper__ - b$lower__) / 2, " ratio to largest se", max(abs(a$upper__ - b$upper__)) / (max(b$upper__ - b$lower__) / (2 * qnorm(.975))), "\n")
  a <- conditional_effects(fn, effects = "x", resolution = 20)[[1]]; b <- conditional_effects(fl, effects = "x", resolution = 20)[[1]]
  cat("mean: est max diff", max(abs(a$estimate__ - b$estimate__)), " upper", max(abs(a$upper__ - b$upper__)), " lower", max(abs(a$lower__ - b$lower__)), " ratio to largest se", max(abs(a$upper__ - b$upper__)) / (max(b$upper__ - b$lower__) / (2 * qnorm(.975))), "\n")
  cat("\nagainst brms at the same parameter values: categorical estimate\n")
  fxb <- suppressWarnings(E$brms_fixed_fit(brms::bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), brms::cumulative(), d, r))
  bb2 <- brms::conditional_effects(fxb, effects = "x", categorical = TRUE, resolution = 5)[[1]]
  cat("brms ce names:", paste(head(names(bb2), 10), collapse = ","), "\n")
  mine <- ce[order(ce$effect2__, ce$x), ]; theirs <- bb2[order(bb2$effect2__, bb2$x), ]
  cat("max |frm - brms| estimate:", max(abs(mine$estimate__ - theirs$estimate__)), " x equal:", isTRUE(all.equal(mine$x, theirs$x)), "\n")
  bb3 <- brms::conditional_effects(fxb, effects = "x", resolution = 5)[[1]]
  a3 <- conditional_effects(r, effects = "x", resolution = 5)[[1]]
  cat("max |frm - brms| mean display:", max(abs(a3$estimate__ - bb3$estimate__)), "\n")
  cat("\nhurdle_cumulative body\n")
  fh <- suppressWarnings(suppressMessages(frm(bf(yh ~ b1 * x, b1 ~ 1, hu ~ 1, nl = TRUE), family = hurdle_cumulative(), data = d)))
  flh <- frm(bf(yh ~ 0 + x, hu ~ 1), family = hurdle_cumulative(), data = d)
  a <- tt(conditional_effects(fh, effects = "x", categorical = TRUE, resolution = 10)[[1]]); b <- tt(conditional_effects(flh, effects = "x", categorical = TRUE, resolution = 10)[[1]])
  if (iserr(a) || iserr(b)) cat(a$msg, b$msg, "\n") else cat("hurdle upper diff", max(abs(a$upper__ - b$upper__)), "\n")
}
if (blk == "pow") {
  cat("nl_pow on the tape against Stan's, y on the tape, x data\n")
  np <- ns$nl_pow
  for (y0 in c(2, 0.5, 1e-8, 0, -1e-8, -1)) {
    x <- c(0, 2, 0, 3)
    tp <- RTMB::MakeTape(function(p) sum(np(x, p)), y0)
    tq <- RTMB::MakeTape(function(p) { xx <- c(2, 3); sum(xx^p) }, y0)
    cat(sprintf("y=%-6g value tape %s, by hand %s | grad tape %s, by hand (Stan: 0 at x = 0) %s\n", y0, format(tp(y0), digits = 10), format(sum(c(ifelse(y0 == 0, 1, if (y0 > 0) 0 else Inf), 2^y0, ifelse(y0 == 0, 1, if (y0 > 0) 0 else Inf), 3^y0)), digits = 10),
                format(as.numeric(tp$jacobian(y0)), digits = 10), format(2^y0 * log(2) + 3^y0 * log(3), digits = 10)))
  }
  cat("sign on the tape at 0:", RTMB::MakeTape(function(p) sign(p), 0)(0), " on -0:", RTMB::MakeTape(function(p) sign(-p), 0)(0), "\n")
  set.seed(9302)
  dz <- data.frame(xz = abs(rnorm(100))); dz$xz[1:5] <- 0; dz$yp <- 2 * dz$xz^0.5 + rnorm(100, 0, 0.3)
  fb <- brms::bf(yp ~ a * xz^b, a + b ~ 1, nl = TRUE); fam <- brms::brmsfamily("gaussian")
  prior <- E$brms_flat_prior(fb, data = dz, family = fam)
  code <- brms::make_stancode(fb, data = dz, family = fam, prior = prior); sdat <- brms::make_standata(fb, data = dz, family = fam, prior = prior)
  m <- E$brms_stan_model(code)
  sf <- suppressMessages(rstan::sampling(m, data = sdat, chains = 0))
  r <- tt(frm(bf(yp ~ a * xz^b, a + b ~ 1, nl = TRUE), data = dz, start = list(beta = c("a_(Intercept)" = 1, "b_(Intercept)" = 1))))
  cat("frm:", if (iserr(r)) r$msg else as.numeric(logLik(r)), format(r$estimates$beta, digits = 7), r$estimates$betad, "\n")
  s <- rstan::optimizing(m, data = sdat, seed = 1, init = list(b_a = array(1), b_b = array(1), sigma = 1), as_vector = FALSE)
  cat("brms optimizing:", s$value, "at", s$par$b_a, s$par$b_b, s$par$sigma, "\n")
  fn <- r$obj$fn; gr <- r$obj$gr
  for (pt in list(c(1, 1), c(2, 0.5), c(1.5, 0.2), c(1, 0))) {
    u <- c(pt, r$estimates$betad)
    sigma <- exp(r$estimates$betad)
    g_stan <- rstan::grad_log_prob(sf, c(pt, log(sigma)), adjust_transform = FALSE)
    lp_stan <- rstan::log_prob(sf, c(pt, log(sigma)), adjust_transform = FALSE)
    cat(sprintf("a=%g b=%g: value frm %.8f Stan %.8f | grad frm %s Stan %s\n", pt[1], pt[2], -fn(u), lp_stan, paste(format(-gr(u)[1:2], digits = 8), collapse = ","), paste(format(g_stan[1:2], digits = 8), collapse = ",")))
  }
  r0 <- tt(frm(bf(yp ~ a * xz^b, a + b ~ 1, nl = TRUE), data = dz)); cat("default start:", if (iserr(r0)) substr(r0$msg, 1, 150) else as.numeric(logLik(r0)), "\n")
  cat("a negative base, a body with x on the tape: ", {v <- tt(frm(bf(yp ~ (a + xz)^2, a ~ 1, nl = TRUE), data = dz)); if (iserr(v)) substr(v$msg, 1, 90) else as.numeric(logLik(v))}, "\n")
  fit <- tt(frm(bf(yp ~ pow(a + xz, 2), a ~ 1, nl = TRUE), data = dz)); cat("pow(tape, 2):", if (iserr(fit)) substr(fit$msg, 1, 90) else as.numeric(logLik(fit)), "\n")
  cat("newdata with zeros:", max(abs(fitted(r, newdata = dz[1:6, ])[, 1] - r$estimates$beta[1] * dz$xz[1:6]^r$estimates$beta[2])), "\n")
}
if (blk == "ext") {
  suppressPackageStartupMessages({ library(frmtmb.survival); library(frmtmb.eam) })
  cat("frmtmb.survival", dirname(find.package("frmtmb.survival")), " eam", dirname(find.package("frmtmb.eam")), "\n")
  ll <- function(r) if (iserr(r)) paste("ERR", substr(r$msg, 1, 160)) else format(as.numeric(logLik(r)), digits = 11)
  set.seed(7101); n <- 300
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  d$time <- rexp(n, exp(0.5 * d$x)); d$cc <- rbinom(n, 1, 0.2)
  cat("cox  linear", ll(tt(frm(bf(time | cens(cc) ~ x), family = cox(), data = d))), " b0 + b1 x", ll(tt(frm(bf(time | cens(cc) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = cox(), data = d))), " b1 x", ll(tt(frm(bf(time | cens(cc) ~ b1 * x, b1 ~ 1, nl = TRUE), family = cox(), data = d))), " exp body", ll(tt(frm(bf(time | cens(cc) ~ b0 + b1 * exp(b2 * x), b0 + b1 + b2 ~ 1, nl = TRUE), family = cox(), data = d, start = list(beta = c("b1_(Intercept)" = 0.5, "b2_(Intercept)" = 0.5))))), "\n")
  cat("rp   linear", ll(tt(frm(bf(time | cens(cc) ~ x), family = royston_parmar(), data = d))), " b0 + b1 x", ll(tt(frm(bf(time | cens(cc) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = royston_parmar(), data = d))), "\n")
  set.seed(7102); dw <- wiener_gng_simulate(400, mu = 1, bs = 1.4, ndt = 0.25, deadline = 1.5); dw$x <- rnorm(nrow(dw))
  cat("gng  linear", ll(tt(frm(bf(rt | dec(responded) ~ x, bias = 0.5), family = wiener_gng(deadline = 1.5), data = dw))), " body", ll(tt(frm(bf(rt | dec(responded) ~ b0 + b1 * x, b0 + b1 ~ 1, bias = 0.5, nl = TRUE), family = wiener_gng(deadline = 1.5), data = dw))), "\n")
  print(as.data.frame(frm_compat("nl", "wiener"))[, "status"]); print(as.data.frame(frm_compat("nl", "cox"))[, "status"]); print(as.data.frame(frm_compat("nl", "royston_parmar"))[, "status"]); print(as.data.frame(frm_compat("nl", "wiener_gng"))[, "status"]); print(as.data.frame(frm_compat("nl", "rdm"))[, "status"]); print(as.data.frame(frm_compat("nl", "gddm"))[, "status"])
  # brms on cox and wiener: L1 by the program's log_prob at frm()'s estimate
  d$cens <- ifelse(d$cc == 1, "right", "none")
  for (nm in c("cox")) {
    fb <- brms::bf(time | cens(cc) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)
    code <- tryCatch(brms::make_stancode(fb, data = d, family = brms::brmsfamily("cox")), error = function(e) conditionMessage(e))
    cat("brms cox body:", if (nchar(code) < 300) code else paste(grep("target \\+=", strsplit(as.character(code), "\n")[[1]], value = TRUE)[1:2], collapse = " | "), "\n")
    fb2 <- brms::bf(time | cens(cc) ~ b1 * x, b1 ~ 1, nl = TRUE)
    code2 <- tryCatch(brms::make_stancode(fb2, data = d, family = brms::brmsfamily("cox")), error = function(e) conditionMessage(e))
    cat("brms cox b1 * x builds:", nchar(code2) > 300, "\n")
    ln <- strsplit(as.character(code2), "\n")[[1]]; cat(paste(grep("Ybhaz|bhaz|mu\\[n\\]|mu =|Intercept", ln, value = TRUE)[1:12], collapse = "\n"), "\n")
  }
  cat("the seam's own check on the row text\n")
  print(substr(as.data.frame(frm_compat("nl", "cox"))$note, 1, 200))
}
