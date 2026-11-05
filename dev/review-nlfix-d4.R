# Fix round 2, items 4, 5, 8: the disc body and the standard errors; messages; bernoulli and allfit.
source("dev/review-nlfix-common.R")
ns <- asNamespace("frmtmb")
set.seed(14001); n <- 300
d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
d$y <- ordered(cut(0.8 * d$x + rlogis(n, scale = exp(-0.3 * d$z)), c(-Inf, -1, 0, 1, Inf), labels = 1:4))
sh <- function(lab, r) {
  if (iserr(r)) { cat(lab, "ERROR", substr(r$msg, 1, 200), "\n"); return(invisible(r)) }
  se <- tt(sqrt(diag(as.matrix(vcov(r)))))
  cat(sprintf("%-36s ll %.5f code %s\n   beta %s betad %s\n   se %s\n   W[%s]\n   M[%s]\n", lab, as.numeric(logLik(r)), r$opt$convergence, paste(format(r$estimates$beta, digits = 4), collapse = ","), paste(format(r$estimates$betad, digits = 4), collapse = ","),
              if (iserr(se)) se$msg else paste(names(se), format(se, digits = 3), collapse = " "), substr(paste(attr(r, "w"), collapse = " || "), 1, 230), substr(paste(attr(r, "m"), collapse = " || "), 1, 230)))
  invisible(r)
}
cat("--- linear disc with an intercept: the family's own check\n")
sh("y ~ x, disc ~ z", tt(frm(bf(y ~ x, disc ~ z), family = cumulative(), data = d)))
sh("y ~ x, disc ~ 0 + z", tt(frm(bf(y ~ x, disc ~ 0 + z), family = cumulative(), data = d)))
cat("--- a disc body with a constant\n")
r1 <- sh("mu linear, disc d0 + d1 z", tt(frm(bf(y ~ x) + nlf(disc ~ d0 + d1 * z) + lf(d0 ~ 1) + lf(d1 ~ 1), family = cumulative(), data = d)))
r2 <- sh("mu body b x, disc d0 + d1 z", tt(frm(bf(y ~ 1) + nlf(mu ~ b * x) + nlf(disc ~ d0 + d1 * z) + lf(b ~ 1) + lf(d0 ~ 1) + lf(d1 ~ 1), family = cumulative(), data = d, start = list(beta = c("b_(Intercept)" = 0.5)))))
r3 <- sh("mu 1 + b x (fixed level), same disc", tt(frm(bf(y ~ 1) + nlf(mu ~ 1 + b * x) + nlf(disc ~ d0 + d1 * z) + lf(b ~ 1) + lf(d0 ~ 1) + lf(d1 ~ 1), family = cumulative(), data = d)))
r4 <- sh("mu exp(b x), same disc", tt(frm(bf(y ~ 1) + nlf(mu ~ exp(b * x)) + nlf(disc ~ d0 + d1 * z) + lf(b ~ 1) + lf(d0 ~ 1) + lf(d1 ~ 1), family = cumulative(), data = d, start = list(beta = c("b_(Intercept)" = 0.5)))))
r5 <- sh("disc d1 z only (no constant)", tt(frm(bf(y ~ x) + nlf(disc ~ d1 * z) + lf(d1 ~ 1), family = cumulative(), data = d)))
if (!iserr(r1)) {
  cat("the flat direction by hand: d0 + c, thresholds and slope scaled by exp(-c); objective at c = 0, 0.5, -0.5\n")
  est <- r1$estimates; nll <- ns$build_objective(r1$frame)
  print(names(est)); print(est[c("beta", "betad", "tau_raw")])
  th <- cumsum(c(est$tau_raw[1], exp(est$tau_raw[-1])))
  for (cc in c(0, 0.5, -0.5)) {
    t1 <- est; t1$betad[1] <- est$betad[1] + cc; t1$beta <- est$beta * exp(-cc)
    th2 <- th * exp(-cc); t1$tau_raw <- c(th2[1], log(diff(th2)))
    cat(" c", cc, "objective", format(as.numeric(nll(t1)), digits = 12), "\n")
  }
  cat("zero variance along a direction of the fit's Hessian?\n")
  print(tt(summary(r1)$fixed))
}
if (!iserr(r3)) {
  cat("mu 1 + b x: is d0 determined? profile the objective along d0 by refits\n")
  for (st in c(-0.5, 0.5)) {
    rr <- tt(frm(bf(y ~ 1) + nlf(mu ~ 1 + b * x) + nlf(disc ~ d0 + d1 * z) + lf(b ~ 1) + lf(d0 ~ 1) + lf(d1 ~ 1), family = cumulative(), data = d, start = list(betad = setNames(c(st, 0), names(r3$estimates$betad)[1:2]))))
    cat(" start d0", st, ":", if (iserr(rr)) rr$msg else paste(format(as.numeric(logLik(rr)), digits = 10), "d0", format(rr$estimates$betad[1], digits = 5)), "\n")
  }
}
cat("brms on a disc with an intercept:\n")
b <- tryCatch(as.character(brms::stancode(brms::bf(y ~ x, disc ~ z), data = d, family = brms::cumulative())), error = function(e) conditionMessage(e)); cat(if (nchar(b) < 400) b else "builds", "\n")
b2 <- tryCatch(as.character(brms::stancode(brms::bf(y ~ x) + brms::nlf(disc ~ d0 + d1 * z) + brms::lf(d0 ~ 1) + brms::lf(d1 ~ 1), data = d, family = brms::cumulative())), error = function(e) conditionMessage(e)); cat(if (nchar(b2) < 400) b2 else "body builds", "\n")
cat("\n--- item 5, class and text\n")
set.seed(11005); dq <- data.frame(x = round(rnorm(n), 3)); dq$yo <- ordered(cut(-0.8 * dq$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = 1:4)); dq$yn <- -0.8 * dq$x + rnorm(n)
cl <- character(0); txt <- character(0)
r <- withCallingHandlers(frm(bf(yo ~ b^2 * x * (1 - b^2), b ~ 1, nl = TRUE), family = cumulative(), data = dq), message = function(m) { cl <<- c(cl, class(m)[1:2]); txt <<- c(txt, conditionMessage(m)); invokeRestart("muffleMessage") })
cat("ll", as.numeric(logLik(r)), " classes", cl, "\n", txt, "\n")
cat("a true one: a + b^2 x, negative slope\n")
txt <- character(0); r <- withCallingHandlers(frm(bf(yn ~ a + b^2 * x, a + b ~ 1, nl = TRUE), data = dq), message = function(m) { txt <<- c(txt, conditionMessage(m)); invokeRestart("muffleMessage") }); cat(as.numeric(logLik(r)), "\n", txt, "\n")
cat("\n--- the saddle warning: its class and text, slope 0\n")
d0 <- d; set.seed(14002); d0$yn <- rnorm(n); d0$y0 <- residuals(lm(yn ~ x, d0)) + 1
wn <- NULL; r <- withCallingHandlers(frm(bf(y0 ~ a + b1 * b2 * x, a + b1 + b2 ~ 1, nl = TRUE), data = d0), warning = function(w) { wn <<- c(wn, list(w)); invokeRestart("muffleWarning") })
cat("ll", as.numeric(logLik(r)), " linear", as.numeric(logLik(frm(y0 ~ x, data = d0))), " warnings", length(wn), "\n"); for (w in wn) cat(class(w)[1:2], "|", conditionMessage(w), "\n")
cat("\n--- item 8: bernoulli, allfit, a located prior\n")
set.seed(9202); dk <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3)); dk$y <- ordered(cut(1.5 * exp(0.6 * dk$x) + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
set.seed(9205); dk$yg <- 1.5 * exp(0.6 * dk$x) + rnorm(n, 0, 0.5); dk$yc <- rpois(n, 1.5 * exp(0.6 * dk$x)); dk$yzi <- ifelse(runif(n) < 0.25, 0L, rpois(n, 1.5 * exp(0.6 * dk$x))); dk$yhu <- ifelse(runif(n) < 0.25, 0L, 1L + rpois(n, 1.5 * exp(0.6 * dk$x))); dk$ypo <- rgamma(n, 3, 3 / (1.5 * exp(0.6 * dk$x))); dk$ybi <- rbinom(n, 1, plogis(1.5 * exp(0.6 * dk$x) - 2))
fb <- bf(ybi ~ b1 * exp(b2 * dk_x), b1 + b2 ~ 1, nl = TRUE)
fb <- bf(ybi ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)
a <- frm(fb, family = bernoulli(), data = dk); cat("default", as.numeric(logLik(a)), "code", a$opt$convergence, "\n")
print(names(formals(frm_allfit)))
af <- tt(frm_allfit(a)); if (iserr(af)) cat(af$msg, "\n") else print(af)
af2 <- tt(frm_allfit(a, start_from_mle = FALSE)); if (iserr(af2)) cat(substr(af2$msg, 1, 200), "\n") else print(af2)
s1 <- tt(frm(fb, family = bernoulli(), data = dk, start = list(beta = c("b1_(Intercept)" = 1)))); cat("start b1 = 1:", as.numeric(logLik(s1)), "\n")
s2 <- tt(frm(fb, family = bernoulli(), data = dk, prior = set_prior("normal(1, 100)", nlpar = "b1"))); cat("located prior at b1 = 1: logLik", as.numeric(logLik(s2)), " beta", format(s2$estimates$beta, digits = 4), " messages:", substr(attr(s2, "m"), 1, 100), "\n")
s3 <- tt(frm(fb, family = bernoulli(), data = dk, method = "x"));
cat("\n--- the help\n")
hl <- readLines("man/par_template.Rd"); i <- grep("second|allfit", hl); cat(hl[sort(unique(c(i, i + 1, i - 1)))], sep = "\n")
hf <- readLines("man/frm.Rd"); i <- grep("second one|allfit\\}", hf); cat(hf[sort(unique(c(i, i + 1, i - 1)))][1:8], sep = "\n")
