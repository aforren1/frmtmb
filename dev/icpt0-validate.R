# Lane icpt0: brms's reserved `Intercept` (`0 + Intercept`) and
# bf(center = FALSE), validated against brms 2.23.0 and against
# independent R log posteriors written from brms's generated Stan code.
#
#   Rscript dev/icpt0-validate.R > dev/icpt0-validate.out 2>&1
#
# Library order: the lane build, then the round's base, deps, conda.
# Seed 1 throughout; every number in dev/icpt0-findings.md comes from the
# output of this script.
.libPaths(c(Sys.getenv("FRMTMB_LIB", "/opt/rlib/lane-icpt0"),
            "/opt/rlib/base", "/opt/rlib/deps", "/opt/r/lib/R/library"))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
cat("brms", as.character(packageVersion("brms")), "\n\n")

set.seed(1)
n <- 80
d <- data.frame(x = rnorm(n, 3), z = rnorm(n),
                f = factor(rep(letters[1:4], n / 4)))
d$y <- 1 + 0.5 * d$x + 0.4 * as.numeric(d$f) + rnorm(n, sd = exp(0.2 * d$z))
d$cnt <- rpois(n, exp(0.2 + 0.2 * d$x))
d$yo <- factor(cut(d$y, quantile(d$y, 0:4 / 4), include.lowest = TRUE,
                   labels = FALSE), ordered = TRUE)

## A. the design: frmtmb's X against brms's make_standata() X --------------
cat("== A. design matrices against brms::make_standata()\n")
frame_of <- function(f) frm(f, data = d, dry_run = "frame")
lp_X <- function(fr, dpar) {
  for (lp in fr$linpreds) if (identical(lp$dpar, dpar)) return(lp$X)
}
cases <- list(
  list("0 + Intercept + x", y ~ 0 + Intercept + x, "mu", "X"),
  list("0 + x + Intercept", y ~ 0 + x + Intercept, "mu", "X"),
  list("-1 + Intercept + x", y ~ -1 + Intercept + x, "mu", "X"),
  list("x + Intercept - 1", y ~ x + Intercept - 1, "mu", "X"),
  list("0 + Intercept + f", y ~ 0 + Intercept + f, "mu", "X"),
  list("0 + Intercept + x * f", y ~ 0 + Intercept + x * f, "mu", "X"),
  list("0 + x * f + Intercept", y ~ 0 + x * f + Intercept, "mu", "X"),
  list("0 + Intercept + x:f", y ~ 0 + Intercept + x:f, "mu", "X"),
  list("0 + Intercept", y ~ 0 + Intercept, "mu", "X"),
  list("sigma ~ 0 + Intercept + z", list(y ~ x, sigma ~ 0 + Intercept + z),
       "sigma", "X_sigma"),
  list("nl a ~ 0 + Intercept + x",
       list(y ~ a * exp(b * z), a ~ 0 + Intercept + x, b ~ 1, nl = TRUE),
       "a", "X_a"),
  list("center = FALSE", list(y ~ x * f, center = FALSE), "mu", "X"))
mkbf <- function(fun, spec) {
  if (inherits(spec, "formula")) fun(spec) else do.call(fun, spec)
}
for (cs in cases) {
  Xb <- brms::make_standata(mkbf(brms::bf, cs[[2]]), data = d)[[cs[[4]]]]
  Xf <- lp_X(frame_of(mkbf(frmtmb::bf, cs[[2]])), cs[[3]])
  cn <- sub("^[(]Intercept[)]$", "Intercept", colnames(Xf))
  same_names <- identical(cn, colnames(Xb))
  diff <- if (same_names) max(abs(as.matrix(Xf) - unclass(Xb))) else NA
  cat(sprintf("%-26s brms cols: %-40s names identical: %s  max|dX| = %s\n",
              cs[[1]], paste(colnames(Xb), collapse = ","), same_names,
              format(diff)))
}

## B. the prior table: class and coef rows against brms::get_prior() ----------
cat("\n== B. b / Intercept rows of the prior table against brms::get_prior()\n")
rows <- function(tab) {
  tab <- as.data.frame(tab)
  tab <- tab[tab$class %in% c("b", "Intercept"), ]
  sort(paste(tab$class, tab$coef, tab$dpar, tab$nlpar, tab$resp, sep = "|"))
}
pcases <- list(
  list("0 + Intercept + x", y ~ 0 + Intercept + x),
  list("0 + Intercept + f", y ~ 0 + Intercept + f),
  list("sigma ~ 0 + Intercept + z", list(y ~ x, sigma ~ 0 + Intercept + z)),
  list("center = FALSE", list(y ~ x, center = FALSE)),
  list("center = FALSE + lf(center)", NULL))
for (cs in pcases) {
  if (is.null(cs[[2]])) {
    fb <- brms::bf(y ~ x, center = FALSE) + brms::lf(sigma ~ z, center = FALSE)
    ff <- frmtmb::bf(y ~ x, center = FALSE) +
      frmtmb::lf(sigma ~ z, center = FALSE)
  } else {
    fb <- mkbf(brms::bf, cs[[2]])
    ff <- mkbf(frmtmb::bf, cs[[2]])
  }
  rb <- rows(brms::get_prior(fb, data = d))
  rf <- rows(get_prior(ff, data = d))
  cat(sprintf("%-28s identical: %s  (%d rows)\n", cs[[1]],
              identical(rb, rf), length(rb)))
  if (!identical(rb, rf)) print(list(brms = rb, frmtmb = rf))
}
fb <- brms::mvbf(brms::bf(y ~ 0 + Intercept + x), brms::bf(cnt ~ x),
                 rescor = FALSE)
ff <- frmtmb::mvbf(frmtmb::bf(y ~ 0 + Intercept + x), frmtmb::bf(cnt ~ x))
cat(sprintf("%-28s identical: %s\n", "mvbf",
            identical(rows(brms::get_prior(fb, data = d)),
                      rows(get_prior(ff, data = d)))))

## C. ML: an identity, checked at full precision ----------------------------
cat("\n== C. ML: 0 + Intercept is the model 1 + x (an identity)\n")
mlpair <- function(lab, f0, f1, fam = gaussian()) {
  a <- frm(f0, family = fam, data = d)
  b <- frm(f1, family = fam, data = d)
  ca <- fixef(a)[, 1]
  cb <- fixef(b)[rownames(fixef(a)), 1]
  cat(sprintf(paste("%-30s logLik %.10f vs %.10f  dlogLik %.3g",
                    "max|dcoef| %.3g  bitwise: %s\n"),
              lab, logLik(a), logLik(b), as.numeric(logLik(a) - logLik(b)),
              max(abs(ca - cb)), identical(unname(ca), unname(cb))))
  invisible(a)
}
mlpair("gaussian 0 + Intercept + x", bf(y ~ 0 + Intercept + x), bf(y ~ x))
mlpair("gaussian 0 + x + Intercept", bf(y ~ 0 + x + Intercept), bf(y ~ x))
mlpair("gaussian f", bf(y ~ 0 + Intercept + x + f), bf(y ~ x + f))
mlpair("sigma ~ 0 + Intercept + z", bf(y ~ x, sigma ~ 0 + Intercept + z),
       bf(y ~ x, sigma ~ z))
mlpair("poisson 0 + Intercept + x*f", bf(cnt ~ 0 + Intercept + x * f),
       bf(cnt ~ x * f), poisson())
mlpair("nl a ~ 0 + Intercept + x",
       bf(y ~ a + exp(b) * z, a ~ 0 + Intercept + x, b ~ 1, nl = TRUE),
       bf(y ~ a + exp(b) * z, a ~ 1 + x, b ~ 1, nl = TRUE))
mlpair("center = FALSE", bf(y ~ x * f, center = FALSE), bf(y ~ x * f))
# and against the classical reference
a <- frm(bf(y ~ 0 + Intercept + x + f), data = d)
l <- lm(y ~ x + f, data = d)
cat(sprintf("%-30s max|coef - lm| %.3g (lm scale %.3g)\n", "vs lm()",
            max(abs(fixef(a)[, 1] - coef(l))), max(abs(coef(l)))))
a <- frm(bf(cnt ~ 0 + Intercept + x * f) + poisson(), data = d)
g <- glm(cnt ~ x * f, family = poisson(), data = d)
gn <- sub("[(]Intercept[)]", "Intercept", names(coef(g)))
cat(sprintf("%-30s max|coef - glm| %.3g  dlogLik %.3g\n", "vs glm()",
            max(abs(fixef(a)[gn, 1] - coef(g))),
            as.numeric(logLik(a) - logLik(g))))

## D. MAP: the prior class, against brms's Stan code written in R -------------
cat("\n== D. MAP against an independent log posterior\n")
# brms 2.23.0, y ~ 0 + Intercept + x with prior(normal(0, s), class = b):
#   lprior += normal_lpdf(b | 0, s);            b = (b_Intercept, b_x)
#   target += normal_id_glm_lpdf(Y | X, 0, b, sigma);   X = [1, x]
# and y ~ 1 + x with prior(normal(m, s), class = Intercept):
#   lprior += normal_lpdf(Intercept | m, s);    Intercept = b0 + mean(x) b1
# sigma has no prior in these fits, and the mode of a flat density does
# not depend on its scale, so the reference optimizes over log sigma.
X <- cbind(1, d$x)
nlp <- function(th, lprior) {
  b <- th[1:2]
  -(sum(dnorm(d$y, drop(X %*% b), exp(th[3]), log = TRUE)) + lprior(b))
}
ref_map <- function(lprior) {
  o <- optim(c(coef(lm(y ~ x, d)), 0), nlp, lprior = lprior,
             method = "BFGS", control = list(reltol = 1e-14, maxit = 1000))
  o$par[1:2]
}
map_row <- function(lab, fit, lprior) {
  est <- fixef(fit)[c("Intercept", "x"), 1]
  ref <- ref_map(lprior)
  se <- fixef(fit)[c("Intercept", "x"), 2]
  cat(sprintf(paste("%-44s frmtmb (%.6f, %.6f)",
                    "reference (%.6f, %.6f)  max|d|/se %.2g\n"),
              lab, est[1], est[2], ref[1], ref[2], max(abs(est - ref) / se)))
}
s <- 0.1
map_row("0 + Intercept + x, class b normal(0, 0.1)",
        frm(bf(y ~ 0 + Intercept + x), data = d,
            prior = set_prior("normal(0, 0.1)", class = "b")),
        function(b) sum(dnorm(b, 0, s, log = TRUE)))
map_row("center = FALSE, class b normal(0, 0.1)",
        frm(bf(y ~ x, center = FALSE), data = d,
            prior = set_prior("normal(0, 0.1)", class = "b")),
        function(b) sum(dnorm(b, 0, s, log = TRUE)))
map_row("1 + x, class b normal(0, 0.1) (slope only)",
        frm(bf(y ~ x), data = d,
            prior = set_prior("normal(0, 0.1)", class = "b")),
        function(b) dnorm(b[2], 0, s, log = TRUE))
map_row("0 + Intercept + x, b coef Intercept normal(0, 0.1)",
        frm(bf(y ~ 0 + Intercept + x), data = d,
            prior = set_prior("normal(0, 0.1)", class = "b",
                              coef = "Intercept")),
        function(b) dnorm(b[1], 0, s, log = TRUE))
map_row("1 + x, class Intercept normal(0, 0.1) (centered)",
        frm(bf(y ~ x), data = d,
            prior = set_prior("normal(0, 0.1)", class = "Intercept")),
        function(b) dnorm(b[1] + mean(d$x) * b[2], 0, s, log = TRUE))
r <- tryCatch(frm(bf(y ~ 0 + Intercept + x), data = d,
                  prior = set_prior("normal(0, 1)", class = "Intercept")),
              error = function(e) conditionMessage(e))
cat("class Intercept on 0 + Intercept:", r, "\n")

## E. ordinal under center = FALSE: thresholds not centered -------------------
cat("\n== E. cumulative(), bf(center = FALSE): threshold prior at tau\n")
# brms: with center = FALSE the ordinal program keeps X uncentered, so a
# class "Intercept" prior reads the thresholds themselves; by default it
# reads tau - mean(x) * b. Reference: cumulative logit written out.
# frmtmb stores (tau_1, log increments), Stan's `ordered` map, and adds
# that map's log-Jacobian (?set_prior, "Ordinal thresholds"), so its
# mode is the mode in the internal coordinates: the reference adds
# sum(log(diff(tau))) as well.
yo <- as.integer(d$yo)
K <- max(yo)
nlp_ord <- function(th, ctr) {
  tau <- th[1:(K - 1)]
  b <- th[K]
  eta <- d$x * b
  cut <- c(-Inf, tau, Inf)
  p <- plogis(cut[yo + 1] - eta) - plogis(cut[yo] - eta)
  lp <- tau - if (ctr) mean(d$x) * b else 0
  if (any(diff(tau) <= 0)) return(Inf)
  -(sum(log(p)) + sum(dnorm(lp, 0, 0.5, log = TRUE)) + sum(log(diff(tau))))
}
for (ctr in c(TRUE, FALSE)) {
  fit <- frm(bf(yo ~ x, center = if (!ctr) FALSE), family = cumulative(),
             data = d,
             prior = set_prior("normal(0, 0.5)", class = "Intercept"))
  fe <- fixef(fit)[, 1]
  o <- optim(c(-1, 0, 1, 0), nlp_ord, ctr = ctr, method = "Nelder-Mead",
             control = list(reltol = 1e-14, maxit = 20000))
  o <- optim(o$par, nlp_ord, ctr = ctr, method = "BFGS",
             control = list(reltol = 1e-14, maxit = 2000))
  cat(sprintf(paste0("center = %-5s frmtmb    %s\n",
                     "               reference %s  max|d| %.2g\n"),
              ctr, paste(sprintf("%.5f", fe), collapse = " "),
              paste(sprintf("%.5f", o$par), collapse = " "),
              max(abs(fe - o$par))))
}

## F. brms's generated Stan code for the prior, as read --------------------
cat("\n== F. brms Stan code: where the class b prior lands\n")
show <- function(lab, ...) {
  sc <- strsplit(brms::make_stancode(..., data = d), "\n")[[1]]
  cat("--", lab, "\n")
  cat(grep("lprior [+]=|Intercept|means_X|vector.K", sc, value = TRUE),
      sep = "\n")
}
show("y ~ 0 + Intercept + x, class b", y ~ 0 + Intercept + x,
     prior = brms::prior(normal(0, 1), class = b))
show("bf(y ~ x, center = FALSE), class b", brms::bf(y ~ x, center = FALSE),
     prior = brms::prior(normal(0, 1), class = b))
show("sigma ~ 0 + Intercept + z, class b dpar sigma",
     brms::bf(y ~ x, sigma ~ 0 + Intercept + z),
     prior = brms::prior(normal(0, 1), class = b, dpar = sigma))
