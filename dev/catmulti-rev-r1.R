# Reviewer, fix round 1: loglik_df() under refcat = NA on hard data
# (transformations, subset, weights, mi(), missing rows, large n, data
# removed), resp_from_stan_name() routing, the SE sentence, the
# two-category transform refusal, conditional_effects on categorical.
# Lane library. Log dev/catmulti-rev-out/r1.txt
.libPaths(c("C:/Users/adf44/source/r/wt-catmulti-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
setwd("C:/Users/adf44/source/r/frmtmb-wt-catmulti")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
source("dev/catmulti-rev-data.R")
sec <- function(t) cat("\n========", t, "========\n")
CAP <- function(expr) {
  w <- character(0); m <- character(0)
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "err")),
    warning = function(cw) { w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning") },
    message = function(cm) { m <<- c(m, conditionMessage(cm)); invokeRestart("muffleMessage") })
  list(value = r, warnings = w, messages = m)
}
show <- function(nm, x) {
  v <- x$value
  if (inherits(v, "err")) cat(nm, ": ERROR:", v, "\n") else { cat(nm, ":\n"); print(v) }
  if (length(x$warnings)) cat("  WARNINGS:", substr(x$warnings, 1, 400), sep = "\n    ")
  if (length(x$messages)) cat("  MESSAGES:", substr(x$messages, 1, 300), sep = "\n    ")
}
dfpair <- function(nm, f_na, f_ref) {
  a <- CAP(logLik(f_na)); b <- logLik(f_ref)
  cat(sprintf("%-28s NA df %s (why: %s) | ref df %d | logLik diff %.2e | AIC %s vs %.4f | msgs %d warns %d\n",
              nm, format(attr(a$value, "df")), attr(a$value, "df_why") %||% "-",
              attr(b, "df"), as.numeric(a$value) - as.numeric(b),
              format(AIC(f_na)), AIC(f_ref), length(a$messages), length(a$warnings)))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
q <- function(expr) suppressWarnings(suppressMessages(expr))

sec("review numbers")
dc <- rev_cat_data(23, 3, n = 300)
fn <- q(frm(bf(y ~ x), family = categorical(refcat = NA), data = dc))
fa <- frm(bf(y ~ x), family = categorical(), data = dc)
dfpair("plain", fn, fa)
show("anova(fa, fn)", CAP(anova(fa, fn)))
d3 <- rev_cat_data(21, 3); set.seed(1); d3$z <- rnorm(nrow(d3))
fn3 <- q(frm(bf(y ~ x + z), family = categorical(refcat = NA), data = d3))
show("drop1 NA", CAP(drop1(fn3, test = "Chisq")))
cat("BIC", BIC(fn), BIC(fa), " df.residual", df.residual(fn), df.residual(fa), "\n")

sec("hard data for the reference frame")
set.seed(5)
n <- 400
dh <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1, w = runif(n, 0.5, 2),
                 keep = runif(n) > 0.1, g = factor(sample(letters[1:10], n, TRUE)))
E <- cbind(0, 0.3 + 0.8 * dh$x, -0.4 - 0.5 * log(dh$x2))
P <- exp(E) / rowSums(exp(E))
dh$y <- factor(apply(P, 1, function(p) sample(c("a", "b", "c"), 1, prob = p)))
dh$xm <- dh$x; dh$xm[sample(n, 30)] <- NA
dh$ymiss <- dh$y; dh$ymiss[sample(n, 15)] <- NA
cases <- list(
  transforms = list(f = bf(y ~ poly(x, 2) + log(x2)), d = dh),
  scale_in_formula = list(f = bf(y ~ scale(x) + I(x2^2)), d = dh),
  subset = list(f = bf(y | subset(keep) ~ x), d = dh),
  weights = list(f = bf(y | weights(w) ~ x), d = dh),
  na_rows = list(f = bf(ymiss ~ xm), d = dh),
  re_cor = list(f = bf(y ~ x + (1 | p | g)), d = dh),
  re_nocor = list(f = bf(y ~ x + (1 || g)), d = dh)
)
for (nm in names(cases)) {
  cs <- cases[[nm]]
  a <- CAP(frm(cs$f, family = categorical(refcat = NA), data = cs$d))
  b <- CAP(frm(cs$f, family = categorical(), data = cs$d))
  if (inherits(a$value, "err") || inherits(b$value, "err")) {
    cat(nm, ": fit ERROR:", if (inherits(a$value, "err")) a$value else b$value, "\n"); next
  }
  dfpair(nm, a$value, b$value)
}
# mi(): a missing predictor modeled in mvbf (multivariate: expect NA with reason)
a <- CAP(frm(bf(y ~ mi(xm)) + bf(xm | mi() ~ x2) + set_rescor(FALSE),
             family = list(categorical(refcat = NA), gaussian()), data = dh))
if (!inherits(a$value, "err")) {
  ll <- CAP(logLik(a$value))
  cat("mi() mvbf: df", format(attr(ll$value, "df")), "why:", attr(ll$value, "df_why") %||% "-", "\n")
  show("AIC mi()", CAP(AIC(a$value)))
} else cat("mi() mvbf fit ERROR:", a$value, "\n")
# data removed from the workspace after the fit
local({
  dloc <- dh
  f <- q(frm(bf(y ~ log(x2)), family = categorical(refcat = NA), data = dloc))
  rm(dloc)
  r <- CAP(logLik(f))
  cat("data removed: df", format(attr(r$value, "df")), "\n")
})
# update() must not carry a stale cached df
f1 <- q(frm(bf(y ~ x), family = categorical(refcat = NA), data = dh))
invisible(logLik(f1))
f2 <- q(update(f1, . ~ . + log(x2)))
f2r <- frm(bf(y ~ x + log(x2)), family = categorical(), data = dh)
dfpair("update adds a term", f2, f2r)
f1c <- f1; f1c$cache <- f1$cache
cat("copy shares cache env:", identical(f1c$cache, f1$cache), "\n")
# own formula -> NA with reason; drop1/anova refuse
fo <- q(frm(bf(y ~ x, mua ~ 1), family = categorical(refcat = NA), data = dh))
ll <- logLik(fo)
cat("own formula: df", format(attr(ll, "df")), "why:", attr(ll, "df_why"), "\n")
show("AIC own formula", CAP(AIC(fo)))
show("drop1 own formula", CAP(drop1(fo)))
# priors
pri <- c(set_prior("normal(0, 2)", class = "Intercept", dpar = "mua"),
         set_prior("normal(0, 2)", class = "Intercept", dpar = "mub"),
         set_prior("normal(0, 2)", class = "Intercept", dpar = "muc"),
         set_prior("normal(0, 2)", class = "b", dpar = "mua"),
         set_prior("normal(0, 2)", class = "b", dpar = "mub"),
         set_prior("normal(0, 2)", class = "b", dpar = "muc"))
fp <- q(frm(bf(y ~ x), family = categorical(refcat = NA), data = dh, prior = pri))
cat("with priors: df", attr(logLik(fp), "df"), "\n")
# multinomial and dirichlet NA
dm <- rev_mult_data(31)
dfpair("multinomial", q(frm(bf(Yc | trials(N) ~ x), family = multinomial(refcat = NA), data = dm)),
       frm(bf(Yc | trials(N) ~ x), family = multinomial(), data = dm))
# timing on large data
set.seed(9)
nb <- 50000
db <- data.frame(x = rnorm(nb), z = rnorm(nb))
Eb <- cbind(0, 0.3 + 0.8 * db$x, -0.4 - 0.5 * db$z)
Pb <- exp(Eb) / rowSums(exp(Eb))
u <- runif(nb); cp <- t(apply(Pb, 1, cumsum))
db$y <- factor(c("a", "b", "c")[1 + rowSums(u > cp[, 1:2])])
tf <- system.time(fb <- q(frm(bf(y ~ x + z), family = categorical(refcat = NA), data = db)))
tl <- system.time(llb <- logLik(fb))
tl2 <- system.time(logLik(fb))
cat(sprintf("n = %d: fit %.2f s, first logLik %.2f s, cached %.3f s, df %d\n",
            nb, tf[["elapsed"]], tl[["elapsed"]], tl2[["elapsed"]], attr(llb, "df")))

sec("SE sentence fires only under refcat = NA")
w_na <- CAP(frm(bf(y ~ x), family = categorical(refcat = NA), data = dh))$warnings
cat("NA fit warnings mention 'refcat = NA':", any(grepl("refcat = NA", w_na)), "\n")
dsep <- data.frame(x = c(rnorm(50, -3), rnorm(50, 3)))
dsep$y <- factor(ifelse(dsep$x > 0, "b", "a")); dsep$y[1:3] <- "c"
w_other <- c(CAP(frm(bf(y ~ x + (1 | g)), family = categorical(), data = dh))$warnings,
             CAP(frm(bf(y ~ x), family = categorical(refcat = "b"), data = dsep))$warnings,
             CAP(frm(bf(Y ~ x), family = dirichlet(refcat = NA), data = local({
               set.seed(2); d <- data.frame(x = rnorm(200)); G <- matrix(rgamma(600, 3), 200)
               d$Y <- G / rowSums(G); colnames(d$Y) <- c("a", "b", "c"); d })))$warnings)
cat("non-NA (and dirichlet NA) warnings:\n"); cat(substr(w_other, 1, 300), sep = "\n")
cat("non-categorical warnings mention 'refcat = NA':",
    any(grepl("refcat = NA", w_other[seq_len(length(w_other))])), "\n")

sec("resp_from_stan_name routing")
set.seed(3)
dr <- data.frame(x = rnorm(200))
dr$Y_a <- 1 + dr$x + rnorm(200); dr$Ya <- -1 + rnorm(200); dr$y.a <- rnorm(200)
show("mvbf(Y_a, Ya) fit", CAP(names(frm(mvbf(bf(Y_a ~ x), bf(Ya ~ x)), data = dr)$spec$responses)))
show("mvbf(Y_a, y.a) fit", CAP(names(frm(mvbf(bf(Y_a ~ x), bf(y.a ~ x)), data = dr)$spec$responses)))
show("brms stancode mvbf(Y_a, Ya)", CAP(substr(brms::make_stancode(
  brms::bf(Y_a ~ x) + brms::bf(Ya ~ x) + brms::set_rescor(FALSE), data = dr), 1, 50)))
dr$Yb_c <- 2 + 0.5 * dr$x + rnorm(200)
fm2 <- frm(mvbf(bf(Y_a ~ x), bf(Yb_c ~ x)), data = dr,
           prior = set_prior("normal(5, 0.01)", class = "b", resp = "Ybc"))
print(round(fixef(fm2), 3))
fm3 <- frm(mvbf(bf(Y_a ~ x), bf(Yb_c ~ x)), data = dr,
           prior = set_prior("normal(5, 0.01)", class = "b", resp = "Yb_c"))
cat("raw and Stan spellings give the same fit:", isTRUE(all.equal(fixef(fm2), fixef(fm3))), "\n")
show("resp = 'nope'", CAP(frm(mvbf(bf(Y_a ~ x), bf(Yb_c ~ x)), data = dr,
                                prior = set_prior("normal(5, 0.01)", class = "b", resp = "nope"))))
bdp <- as.data.frame(brms::default_prior(brms::bf(Y_a ~ x) + brms::bf(Yb_c ~ x) + brms::set_rescor(FALSE), data = dr))
cat("brms resp names:", unique(bdp$resp), "\n")
show("brms validate_prior resp = 'Yb_c'", CAP(nrow(brms::validate_prior(
  brms::set_prior("normal(5, 0.01)", class = "b", resp = "Yb_c"),
  brms::bf(Y_a ~ x) + brms::bf(Yb_c ~ x) + brms::set_rescor(FALSE), data = dr))))
show("brms validate_prior resp = 'Ybc'", CAP(nrow(brms::validate_prior(
  brms::set_prior("normal(5, 0.01)", class = "b", resp = "Ybc"),
  brms::bf(Y_a ~ x) + brms::bf(Yb_c ~ x) + brms::set_rescor(FALSE), data = dr))))

sec("two-category transform, ce on categorical")
d2 <- dc; d2$y <- factor(ifelse(d2$y == "a", "a", "b"))
f2c <- q(frm(bf(y ~ x), family = categorical(), data = d2))
ds2 <- q(frm_sample(f2c, chains = 1, iter = 300, warmup = 150, refresh = 0, seed = 9))
show("2cat transform", CAP(posterior_linpred(ds2, transform = TRUE)))
show("2cat dpar mub", CAP(dim(posterior_linpred(ds2, dpar = "mub"))))
show("2cat transform dpar mub range", CAP(range(posterior_linpred(ds2, dpar = "mub", transform = TRUE))))
ce <- CAP(conditional_effects(fa, "x", categorical = TRUE)[[1]])
if (!inherits(ce$value, "err")) {
  cv <- ce$value
  nd <- cv[cv$cats__ == "b", , drop = FALSE]
  ff <- fitted(fa, newdata = nd)
  cat("ce vs fitted(newdata) estimate rel:", max(abs(nd$estimate__ - ff[, "Estimate", "P(Y = b)"]) /
                                                 ff[, "Estimate", "P(Y = b)"]),
      " se rel:", max(abs(nd$se__ - ff[, "Est.Error", "P(Y = b)"]) / ff[, "Est.Error", "P(Y = b)"]), "\n")
  cat("ce bounds within [0, 1]:", all(cv$lower__ >= 0 & cv$upper__ <= 1), "\n")
} else show("ce", ce)
show("ce categorical band profile", CAP(conditional_effects(fa, "x", categorical = TRUE, band = "profile")))
show("ce categorical method predict", CAP(conditional_effects(fa, "x", categorical = TRUE, method = "predict")))
fre <- frm(bf(y ~ x + (1 | g)), family = categorical(refcat = "b"), data = dh)
ce2 <- CAP(conditional_effects(fre, "x", categorical = TRUE)[[1]])
cat("ce with RE refcat b:", if (inherits(ce2$value, "err")) ce2$value else
  paste(dim(ce2$value), collapse = " x "), "levels", levels(ce2$value$cats__), "\n")
sec("print and nits")
print(multinomial(refcat = "y2")); print(categorical(refcat = NA))
show("integer refusal", CAP(frm(bf(Yc | trials(N) ~ x), family = multinomial(),
                               data = within(dm, Yc[1, 1] <- 0.5))))
cat("DONE\n")
