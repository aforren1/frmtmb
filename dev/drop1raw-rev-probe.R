# Reviewer probes of lane drop1raw: update() semantics, rows and values
# of fit_refit_data(), offsets, bases under update(), priors.
# Arg: lane | base. Log: dev/drop1raw-rev-out/probe-<arm>.txt
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
tr <- function(lab, expr) {
  r <- tryCatch(q(expr), error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-58s ERROR %s\n", lab,
                substr(gsub("\n", " ", conditionMessage(r)), 1, 160)))
  } else if (inherits(r, "frmtmb_fit")) {
    cat(sprintf("%-58s nobs %d n_obs %d logLik %.10f\n", lab, nobs(r),
                r$frame$n_obs, as.numeric(logLik(r))))
  } else {
    cat(sprintf("%-58s %s\n", lab, paste(format(r), collapse = " ")))
  }
  invisible(r)
}
ll <- function(f) as.numeric(logLik(f))
set.seed(5)
n <- 200
d <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1, x3 = rnorm(n),
                g = factor(sample(letters[1:8], n, TRUE)),
                expo = runif(n, 0.5, 3), wt = runif(n, 0.5, 2))
d$z <- 1 + d$x + log(d$x2) + rnorm(8)[d$g] + rnorm(n)
d$cnt <- rpois(n, d$expo * exp(0.3 + 0.4 * d$x + 0.3 * log(d$x2)))
d$x3[c(3, 17, 40)] <- NA
d0 <- d

cat("\n== 1. update(): which data, which rows ==\n")
fit <- q(frm(bf(z ~ x + log(x2) + x3), data = d))
tr("fit z ~ x + log(x2) + x3 (3 NA in x3)", fit)
tr("update(fit, ~ . - x3)", update(fit, ~ . - x3))
tr("lm(z ~ x + log(x2)) nobs for reference", nobs(lm(z ~ x + log(x2), d)))
d$x2 <- d$x2 * 10
d <- d[-(1:10), ]
d$w <- rnorm(nrow(d))
tr("after d changed: update(fit)", update(fit))
tr("after d changed: update(fit, family = student())",
   update(fit, family = student()))
tr("after d changed: update(fit, ~ . + w)  [w only in d]",
   update(fit, ~ . + w))
wv <- rnorm(n)
tr("after d changed: update(fit, ~ . + wv) [wv global, 200]",
   update(fit, ~ . + wv))
tr("frm on the changed d, z ~ x + log(x2) + x3 + w (reference)",
   frm(bf(z ~ x + log(x2) + x3 + w), data = d))
tr("update(fit, data = d0)", update(fit, data = d0))
tr("update(fit, newdata = d0[1:50, ])", update(fit, newdata = d0[1:50, ]))
tr("update(fit, na.action = na.exclude)",
   update(fit, na.action = na.exclude))
d <- d0
fx <- q(frm(bf(z ~ x + log(x2) + x3), data = d, na.action = na.exclude))
ux <- q(update(fx, ~ . - x))
tr("na.exclude fit: length(residuals(fit))", length(residuals(fx)))
tr("na.exclude fit: length(residuals(update(fit, ~ . - x)))",
   length(residuals(ux)))
tr("na.exclude: lm length(residuals(update(lm, ~ . - x)))",
   length(residuals(update(lm(z ~ x + log(x2) + x3, d,
                              na.action = na.exclude), ~ . - x))))

cat("\n== 2. update() and bases when na.action dropped rows ==\n")
for (f in list(z ~ splines::ns(x2, 3) + x3, z ~ poly(x2, 2) + x3,
               z ~ scale(x) + x3, z ~ splines::ns(x2, 3) * g + x3)) {
  fb <- q(frm(bf(f), data = d))
  u <- q(update(fb))
  us <- q(update(fb, family = student()))
  fs <- q(frm(bf(f), family = student(), data = d))
  cat(sprintf("%-34s fit %.8f | update(fit) %.8f | diff %.2e\n",
              deparse1(f), ll(fb), ll(u), ll(u) - ll(fb)))
  cat(sprintf("%-34s student: frm %.8f | update(family=) %.8f\n", "",
              ll(fs), ll(us)))
  cat(sprintf("%-34s coef 2 fit %.6f update %.6f\n", "",
              unname(fixef(fb)[2, 1]), unname(fixef(u)[2, 1])))
  dt <- q(drop1(fb))
  u3 <- q(update(fb, ~ . - x3))
  cat(sprintf("%-34s drop1 AIC(-x3) %.8f | AIC(update(fit, ~ . - x3)) %.8f\n",
              "", dt["x3", "AIC"], AIC(u3)))
}

cat("\n== 3. offsets, weights ==\n")
fo <- q(frm(bf(cnt ~ x + log(x2) + offset(log(expo))), family = poisson(),
            data = d))
tr("poisson + offset(log(expo)) fit", fo)
tr("names(raw_vars)", paste(names(fo$frame$raw_vars), collapse = ","))
tr("drop1(fo)", drop1(fo)$AIC)
g0 <- glm(cnt ~ x + log(x2) + offset(log(expo)), poisson, d)
tr("glm drop1 AIC", drop1(g0)$AIC)
tr("update(fo, ~ . - x)", update(fo, ~ . - x))
fo_gone <- local({
  dg <- d
  q(frm(bf(cnt ~ x + log(x2) + offset(log(expo))), family = poisson(),
        data = dg))
})
tr("offset fit, data gone: drop1()", drop1(fo_gone)$AIC)
tr("offset fit, data gone: update(fit)", update(fo_gone))
expo <- rep(1, n)
tr("offset fit, data gone, global expo = 1: drop1()",
   drop1(fo_gone)$AIC)
tr("offset fit, data gone, global expo = 1: update(fit)", update(fo_gone))
rm(expo)
fw <- q(frm(bf(z | weights(wt) ~ x + log(x2) + x3), data = d))
tr("weights(wt) fit", fw)
tr("weights(wt): update(fit)", update(fw))
tr("weights(wt): drop1 AIC", drop1(fw)$AIC)

cat("\n== 4. fit_refit_data(): rows, values, classes ==\n")
frd <- getFromNamespace("fit_refit_data", "frmtmb")
set.seed(7)
m <- 120
e <- data.frame(x = rnorm(m), x2 = rexp(m) + 0.1,
                f = factor(sample(c("a", "b", "c"), m, TRUE),
                           levels = c("a", "b", "c", "unused")),
                ch = sample(c("p", "q"), m, TRUE),
                dt = as.Date("2020-01-01") + sample(0:400, m),
                unused_na = NA_real_, stringsAsFactors = FALSE)
e$z <- 1 + e$x + log(e$x2) + rnorm(m)
e$s <- rbinom(m, 10, 0.4); e$fl <- 10 - e$s
e$x2[c(5, 9)] <- NA
e$unused_na[1:20] <- 1
rownames(e) <- paste0("r", seq_len(m))
cmp <- function(lab, f, data, ...) {
  fit <- tryCatch(q(frm(bf(f), data = data, ...)), error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("%-40s FIT ERROR %s\n", lab,
                substr(conditionMessage(fit), 1, 120)))
    return(invisible())
  }
  rd <- frd(fit)
  keep <- rownames(data) %in% rownames(rd)
  vars <- intersect(names(data), names(rd))
  same <- vapply(vars, function(v) {
    identical(unname(rd[[v]]), unname(data[[v]][keep])) ||
      (is.factor(data[[v]]) && identical(as.character(rd[[v]]),
                                         as.character(data[[v]][keep])))
  }, NA)
  cls <- vapply(vars, function(v) class(rd[[v]])[1], "")
  cat(sprintf("%-40s rows %d (expect %d) | rownames ok %s | equal: %s\n",
              lab, nrow(rd), sum(keep),
              identical(rownames(rd), rownames(data)[keep]),
              paste(sprintf("%s=%s(%s)", vars, same, cls), collapse = " ")))
  u <- tryCatch(q(update(fit)), error = function(e) e)
  dr <- tryCatch(q(drop1(fit)), error = function(e) e)
  cat(sprintf("%-40s update ll diff %s | drop1 %s\n", "",
              if (inherits(u, "error")) conditionMessage(u) else
                format(ll(u) - ll(fit)),
              if (inherits(dr, "error")) substr(conditionMessage(dr), 1, 90)
              else "ok"))
}
cmp("na.omit, NA in x2 only", z ~ x + log(x2), e)
cmp("NAs in an unused column", z ~ x + log(x2) + f, e)
cmp("factor with unused level", z ~ log(x2) + f, e)
cmp("factor unused level, drop_unused_levels=F", z ~ log(x2) + f, e,
    drop_unused_levels = FALSE)
cmp("character column", z ~ log(x2) + ch, e)
cmp("Date column via as.numeric", z ~ log(x2) + as.numeric(dt), e)
cmp("poly(x2, 2, raw = TRUE) matrix column",
    z ~ poly(x2, 2, raw = TRUE) + x, e)
cmp("binomial s | trials(10) ~ log(x2)", s | trials(10) ~ log(x2) + x, e,
    family = binomial())
cmp("cbind(s, fl) ~ log(x2) binomial", cbind(s, fl) ~ log(x2) + x, e,
    family = binomial())
e$sb <- e$x > -0.5
cmp("subset(sb) in the formula", z | subset(sb) ~ log(x2) + x, e)
fl <- q(frm(bf(z ~ log(x2) + f), data = e))
rd <- frd(fl)
cat("factor levels in refit data:", paste(levels(rd$f), collapse = ","),
    "| fit levels:", paste(levels(fl$frame$data_frame$f), collapse = ","),
    "\n")

cat("\n== 5. multivariate, rows per response ==\n")
e$z2 <- e$x + rnorm(m)
e$z2[c(11, 12)] <- NA
e$s1 <- e$x > -1
mv <- tryCatch(q(frm(mvbf(bf(z | subset(s1) ~ log(x2)), bf(z2 ~ x),
                          rescor = FALSE), data = e)),
               error = function(e) e)
if (inherits(mv, "error")) {
  cat("mv fit ERROR", conditionMessage(mv), "\n")
} else {
  tr("mvbf subset + NA in z2: fit", mv)
  tr("nobs per resp", c(nobs(mv, resp = "z"), nobs(mv, resp = "z2")))
  tr("update(mv)", update(mv))
  tr("nrow fit_refit_data(mv, 'data')", nrow(frd(mv, "data")))
}
mv2 <- tryCatch(q(frm(mvbf(bf(z ~ log(x2)), bf(z2 ~ x), rescor = FALSE),
                      data = e)), error = function(e) e)
if (!inherits(mv2, "error")) {
  tr("mvbf NA in z2 (no subset): fit", mv2)
  tr("update(mv2)", update(mv2))
}

cat("\n== 6. priors: logLik of a MAP fit, drop1 with priors ==\n")
dd <- d[!is.na(d$x3), ]
fp <- q(frm(bf(z ~ x + log(x2) + x3), data = dd,
            prior = set_prior("normal(0, 1)", class = "b")))
f0 <- q(frm(bf(z ~ x + log(x2) + x3), data = dd))
b <- fixef(fp)[, 1]
X <- model.matrix(~ x + log(x2) + x3, dd)
sig <- if ("sigma" %in% names(fp$estimates)) NA else NA
mu <- drop(X %*% b[colnames(X) %in% colnames(X)][seq_len(ncol(X))])
sd_hat <- tryCatch(as.numeric(summary(fp)$spec_pars[1, 1]),
                   error = function(e) NA)
cat("MAP logLik", sprintf("%.8f", ll(fp)), "| ML logLik",
    sprintf("%.8f", ll(f0)), "\n")
print(q(fixef(fp)))
print(tryCatch(q(VarCorr(fp)), error = function(e) conditionMessage(e)))
cat("drop1 on the MAP fit:\n")
print(tryCatch(q(drop1(fp, test = "Chisq")), error = function(e) e))
fps <- q(frm(bf(z ~ scale(x) + log(x2) + x3), data = dd,
             prior = set_prior("normal(0, 0.1)", class = "b",
                               coef = "scalex")))
cat("prior on the KEPT term scale(x), drop1 messages:\n")
withCallingHandlers(
  tryCatch(drop1(fps), error = function(e) cat("ERROR", conditionMessage(e), "\n")),
  message = function(m) {
    cat("MESSAGE:", conditionMessage(m)); invokeRestart("muffleMessage")
  }, warning = function(w) invokeRestart("muffleWarning"))
fps_ns <- q(frm(bf(z ~ scale(x) + log(x2) + x3), data = dd,
                prior = set_prior("normal(0, 0.1)", class = "b",
                                  coef = "scalex")))
u_x3 <- q(update(fps_ns, ~ . - x3))
cat("update(~ . - x3) keeps the scale(x) prior: coef",
    sprintf("%.6f", fixef(u_x3)["scalex", 1]), "vs drop1 refit by hand: ")
pri <- fps_ns$prior
fb <- getFromNamespace("freeze_bases", "frmtmb")
nb <- fps_ns$bform
nb$formula <- update(nb$formula, . ~ . - x3)
nb <- fb(nb, fps_ns$frame$predvar_map)
cat(deparse1(nb$formula), "\n")
