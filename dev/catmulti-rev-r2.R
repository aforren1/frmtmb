# Reviewer, fix round 2: transformed-term df, update(), messages under
# NA, the SE predicate, and drop1() with transformed terms on both
# builds. Arg: lane | base. Log dev/catmulti-rev-out/r2-<arm>.txt
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-catmulti-lib", "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
setwd("C:/Users/adf44/source/r/frmtmb-wt-catmulti")
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
`%||%` <- function(a, b) if (is.null(a)) b else a
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(5); n <- 400
dh <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1, g = factor(sample(letters[1:8], n, TRUE)))
E <- cbind(0, 0.3 + 0.8 * dh$x, -0.4 - 0.5 * log(dh$x2)); P <- exp(E) / rowSums(exp(E))
dh$y <- factor(apply(P, 1, function(p) sample(c("a", "b", "c"), 1, prob = p)))
dh$yc <- as.character(dh$y)
dh$z <- 1 + dh$x + log(dh$x2) + rnorm(n)
if (arm == "lane") {
  fs <- list(log = y ~ log(x2), poly = y ~ poly(x, 2) + log(x2), scale = y ~ scale(x) + I(x2^2),
             ns = y ~ splines::ns(x, 3), interaction = y ~ x * log(x2))
  for (nm in names(fs)) {
    fn <- q(frm(bf(fs[[nm]]), family = categorical(refcat = NA), data = dh))
    fa <- q(frm(bf(fs[[nm]]), family = categorical(), data = dh))
    cat(sprintf("%-12s df %s vs %d  AIC %.4f vs %.4f  anova Chi Df %s\n", nm, format(attr(logLik(fn), "df")),
                attr(logLik(fa), "df"), AIC(fn), AIC(fa),
                format(tryCatch(anova(fa, fn)[2, "Chi Df"], error = function(e) conditionMessage(e)))))
  }
  f1 <- q(frm(bf(y ~ x), family = categorical(refcat = NA), data = dh))
  f2 <- q(update(f1, . ~ . + log(x2)))
  f3 <- q(update(f1, . ~ . + log(x2), newdata = dh[1:300, ]))
  r2 <- q(frm(bf(y ~ x + log(x2)), family = categorical(), data = dh))
  r3 <- q(frm(bf(y ~ x + log(x2)), family = categorical(), data = dh[1:300, ]))
  cat("update no newdata: df", attr(logLik(f2), "df"), "ref", attr(logLik(r2), "df"),
      "AIC", AIC(f2), AIC(r2), "\n")
  cat("update newdata 300 rows: df", attr(logLik(f3), "df"), "ref", attr(logLik(r3), "df"),
      "nobs", nobs(f3), "AIC", AIC(f3), AIC(r3), "\n")
  fr <- q(frm(bf(y ~ x + (1 | p | g)), family = categorical(refcat = "b"), data = dh))
  cat("non-NA fit stores refcat_first_outer:", !is.null(fr$frame[["refcat_first_outer"]]),
      "| names in frame identical to a gaussian fit's set? ",
      identical(sort(names(fr$frame)), sort(names(q(frm(bf(y ~ x + (1 | p | g)), family = categorical(), data = dh))$frame))), "\n")
  fu <- q(update(fr, family = categorical(refcat = NA)))
  cat("update to NA: df", attr(logLik(fu), "df"), "stored", !is.null(fu$frame$refcat_first_outer$lens), "\n")
  # messages: a character response is coerced once, not twice
  cm <- function(rc) { m <- character(0); withCallingHandlers(suppressWarnings(
    frm(bf(yc ~ x), family = categorical(refcat = rc), data = dh)),
    message = function(c) { m <<- c(m, conditionMessage(c)); invokeRestart("muffleMessage") }); m }
  cat("messages refcat NA:", length(cm(NA)), " refcat default:", length(cm(NULL)), "\n")
  # SE predicate
  w <- function(expr) { out <- character(0); withCallingHandlers(suppressMessages(expr),
    warning = function(cw) { out <<- c(out, conditionMessage(cw)); invokeRestart("muffleWarning") }); out }
  wn <- w(frm(bf(y ~ x), family = categorical(refcat = NA), data = dh))
  cat("NA plain: sentence", any(grepl("Under refcat = NA", wn)), "\n")
  pri <- c(set_prior("normal(0, 1)", class = "Intercept", dpar = "mua"), set_prior("normal(0, 1)", class = "Intercept", dpar = "mub"),
           set_prior("normal(0, 1)", class = "Intercept", dpar = "muc"), set_prior("normal(0, 1)", class = "b", dpar = "mua"),
           set_prior("normal(0, 1)", class = "b", dpar = "mub"), set_prior("normal(0, 1)", class = "b", dpar = "muc"))
  dz <- dh; dz$g <- factor(sample(letters[1:8], n, TRUE))
  wp <- w(frm(bf(y ~ x + (1 | g)), family = categorical(refcat = NA), data = dz, prior = pri))
  cat("NA with priors and (1 | g) (sd likely at 0):", length(wp), "warnings; sentence", any(grepl("Under refcat = NA", wp)), "\n")
  cat(substr(wp, 1, 250), sep = "\n")
}
cat("\n== drop1 with transformed terms ==\n")
for (f in list(z ~ x + log(x2), z ~ poly(x, 2), z ~ scale(x) + x2, z ~ I(x^2) + x2, z ~ x + x2,
               z ~ x + splines::ns(x2, 3), z ~ x + (1 | g))) {
  fit <- q(frm(bf(f), data = dh))
  r <- tryCatch({ d <- drop1(fit, test = "Chisq"); paste("ok, rows", paste(rownames(d), collapse = " | ")) },
                error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 110)))
  cat(sprintf("%-28s %s\n", deparse(f), r))
}
fit <- q(frm(bf(z ~ x + log(x2)), data = dh))
r <- tryCatch(rownames(drop1(fit, scope = ~ x)), error = function(e) paste("ERROR:", conditionMessage(e)))
cat("drop1(scope = ~ x) on z ~ x + log(x2):", r, "\n")
