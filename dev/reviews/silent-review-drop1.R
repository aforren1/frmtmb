# Review of item 4235. Usage: Rscript silent-review-drop1.R lane|base
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r20",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-silent-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
sh <- function(lab, e) {
  r <- tryCatch(withCallingHandlers(e, message = function(m) {
    cat("  [message]", substr(conditionMessage(m), 1, 200), "\n")
    invokeRestart("muffleMessage")
  }), error = function(err) {
    cat("  ", lab, "ERROR:", substr(conditionMessage(err), 1, 300), "\n")
    NULL
  })
  if (!is.null(r)) {
    cat("==", lab, "\n")
    print(r)
  }
  invisible(r)
}
ll <- function(f) as.numeric(logLik(f))
set.seed(11)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), k = sample(1:4, n, TRUE),
                f = factor(sample(c("a", "b"), n, TRUE)),
                g = factor(sample(letters[1:8], n, TRUE)),
                o = runif(n, 0, 0.5))
d$xo <- d$x + rnorm(n, 0, 0.3)
d$sdx <- 0.3
d$y <- 1 + 0.6 * d$k + 0.4 * d$z + sin(2 * d$x) + rnorm(8)[d$g] +
  d$o + rnorm(n, 0, 0.5)

# gp(): row vs hand refit
fg <- frm(y ~ gp(x) + z, data = d)
tg <- sh("gp", drop1(fg, test = "Chisq"))
if (!is.null(tg) && "gp(x)" %in% rownames(tg)) {
  h <- frm(y ~ z, data = d)
  cat(sprintf("GP row LRT %.6f hand %.6f Df %s hand %s\n", tg["gp(x)", "LRT"],
              2 * (ll(fg) - ll(h)), tg["gp(x)", "Df"],
              attr(logLik(fg), "df") - attr(logLik(h), "df")))
}
# cs() in an ordinal family
set.seed(3)
d$yo <- cut(d$z + 0.8 * d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf),
            labels = FALSE)
fc <- frm(yo ~ cs(x) + z, family = acat(), data = d)
tc <- sh("cs", drop1(fc, test = "Chisq"))
if (!is.null(tc) && "cs(x)" %in% rownames(tc)) {
  h <- frm(yo ~ z, family = acat(), data = d)
  cat(sprintf("CS row LRT %.6f hand %.6f Df %s hand %s\n", tc["cs(x)", "LRT"],
              2 * (ll(fc) - ll(h)), tc["cs(x)", "Df"],
              attr(logLik(fc), "df") - attr(logLik(h), "df")))
  h2 <- frm(yo ~ cs(x), family = acat(), data = d)
  cat(sprintf("CS z row LRT %.6f hand %.6f\n", tc["z", "LRT"],
              2 * (ll(fc) - ll(h2))))
}
# marginality
sh("x * mo(k)", drop1(frm(y ~ x * mo(k), data = d)))
sh("s(x) + s(x, by = f)", drop1(frm(y ~ s(x) + s(x, by = f) + f, data = d)))
sh("cs(x) * z acat", drop1(frm(yo ~ cs(x) * z, family = acat(), data = d)))
# offset
fo <- frm(y ~ z + offset(o), data = d)
to <- sh("offset", drop1(fo, test = "Chisq"))
if (!is.null(to)) {
  h <- frm(y ~ 1 + offset(o), data = d)
  cat(sprintf("OFFSET z row LRT %.6f hand %.6f\n", to["z", "LRT"],
              2 * (ll(fo) - ll(h))))
}
# distributional
sh("sigma ~ x", drop1(frm(bf(y ~ z + x, sigma ~ x), data = d)))
# nonlinear
d$yn <- 2 * exp(0.3 * d$x) + rnorm(n, 0, 0.3)
sh("nl", drop1(frm(bf(yn ~ a * exp(b * x), a + b ~ 1, nl = TRUE), data = d,
                   prior = c(prior(normal(1, 2), nlpar = "a"),
                             prior(normal(0, 1), nlpar = "b")))))
# nonlinear with a linear nlpar formula holding terms
sh("nl a ~ z", drop1(frm(bf(yn ~ a * exp(b * x), a ~ z, b ~ 1, nl = TRUE),
                        data = d,
                        prior = c(prior(normal(1, 2), nlpar = "a"),
                                  prior(normal(0, 1), nlpar = "b")))))
# add1 of a new random-effect term
f1 <- frm(y ~ x, data = d)
sh("add1 (1|g)", add1(f1, ~ . + (1 | g), test = "Chisq"))
sh("add1 mo(k) + gp(z)", add1(f1, ~ . + mo(k) + z, test = "Chisq"))
# step() over a table with an NA me() row
fm <- frm(y ~ me(xo, sdx) + z + k, data = d)
st <- sh("step", step(fm, trace = 0))
if (!is.null(st)) cat("STEP kept:", deparse(st$bform$formula), "\n")
st2 <- sh("step mo", step(frm(y ~ mo(k) + z + f, data = d), trace = 0))
if (!is.null(st2)) cat("STEP2 kept:", deparse(st2$bform$formula), "\n")
# the default scope picks up a term spelled with a string argument
sh("s(x, bs = 'cr')", drop1(frm(y ~ s(x, bs = "cr") + z, data = d)))
# a term from an extension-registered special in the main formula? the
# me() of a variable used elsewhere
sh("me + me:z", drop1(frm(y ~ me(xo, sdx) * z, data = d), test = "Chisq"))
