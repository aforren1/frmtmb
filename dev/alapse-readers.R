# Lane alapse, backlog 4982: what every reader of the covariance says on
# an asym_laplace() fit without random effects (seed 5, gaussian data),
# a zero_inflated_asym_laplace() fit with zi ~ x (seed 9), and a
# distributional fit (sigma ~ z, quantile ~ z).
#   Rscript dev/alapse-readers.R <base|lane>
lib_arg <- commandArgs(TRUE)[1]
source("dev/alapse-common.R")
say <- function(lab, r) {
  if (inherits(r$v, "error")) cat("  ", lab, " ERROR: ", conditionMessage(r$v),
                                  "\n", sep = "")
  for (i in seq_along(r$w)) {
    cat("  ", lab, " WARN [", r$wc[i], "]: ", r$w[i], "\n", sep = "")
  }
  for (i in seq_along(r$m)) {
    cat("  ", lab, " MSG [", r$mc[i], "]: ", r$m[i], "\n", sep = "")
  }
}
set.seed(5)
d <- data.frame(x = stats::rnorm(100), z = stats::rnorm(100),
                g = factor(sample(letters[1:3], 100, TRUE)))
d$y <- 1 + 0.5 * d$x + stats::rnorm(100)
set.seed(9)
dz <- data.frame(x = stats::rnorm(200))
dz$y <- ifelse(stats::runif(200) < plogis(-1.5 + dz$x), 0,
               brms::rasym_laplace(200, 1 + 0.5 * dz$x, 1, 0.5))
cases <- list(
  plain = list(fo = y ~ x + g, fam = asym_laplace(), d = d),
  dist = list(fo = bf(y ~ x, sigma ~ z, quantile ~ z), fam = asym_laplace(),
              d = d),
  zi = list(fo = bf(y ~ x, zi ~ x), fam = zero_inflated_asym_laplace(),
            d = dz))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  cat("==", nm, "\n")
  r <- cap(frm(cs$fo, family = cs$fam, data = cs$d))
  say("frm", r)
  if (inherits(r$v, "error")) next
  f <- r$v
  cat("  opt:", f$opt$convergence, f$opt$message, "\n")
  if (!is.null(f$opt$ald_cert)) str(f$opt$ald_cert)
  s <- cap(summary(f))
  say("summary", s)
  out <- capture.output(print(s$v))
  cat(paste0("  | ", out), sep = "\n")
  v <- cap(vcov(f))
  say("vcov", v)
  print(v$v)
  ci <- cap(confint(f))
  say("confint", ci)
  print(ci$v)
  fe <- cap(fixef(f))
  say("fixef", fe)
  print(fe$v)
  ft <- cap(fitted(f))
  say("fitted", ft)
  if (is.matrix(ft$v)) {
    cat(sprintf("  fitted: %d of %d SEs NaN\n",
                sum(is.nan(ft$v[, "Est.Error"])), nrow(ft$v)))
  }
  if (nm == "zi") {
    pz <- cap(fitted(f, dpar = "zi"))
    say("fitted zi", pz)
    if (is.matrix(pz$v)) {
      cat(sprintf("  fitted zi: %d of %d SEs NaN\n",
                  sum(is.nan(pz$v[, "Est.Error"])), nrow(pz$v)))
    }
  }
  pr <- cap(predict(f, newdata = cs$d[1:3, ], se.fit = TRUE))
  say("predict", pr)
  str(pr$v)
  dg <- cap(diagnose(f))
  say("diagnose", dg)
  out <- capture.output(print(dg$v))
  cat(paste0("  | ", out), sep = "\n")
  if (requireNamespace("emmeans", quietly = TRUE) && nm == "plain") {
    em <- cap(emmeans::emmeans(f, ~ g))
    say("emmeans", em)
    print(em$v)
  }
  if (requireNamespace("marginaleffects", quietly = TRUE)) {
    me <- cap(marginaleffects::avg_slopes(f, variables = "x"))
    say("avg_slopes", me)
    print(me$v)
  }
  pf <- cap(confint(f, method = "profile", parm = "x"))
  say("confint profile", pf)
  print(pf$v)
  bt <- cap(frm_bootstrap(f, nsim = 20, seed = 1))
  say("bootstrap", bt)
  print(bt$v)
}
cat("DONE\n")
