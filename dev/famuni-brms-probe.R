# Probe brms 2.23.0 for the three families this lane adds: links, dpars,
# default priors, Stan code with distributional formulas, cens/trunc and
# mixture support. Output: dev/famuni-out/brms-probe.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(1)
n <- 60
d <- data.frame(x = rnorm(n), g = gl(6, 10))
d$yc <- 2 + rexp(n)
d$yd <- rpois(n, 3)
d$cen <- sample(c("none", "right", "left"), n, TRUE)
d$lb <- 0.5
for (fam in c("frechet", "gen_extreme_value", "discrete_weibull")) {
  cat("\n\n================ ", fam, "\n")
  f <- get(fam, asNamespace("brms"))
  cat("formals:\n"); print(formals(f))
  fo <- f()
  str(unclass(fo)[c("family", "link", "dpars", "type", "ybounds",
                    "closed", "specials")])
  fi <- get(paste0(".family_", fam), asNamespace("brms"))()
  cat("links allowed:", fi$links, "\n")
  cat("specials:", fi$specials, "\n")
  cat("ad:", fi$ad, "\n")
  y <- if (fam == "discrete_weibull") "yd" else "yc"
  others <- setdiff(fo$dpars, "mu")
  form <- as.formula(paste(y, "~ x + (1 | g)"))
  dl <- lapply(others, function(dp) as.formula(paste(dp, "~ x")))
  bfo <- do.call(bf, c(list(form), dl))
  cat("\n-- default_prior, distributional\n")
  print(default_prior(bfo, data = d, family = fo))
  cat("\n-- default_prior, intercept-only dpars\n")
  print(default_prior(bf(form), data = d, family = fo))
  cat("\n-- stancode, distributional\n")
  cat(stancode(bfo, data = d, family = fo), "\n")
  cat("\n-- stancode, intercept-only dpars\n")
  cat(stancode(bf(form), data = d, family = fo), "\n")
  for (adt in c("cens(cen)", "trunc(lb = 0.5)", "weights(x^2)")) {
    yy <- paste0(y, " | ", adt, " ~ x")
    cat("\n-- addition", adt, ":\n")
    r <- tryCatch(stancode(bf(as.formula(yy)), data = d, family = fo),
                  error = function(e) paste("ERROR:", conditionMessage(e)))
    cat(grep("lpdf|lpmf|lcdf|lccdf|ERROR|target", strsplit(r, "\n")[[1]],
             value = TRUE), sep = "\n")
  }
  cat("\n-- mixture with itself:\n")
  r <- tryCatch({
    sc <- stancode(bf(form), data = d, family = mixture(fo, fo))
    grep("lpdf|lpmf", strsplit(sc, "\n")[[1]], value = TRUE)
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(r, sep = "\n")
  cat("\n-- links for mu tried:\n")
  for (lk in c("log", "identity", "logit", "softplus", "inverse", "probit",
               "cloglog", "sqrt")) {
    r <- tryCatch({f(link = lk); "ok"},
                  error = function(e) conditionMessage(e))
    cat(lk, ":", r, "\n")
  }
}
cat("\n-- wiener\n")
fi <- brms:::.family_wiener()
str(fi)
print(formals(brms::wiener))
