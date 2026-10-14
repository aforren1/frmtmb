# Probe brms 2.23.0's wiener family without compiling: links accepted,
# Stan code under scalar and distributional spellings, default priors.
# Output: dev/famuni-out/wiener-brms.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
cat("rstan", format(packageVersion("rstan")), " stanc",
    rstan::stan_version(), "\n")
set.seed(20261013)
n <- 80
d <- data.frame(x = rnorm(n))
sim <- RWiener::rwiener(n, alpha = 1.4, tau = 0.25, beta = 0.5, delta = 0.8)
d$rt <- sim$q
d$d <- as.integer(sim$resp == "upper")
cat("min rt", min(d$rt), "\n")

cat("\n== formals(brms::wiener)\n"); print(formals(brms::wiener))
cat("\n== .family_wiener\n"); str(brms:::.family_wiener())

cat("\n== links accepted by wiener()\n")
lks <- c("identity", "log", "softplus", "squareplus", "logit", "inverse",
         "probit", "cloglog", "sqrt", "tan_half", "logm1")
for (arg in c("link", "link_bs", "link_ndt", "link_bias")) {
  for (lk in lks) {
    r <- tryCatch({do.call(brms::wiener, setNames(list(lk), arg)); "ok"},
                  error = function(e) paste("ERR", conditionMessage(e)))
    cat(sprintf("%-9s %-10s %s\n", arg, lk, substr(r, 1, 70)))
  }
}

show <- function(label, form, fam = wiener()) {
  cat("\n\n================", label, "\n")
  sc <- stancode(form, family = fam, data = d)
  cat(sc, "\n")
  cat("\n-- default_prior\n")
  print(default_prior(form, family = fam, data = d))
}
show("scalar dpars", bf(rt | dec(d) ~ x))
show("bs ~ x, ndt ~ 1, bias ~ 1",
     bf(rt | dec(d) ~ x, bs ~ x, ndt ~ 1, bias ~ 1))
show("intercept-only mu, bs ~ 1, ndt ~ 1, bias ~ 1",
     bf(rt | dec(d) ~ 1, bs ~ 1, ndt ~ 1, bias ~ 1))
show("bias = 0.5 fixed, ndt ~ x",
     bf(rt | dec(d) ~ x, ndt ~ x, bias = 0.5))
show("mu link softplus, link_ndt identity",
     bf(rt | dec(d) ~ x), wiener(link = "softplus", link_ndt = "identity"))

cat("\n\n== addition terms brms refuses for wiener\n")
d$cen <- "none"; d$w <- 1
for (adt in c("cens(cen)", "trunc(ub = 5)", "weights(w)", "mi()")) {
  f <- as.formula(paste("rt | dec(d) +", adt, "~ x"))
  r <- tryCatch({stancode(f, family = wiener(), data = d); "ok"},
                error = function(e) paste("ERR", conditionMessage(e)))
  cat(adt, ":", substr(r, 1, 160), "\n")
}
cat("\n== mixture(wiener, wiener)\n")
r <- tryCatch({
  sc <- stancode(bf(rt | dec(d) ~ 1), family = mixture(wiener, wiener),
                 data = d)
  grep("wiener", strsplit(sc, "\n")[[1]], value = TRUE)
}, error = function(e) paste("ERR", conditionMessage(e)))
cat(r, sep = "\n")
cat("\n== missing dec()\n")
r <- tryCatch({stancode(bf(rt ~ x), family = wiener(), data = d); "ok"},
              error = function(e) paste("ERR", conditionMessage(e)))
cat(r, "\n")
cat("\n== dec() given a factor\n")
d$fd <- factor(ifelse(d$d == 1, "upper", "lower"))
r <- tryCatch({sd <- standata(bf(rt | dec(fd) ~ x), family = wiener(),
                              data = d)
  paste(identical(as.integer(sd$dec), d$d), "head", paste(head(sd$dec), collapse = ","))},
  error = function(e) paste("ERR", conditionMessage(e)))
cat(r, "\n")
cat("\n== dec() given 'lower'/'upper' character\n")
d$cd <- ifelse(d$d == 1, "upper", "lower")
r <- tryCatch({sd <- standata(bf(rt | dec(cd) ~ x), family = wiener(),
                              data = d)
  paste(identical(as.integer(sd$dec), d$d))},
  error = function(e) paste("ERR", conditionMessage(e)))
cat(r, "\n")
