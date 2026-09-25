.libPaths(c("/opt/rlib/deps", "/opt/r/lib/R/library"))
suppressMessages(library(brms))
set.seed(1)
n <- 40
d <- data.frame(y = rnorm(n), x = rnorm(n), z = rnorm(n),
                f = factor(rep(c("a", "b", "c", "d"), 10)),
                g = factor(rep(1:5, 8)),
                yo = factor(sample(1:4, n, TRUE), ordered = TRUE),
                xm = factor(sample(1:3, n, TRUE), ordered = TRUE),
                y2 = rnorm(n))
try1 <- function(lab, expr) {
  cat("\n#####", lab, "\n")
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("WARNING:", conditionMessage(w), "\n"); invokeRestart("muffleWarning")
  }), error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
  invisible(NULL)
}
sd_X <- function(f, data = d, fam = gaussian(), nm = "X") {
  s <- make_standata(f, data = data, family = fam)
  X <- s[[nm]]
  print(dim(X)); print(colnames(X)); print(head(X, 4))
  invisible(s)
}
pr <- function(f, data = d, fam = gaussian()) {
  p <- get_prior(f, data = data, family = fam)
  print(as.data.frame(p)[, c("prior", "class", "coef", "dpar", "nlpar",
                             "resp")])
}

try1("0 + Intercept + x", { sd_X(y ~ 0 + Intercept + x); pr(y ~ 0 + Intercept + x) })
try1("1 + x", { sd_X(y ~ 1 + x); pr(y ~ 1 + x) })
try1("0 + x + Intercept (order)", sd_X(y ~ 0 + x + Intercept))
try1("0 + Intercept + f", sd_X(y ~ 0 + Intercept + f))
try1("0 + f", sd_X(y ~ 0 + f))
try1("0 + Intercept + x:f", sd_X(y ~ 0 + Intercept + x:f))
try1("0 + Intercept + x*f", sd_X(y ~ 0 + Intercept + x * f))
try1("0 + Intercept + f + z:f", sd_X(y ~ 0 + Intercept + f + x:f))
try1("0 + Intercept:x", sd_X(y ~ 0 + Intercept:x))
try1("0 + Intercept alone", { sd_X(y ~ 0 + Intercept); pr(y ~ 0 + Intercept) })
try1("lowercase intercept", sd_X(y ~ 0 + intercept + x))
try1("1 + Intercept + x", sd_X(y ~ 1 + Intercept + x))
try1("Intercept + x", sd_X(y ~ Intercept + x))
try1("-1 + Intercept + x", sd_X(y ~ -1 + Intercept + x))
try1("x + Intercept - 1", sd_X(y ~ x + Intercept - 1))
dI <- d; dI$Intercept <- 2
try1("data col Intercept=2, 0+Intercept", sd_X(y ~ 0 + Intercept + x, data = dI))
try1("data col Intercept=2, 1 + x (unused)", sd_X(y ~ 1 + x, data = dI))
try1("data col Intercept=2, x + Intercept", sd_X(y ~ x + Intercept, data = dI))
dI1 <- d; dI1$Intercept <- 1
try1("data col Intercept=1, 0+Intercept", sd_X(y ~ 0 + Intercept + x, data = dI1))
try1("data col Intercept=2, 1 + x, sigma ~ 0 + Intercept",
     sd_X(bf(y ~ 1 + x, sigma ~ 0 + Intercept + z), data = dI))
try1("data col Intercept=2, 0 + Intercept + x (x also named)",
     sd_X(bf(y ~ 0 + Intercept + x), data = dI))
try1("dpar sigma ~ 0 + Intercept + z", {
  sd_X(bf(y ~ x, sigma ~ 0 + Intercept + z), nm = "X_sigma")
  pr(bf(y ~ x, sigma ~ 0 + Intercept + z)) })
try1("mu 1+x, sigma ~ 1 + Intercept + z", sd_X(bf(y ~ 1 + x, sigma ~ 1 + Intercept + z), nm = "X_sigma"))
try1("mu 0+Intercept, sigma ~ 1 + Intercept + z", sd_X(bf(y ~ 0 + Intercept + x, sigma ~ 1 + Intercept + z), nm = "X_sigma"))
try1("nl a ~ 0 + Intercept + x", {
  f <- bf(y ~ a * exp(b * z), a ~ 0 + Intercept + x, b ~ 1, nl = TRUE)
  sd_X(f, nm = "X_a"); pr(f) })
try1("nl a ~ 1 + x", {
  f <- bf(y ~ a * exp(b * z), a ~ 1 + x, b ~ 1, nl = TRUE)
  sd_X(f, nm = "X_a"); pr(f) })
try1("center = FALSE", { sd_X(bf(y ~ x, center = FALSE)); pr(bf(y ~ x, center = FALSE)) })
try1("center = FALSE, sigma", { pr(bf(y ~ x, sigma ~ z, center = FALSE)) })
try1("center = FALSE, 0 + x", { pr(bf(y ~ 0 + x, center = FALSE)) })
try1("center = FALSE factor", sd_X(bf(y ~ f, center = FALSE)))
try1("ordinal 0 + Intercept + x", {
  sd_X(yo ~ 0 + Intercept + x, fam = cumulative()); pr(yo ~ 0 + Intercept + x, fam = cumulative()) })
try1("ordinal stancode", {
  cat(grep("Intercept|b\\b|Xc|means_X", strsplit(make_stancode(yo ~ 0 + Intercept + x, data = d, family = cumulative()), "\n")[[1]], value = TRUE), sep = "\n") })
try1("ordinal center=FALSE", pr(bf(yo ~ x, center = FALSE), fam = cumulative()))
try1("cs(): 0 + Intercept + cs(x) acat", { sd_X(yo ~ 0 + Intercept + cs(x), fam = acat(), nm = "Xcs"); pr(yo ~ 0 + Intercept + cs(x), fam = acat()) })
try1("cs(Intercept)", { s <- make_standata(yo ~ 0 + Intercept + cs(Intercept) + x, data = d, family = acat()); print(names(s)); print(head(s$Xcs)) ; pr(yo ~ 0 + cs(Intercept) + x, fam = acat())})
try1("mo()", { sd_X(y ~ 0 + Intercept + mo(xm)); pr(y ~ 0 + Intercept + mo(xm)) })
try1("smooth", { s <- sd_X(y ~ 0 + Intercept + s(x)); print(names(s)); pr(y ~ 0 + Intercept + s(x)) })
try1("smooth 1 + s(x)", { s <- sd_X(y ~ 1 + s(x)); print(names(s)) })
try1("RE (0 + Intercept | g)", { s <- make_standata(y ~ 0 + Intercept + x + (0 + Intercept | g), data = d); print(s$Z_1_1[1:4]); pr(y ~ 0 + Intercept + x + (0 + Intercept | g)) })
try1("RE (1 + Intercept | g) with 1+x", { s <- make_standata(y ~ 1 + x + (0 + Intercept | g), data = d); print(names(s)) })
try1("mv", { f <- mvbf(bf(y ~ 0 + Intercept + x), bf(y2 ~ x), rescor = FALSE); pr(f) })
try1("Intercept in mu and a data col not 1 but mu has 1",
     sd_X(bf(y ~ 1 + x), data = dI))
try1("0 + Intercept + Intercept:x", sd_X(y ~ 0 + Intercept + Intercept:x))
try1("0 + Intercept + I(Intercept*2)", sd_X(y ~ 0 + Intercept + I(2*Intercept)))
try1("0 + Intercept + poly(x,2)", sd_X(y ~ 0 + Intercept + poly(x, 2)))
try1("zi dpar", pr(bf(y ~ x, hu ~ 0 + Intercept + z), fam = hurdle_lognormal()))
