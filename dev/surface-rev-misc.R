# Reviewer of lane surface: claims 6 (categorical message), 8 (fitted()
# with no fixed disc column), 2 (plot() of a fit refusals).
#   Rscript dev/surface-rev-misc.R lane|base
source("dev/surface-rev-env.R")
arm <- rev_env(commandArgs(TRUE)[1])
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, format(packageVersion("frmtmb")), "\n")
grDevices::pdf(NULL)

## 6: method = "predict" on categorical-like families
set.seed(3)
n <- 200
d <- data.frame(x = rnorm(n))
eta2 <- 0.5 + d$x; eta3 <- -0.3 - d$x
pr <- cbind(1, exp(eta2), exp(eta3)); pr <- pr / rowSums(pr)
d$y <- factor(apply(pr, 1, function(p) sample(c("a", "b", "c"), 1,
                                              prob = p)))
fc <- frm(bf(y ~ x), family = categorical(), data = d)
cat("categorical dpars:", names(fc$spec$responses[[1]]$dpars), "\n")
rev_show("6 categorical CE predict",
         conditional_effects(fc, "x", method = "predict"))
dd <- d
dd$Y <- t(rmultinom(n, 10, c(0.2, 0.3, 0.5)))[, 1:3]
dd$Y <- t(sapply(seq_len(n), function(i) rmultinom(1, 10, pr[i, ])))
colnames(dd$Y) <- c("a", "b", "c")
dd$size <- 10
fmn <- rev_show("fit multinomial",
                frm(bf(Y | trials(size) ~ x), family = multinomial(),
                    data = dd))
if (!inherits(fmn, "rev_err")) {
  rev_show("6 multinomial CE predict",
           conditional_effects(fmn, "x", method = "predict"))
}
dd$P <- pr * 0.98 + 0.01 / 1.5
dd$P <- dd$P / rowSums(dd$P)
fdi <- rev_show("fit dirichlet",
                frm(bf(P ~ x), family = dirichlet(), data = dd))
if (!inherits(fdi, "rev_err")) {
  rev_show("6 dirichlet CE predict",
           conditional_effects(fdi, "x", method = "predict"))
}
set.seed(5)
do <- data.frame(x = rnorm(n))
do$y <- cut(do$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = FALSE)
fo <- frm(bf(y ~ x), family = cumulative(), data = do)
rev_show("6 cumulative CE predict (control)",
         conditional_effects(fo, "x", method = "predict"))

## 8: fitted() with no fixed disc column
fo2 <- rev_show("8 fit cumulative disc ~ 0 + (1 | g)", {
  do$g <- factor(sample(letters[1:6], n, TRUE))
  frm(bf(y ~ x, disc ~ 0 + (1 | g)), family = cumulative(), data = do)
})
if (!inherits(fo2, "rev_err")) {
  f8 <- rev_show("8 fitted()", fitted(fo2))
  rev_show("8 predict()", predict(fo2))
  rev_show("8 fitted(newdata)", fitted(fo2, newdata = do[1:5, ]))
  rev_show("8 simulate newdata", simulate(fo2, newdata = do[1:5, ]))
  rev_show("8 conditional_effects", conditional_effects(fo2, "x"))
  if (is.array(f8)) cat("   fitted dim", dim(f8), "\n")
}
fo3 <- rev_show("8 fit cumulative disc ~ 1 + (1 | g) (control)",
  frm(bf(y ~ x, disc ~ 1 + (1 | g)), family = cumulative(), data = do))
if (!inherits(fo3, "rev_err")) rev_show("8 fitted() control", fitted(fo3))

## 2: plot() of a fit, every brms argument
set.seed(1)
d1 <- data.frame(x = rnorm(40)); d1$y <- d1$x + rnorm(40)
f1 <- frm(bf(y ~ x), family = gaussian(), data = d1)
for (a in c("pars", "combo", "nvariables", "N", "variable", "regex",
            "fixed", "bins", "theme", "plot", "newpage")) {
  args <- list(f1, 1)
  names(args) <- c("x", a)
  rev_show(paste0("2 plot(fit, ", a, " = 1)"), do.call(plot, args))
}
rev_show("2 plot(fit, col = 2) [graphical par]", plot(f1, col = 2))
rev_show("2 plot(fit, ask = FALSE) control", plot(f1, ask = FALSE))
