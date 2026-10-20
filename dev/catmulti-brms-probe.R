# What brms 2.23.0 does for categorical() / multinomial() refcat, the bare
# multinomial(), and the simplex families in mvbf(): stancode, standata,
# default_prior, family objects. No compile.
# Output: dev/catmulti-out/brms-probe.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(1)
n <- 60
dd <- data.frame(x = rnorm(n), z = rnorm(n), g = gl(6, 10))
dd$f <- factor(sample(c("a", "b", "c"), n, TRUE))
dd$f2 <- factor(sample(c("p.q", "r_s", "t"), n, TRUE))
P <- matrix(rgamma(n * 3, 2), n, 3); P <- P / rowSums(P)
Yc <- t(sapply(seq_len(n), function(i) rmultinom(1, 10, P[i, ])))
colnames(Yc) <- c("y1", "y2", "y3")
dd$y1 <- Yc[, 1]; dd$y2 <- Yc[, 2]; dd$y3 <- Yc[, 3]
dd$N <- 10L
dd$Yc <- Yc
Yu <- Yc; colnames(Yu) <- NULL
dd$Yu <- Yu
G <- matrix(rgamma(n * 3, 2), n, 3); Y <- G / rowSums(G)
colnames(Y) <- c("a", "b", "c")
dd$Y <- Y
sec <- function(t) cat("\n========", t, "========\n")
try_show <- function(expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)),
                warning = function(w) paste("WARNING:", conditionMessage(w)))
  print(r)
}
dp_names <- function(form, fam) {
  bt <- brms:::validate_formula(form, data = dd, family = fam)
  print(bt$family$dpars)
  print(bt$family$refcat)
  print(bt$family$cats)
}

sec("family objects")
print(args(categorical)); print(args(multinomial))
str(categorical(refcat = "b"))
str(categorical(refcat = NA))
str(multinomial(refcat = NA))
try_show(categorical(refcat = c("a", "b")))
try_show(categorical(refcat = 2))
str(categorical(refcat = 2))

sec("categorical dpars by refcat")
dp_names(bf(f ~ x), categorical())
dp_names(bf(f ~ x), categorical(refcat = "b"))
dp_names(bf(f ~ x), categorical(refcat = NA))
try_show(dp_names(bf(f ~ x), categorical(refcat = "zz")))
dp_names(bf(f2 ~ x), categorical(refcat = "r_s"))

sec("categorical refcat = b stancode")
cat(stancode(bf(f ~ x), data = dd, family = categorical(refcat = "b")))
sec("categorical refcat = b standata Y")
sd <- standata(bf(f ~ x), data = dd, family = categorical(refcat = "b"))
print(names(sd)); print(head(sd$Y)); print(head(dd$f)); print(sd$ncat)
sec("categorical refcat = NA stancode")
cat(stancode(bf(f ~ x), data = dd, family = categorical(refcat = NA)))
sec("categorical default priors")
print(default_prior(bf(f ~ x), data = dd, family = categorical()))
print(default_prior(bf(f ~ x), data = dd, family = categorical(refcat = "b")))
print(default_prior(bf(f ~ x), data = dd, family = categorical(refcat = NA)))
sec("categorical refcat=b with own formula for mua")
try_show(cat(stancode(bf(f ~ x, mua ~ z), data = dd,
                      family = categorical(refcat = "b"))))
try_show(cat(stancode(bf(f ~ x, mub ~ z), data = dd,
                      family = categorical(refcat = "b"))))
sec("categorical refcat=b with (1 | ID | g)")
print(default_prior(bf(f ~ x + (1 | ID | g)), data = dd,
                    family = categorical(refcat = "b")))

sec("multinomial bare, cbind")
dp_names(bf(cbind(y1, y2, y3) | trials(N) ~ x), multinomial())
dp_names(bf(Yc | trials(N) ~ x), multinomial())
dp_names(bf(Yu | trials(N) ~ x), multinomial())
dp_names(bf(Yc | trials(N) ~ x), multinomial(refcat = "y2"))
dp_names(bf(Yc | trials(N) ~ x), multinomial(refcat = NA))
try_show(dp_names(bf(Yc | trials(N) ~ x), multinomial(refcat = "zz")))
sec("multinomial stancode refcat = y2")
cat(stancode(bf(Yc | trials(N) ~ x), data = dd,
             family = multinomial(refcat = "y2")))
sec("multinomial without trials")
try_show(cat(stancode(bf(Yc ~ x), data = dd, family = multinomial())))
sd <- try(standata(bf(Yc ~ x), data = dd, family = multinomial()))
if (!inherits(sd, "try-error")) {print(names(sd)); print(head(sd$trials))}
sec("multinomial default priors")
print(default_prior(bf(Yc | trials(N) ~ x), data = dd, family = multinomial()))
print(default_prior(bf(Yc | trials(N) ~ x), data = dd,
                    family = multinomial(refcat = "y2")))
print(default_prior(bf(Yc | trials(N) ~ x), data = dd,
                    family = multinomial(refcat = NA)))
sec("multinomial one column / non-integer")
dd$Y1 <- matrix(dd$y1, ncol = 1)
try_show(standata(bf(Y1 | trials(N) ~ x), data = dd, family = multinomial()))
dd$Yh <- dd$Yc + 0.5
try_show(standata(bf(Yh | trials(N) ~ x), data = dd, family = multinomial()))
dd$Ng <- 11L
try_show(standata(bf(Yc | trials(Ng) ~ x), data = dd, family = multinomial()))

sec("dirichlet refcat NA stancode")
cat(stancode(bf(Y ~ x), data = dd, family = dirichlet(refcat = NA)))
sec("logistic_normal refcat NA")
try_show(logistic_normal(refcat = NA))

sec("mvbf: dirichlet + gaussian")
try_show(cat(stancode(mvbf(bf(Y ~ x), bf(z ~ x, family = gaussian())),
                      data = dd, family = dirichlet())))
mvf1 <- bf(Y ~ x, family = dirichlet()) + bf(z ~ x, family = gaussian()) +
  set_rescor(FALSE)
try_show(cat(stancode(mvf1, data = dd)))
try_show(print(default_prior(mvf1, data = dd)))
sec("mvbf: logistic_normal + gaussian")
mvf2 <- bf(Y ~ x, family = logistic_normal()) +
  bf(z ~ x, family = gaussian()) + set_rescor(FALSE)
try_show(cat(stancode(mvf2, data = dd)))
try_show(print(default_prior(mvf2, data = dd)))
sec("mvbf: logistic_normal + logistic_normal")
dd$W <- Y[, c(2, 3, 1)]; colnames(dd$W) <- c("u", "v", "w")
mvf3 <- bf(Y ~ x, family = logistic_normal()) +
  bf(W ~ x, family = logistic_normal()) + set_rescor(FALSE)
try_show(cat(stancode(mvf3, data = dd)))
try_show(print(default_prior(mvf3, data = dd)))
sec("mvbf: dirichlet + dirichlet")
mvf4 <- bf(Y ~ x, family = dirichlet()) +
  bf(W ~ x, family = dirichlet(refcat = "w")) + set_rescor(FALSE)
try_show(cat(stancode(mvf4, data = dd)))
sec("mvbf: categorical + gaussian")
mvf5 <- bf(f ~ x, family = categorical(refcat = "c")) +
  bf(z ~ x, family = gaussian()) + set_rescor(FALSE)
try_show(print(default_prior(mvf5, data = dd)))
