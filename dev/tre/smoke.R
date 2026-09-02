# Quick end-to-end smoke test of the t-random-effect implementation.
# Run: Rscript dev/tre/smoke.R
sink("dev/tre/smoke.txt")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
set.seed(1)
G <- 40; n <- 8
d <- data.frame(x = rnorm(G * n), g = factor(rep(1:G, each = n)))
b <- rt(G, df = 5)
d$y <- 1 + 0.5 * d$x + b[d$g] + rnorm(G * n)

cat("== scalar t intercept ==\n")
f <- frm(bf(y ~ x + (1 | gr(g, dist = "student"))) + gaussian(), data = d)
print(summary(f))
cat("\ncovstruct: ", f$frame$re_blocks[[1]]$covstruct,
    "  dist_nu: ", f$frame$re_blocks[[1]]$dist_nu, "\n")
print(VarCorr(f))
cat("logLik ", as.numeric(logLik(f)), " AIC ", AIC(f), "\n")

cat("\n== dist_nu = 3 ==\n")
f3 <- frm(bf(y ~ x + (1 | gr(g, dist = "student", dist_nu = 3))) +
            gaussian(), data = d)
print(VarCorr(f3))

cat("\n== gaussian limit ==\n")
fg <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d)
fbig <- frm(bf(y ~ x + (1 | gr(g, dist = "student",
                               dist_nu = 1e8))) + gaussian(), data = d)
cat("gaussian theta ", f$estimates$theta, "\n")
cat("gauss  : ", format(c(fg$estimates$beta, fg$estimates$theta),
                        digits = 12), "\n")
cat("nu=1e8 : ", format(c(fbig$estimates$beta, fbig$estimates$theta),
                        digits = 12), "\n")
cat("max abs diff: ",
    max(abs(c(fg$estimates$beta, fg$estimates$theta) -
              c(fbig$estimates$beta, fbig$estimates$theta))), "\n")
cat("logLik diff : ",
    abs(as.numeric(logLik(fg)) - as.numeric(logLik(fbig))), "\n")

cat("\n== correlated slopes ==\n")
fc <- frm(bf(y ~ x + (x | gr(g, dist = "student"))) + gaussian(),
          data = d)
print(VarCorr(fc))
cat("covstruct ", fc$frame$re_blocks[[1]]$covstruct, "\n")
fcg <- frm(bf(y ~ x + (x | gr(g, dist = "student", dist_nu = 1e8))) +
             gaussian(), data = d)
fcn <- frm(bf(y ~ x + (x | g)) + gaussian(), data = d)
cat("corr-slope gaussian limit max diff: ",
    max(abs(c(fcg$estimates$beta, fcg$estimates$theta) -
              c(fcn$estimates$beta, fcn$estimates$theta))), "\n")

cat("\n== diag ==\n")
fd <- frm(bf(y ~ x + diag(x | gr(g, dist = "student"))) + gaussian(),
          data = d)
cat("covstruct ", fd$frame$re_blocks[[1]]$covstruct, "\n")
print(VarCorr(fd))
fdg <- frm(bf(y ~ x + diag(x | gr(g, dist = "student",
                                  dist_nu = 1e8))) + gaussian(), data = d)
fdn <- frm(bf(y ~ x + diag(x | g)) + gaussian(), data = d)
cat("diag gaussian limit max diff: ",
    max(abs(c(fdg$estimates$beta, fdg$estimates$theta) -
              c(fdn$estimates$beta, fdn$estimates$theta))), "\n")

cat("\n== downstream ==\n")
print(head(ranef(f)[[1]], 3))
print(head(fitted(f), 3))
print(head(predict(f), 3))
print(confint(f, parm = "theta_1"))
cat("simulate: ")
print(dim(simulate(f, nsim = 2)))
nd <- data.frame(x = 0, g = factor("newlevel", levels = "newlevel"))
print(predict(f, newdata = nd, allow_new_levels = TRUE, se.fit = TRUE))

cat("\n== quadrature ==\n")
fq <- frm(bf(y ~ x + (1 | gr(g, dist = "student"))) + gaussian(),
          data = d, quadrature = TRUE)
cat("Laplace  : ", format(c(f$estimates$beta, f$estimates$theta),
                          digits = 10), " ll ",
    as.numeric(logLik(f)), "\n")
cat("quadrature: ", format(c(fq$estimates$beta, fq$estimates$theta),
                           digits = 10), " ll ",
    as.numeric(logLik(fq)), "\n")

cat("\n== REML ==\n")
fr <- frm(bf(y ~ x + (1 | gr(g, dist = "student"))) + gaussian(),
          data = d, REML = TRUE)
cat("ML theta   ", f$estimates$theta, "\n")
cat("REML theta ", fr$estimates$theta, "\n")

cat("\n== guards ==\n")
gg <- function(e) tryCatch({ eval(e); "NO ERROR" },
                           error = function(x) conditionMessage(x))
A <- diag(G); dimnames(A) <- list(levels(d$g), levels(d$g))
cat("1 ", gg(quote(frm(bf(y ~ (1 | gr(g, dist = "poisson"))) +
                         gaussian(), data = d))), "\n")
cat("2 ", gg(quote(frm(bf(y ~ (1 | gr(g, dist = "student",
                                      dist_nu = 1))) + gaussian(),
                        data = d))), "\n")
cat("3 ", gg(quote(frm(bf(y ~ (1 | gr(g, cov = A, dist = "student"))) +
                         gaussian(), data = d, data2 = list(A = A)))), "\n")
cat("4 ", gg(quote(frm(bf(y ~ (1 | gr(g, dist_nu = 4))) + gaussian(),
                        data = d))), "\n")
cat("5 ", gg(quote(frm(bf(y ~ cs(x | gr(g, dist = "student"))) +
                         gaussian(), data = d))), "\n")
cat("6 ", gg(quote(frm(bf(y ~ (1 | q | gr(g, dist = "student")) +
                            (0 + x | q | gr(g, dist = "student"))) +
                         gaussian(), data = d))), "\n")
cat("7 ", gg(quote(frm(bf(y ~ (1 | gr(g, dist = "gaussian"))) +
                         gaussian(), data = d))), "\n")
cat("8 ", gg(quote(frm(bf(y ~ us_t(x | g)) + gaussian(), data = d))),
    "\n")
sink()
cat("done\n")
