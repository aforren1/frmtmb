# Reviewer, check of fix round 2: does every legitimate refit still use
# the channel (the frozen bases make ns() + NA rows the witness: without
# the channel update(fit) is 1.49e-2 off the fit), and the data_vars
# rule. Arg: lane. Log: dev/drop1raw-rev-out/r2-lane.txt (seed 5, n = 200)
.libPaths(c("C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
ll <- function(f) as.numeric(logLik(f))
set.seed(5)
n <- 200
d <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1, x3 = rnorm(n),
                g = factor(sample(letters[1:8], n, TRUE)))
d$z <- 1 + d$x + log(d$x2) + rnorm(8)[d$g] + rnorm(n)
d$cnt <- rpois(n, exp(0.3 + 0.3 * d$x))
d$x3[c(3, 17, 40)] <- NA
d$w <- rnorm(n)
chk <- function(lab, fit, u) {
  if (inherits(u, "error")) {
    cat(sprintf("%-52s ERROR %s\n", lab, substr(conditionMessage(u), 1, 120)))
    return(invisible())
  }
  cat(sprintf("%-52s fit %.10f | refit %.10f | diff %.2e | na_action %s\n",
              lab, ll(fit), ll(u), ll(u) - ll(fit),
              paste(unclass(u$frame$na_action), collapse = ",")))
}
tu <- function(e) tryCatch(q(e), error = function(err) err)

cat("\n== 1. spellings: the channel must be used (diff 0) ==\n")
f1 <- q(frm(bf(z ~ splines::ns(x2, 3) + x3), data = d))
chk("bf(), update(fit)", f1, tu(update(f1)))
f2 <- tu(frm(z ~ splines::ns(x2, 3) + x3, data = d))
if (!inherits(f2, "error")) chk("plain formula, update(fit)", f2, tu(update(f2)))
f3 <- tu(frm(bf(z ~ splines::ns(x2, 3) + x3), family = "gaussian", data = d))
if (!inherits(f3, "error")) chk("family = \"gaussian\", update(fit)", f3,
                                tu(update(f3)))
f4 <- tu(frm(bf(z ~ splines::ns(x2, 3) + x3), family = gaussian, data = d))
if (!inherits(f4, "error")) chk("family = gaussian (function), update(fit)",
                                f4, tu(update(f4)))
mk <- function(dat, fam) {
  ff <- z ~ splines::ns(x2, 3) + x3
  frm(bf(ff), family = fam, data = dat)
}
f5 <- q(mk(d, gaussian()))
chk("formula from a function frame, update(fit)", f5, tu(update(f5)))
chk("  update(fit, family = student()) vs frm", q(frm(bf(z ~ splines::ns(x2,
    3) + x3), family = student(), data = d)), tu(update(f5,
    family = student())))
fam_s <- student()
chk("  update(fit, family = fam_s) (a symbol)", q(frm(bf(z ~ splines::ns(x2,
    3) + x3), family = student(), data = d)), tu(update(f5, family = fam_s)))
chk("  update(fit, bf(z ~ splines::ns(x2, 3) + x3)) complete",
    f1, tu(update(f5, bf(z ~ splines::ns(x2, 3) + x3))))
chk("  update(fit, z ~ splines::ns(x2, 3) + x3) complete",
    f1, tu(update(f5, z ~ splines::ns(x2, 3) + x3)))
chk("  update(fit, . ~ .)", f1, tu(update(f5, . ~ .)))
chk("  update(fit, bf(. ~ ., sigma ~ 1))", f1,
    tu(update(f5, bf(. ~ ., sigma ~ 1))))
# drop1 AIC(-x3) must equal AIC(update(~ . - x3)), as in round 1
dt <- q(drop1(f5)); u3 <- q(update(f5, ~ . - x3))
cat(sprintf("drop1 AIC(-x3) %.8f | AIC(update(~ . - x3)) %.8f\n",
            dt["x3", "AIC"], AIC(u3)))
# step via evaluate = FALSE
s <- q(step(q(frm(bf(z ~ splines::ns(x2, 3) + x3 + w), data = d)),
            trace = 0))
cat("step:", deparse1(formula(s)), sprintf("logLik %.10f nobs %d", ll(s),
                                           nobs(s)), "\n")
fp <- q(frm(bf(cnt ~ splines::ns(x2, 3) + x3), family = poisson(), data = d))
chk("poisson(), update(fit)", fp, tu(update(fp)))
fpm <- tu(q(frm(bf(cnt ~ splines::ns(x2, 3) + x3), family = poisson, data = d)))
if (!inherits(fpm, "error")) chk("poisson (function), update(fit)", fpm,
                                 tu(update(fpm)))

cat("\n== 2. data_vars ==\n")
x9 <- rnorm(n)
fe <- q(frm(bf(z ~ x + log(x2) + x9), data = d))
cat("data_vars of a fit reading x9 from the environment:",
    paste(fe$frame$data_vars, collapse = ","), "\n")
chk("update(fe, ~ . + w), x9 not in d", fe, tu(update(fe, ~ . + w)))
in_fun <- function(fit, dloc) { d <- dloc; tu(update(fit, ~ . + w)) }
loc3 <- d[, c("z", "x", "x3", "g", "w")]
chk("local d without x2", fe, in_fun(fe, loc3))
x2 <- rexp(n) + 5
chk("local d without x2, global x2", fe, in_fun(fe, loc3))
rm(x2)
dl <- as.list(d)
fl <- q(frm(bf(z ~ x + log(x2) + x3), data = dl))
chk("named list data: update(fl)", fl, tu(update(fl)))
chk("named list data: update(fl, ~ . + w)", fl, tu(update(fl, ~ . + w)))
A <- diag(8); dimnames(A) <- list(letters[1:8], letters[1:8])
fd <- tu(q(frm(bf(z ~ x + log(x2) + (1 | gr(g, cov = A))), data = d,
               data2 = list(A = A))))
if (inherits(fd, "error")) cat("data2 fit ERROR", conditionMessage(fd), "\n") else {
  chk("data2 fit: update(fd)", fd, tu(update(fd)))
  chk("data2 fit: update(fd, ~ . + w)", fd, tu(update(fd, ~ . + w)))
  cat("data2 fit data_vars:", paste(fd$frame$data_vars, collapse = ","), "\n")
}
