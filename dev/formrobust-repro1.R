# Item 1 probe: addition terms given expressions. Before arm: set
# FORMROBUST_LIB="" to run against rellib-r4.
LIB <- Sys.getenv("FORMROBUST_LIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
set.seed(1)
n <- 60
d <- data.frame(x = rnorm(n), wt = runif(n, 0.5, 2), time = runif(n, 1, 3),
                n = rpois(n, 5) + 2L, s = runif(n, 0.2, 0.6),
                g = gl(6, 10))
d$yc <- rpois(n, exp(0.2 + 0.3 * d$x) * 2 * d$time)
d$yb <- rbinom(n, d$n + 1, plogis(0.3 * d$x))
d$y <- rnorm(n, 0.5 * d$x, 1)
d$c <- sample(c(0, 1), n, TRUE, prob = c(0.8, 0.2))
d$lb <- -3
try1 <- function(lab, expr) {
  r <- tryCatch({
    f <- expr
    p <- tryCatch(dim(predict(f, newdata = d[1:5, ])), error = function(e)
      paste("PREDICT ERR:", conditionMessage(e)))
    ft <- tryCatch(dim(fitted(f, newdata = d[1:5, ])), error = function(e)
      paste("FITTED ERR:", conditionMessage(e)))
    sm <- tryCatch(dim(as.data.frame(simulate(f, nsim = 2, seed = 1,
                                               newdata = d[1:5, ]))),
                   error = function(e) paste("SIM ERR:", conditionMessage(e)))
    ce <- tryCatch(class(conditional_effects(f, effects = "x"))[1],
                   error = function(e) paste("CE ERR:", conditionMessage(e)))
    paste("OK ll=", format(as.numeric(logLik(f)), digits = 10),
          "| pred", paste(p, collapse = "x"), "| fit",
          paste(ft, collapse = "x"), "| sim", paste(sm, collapse = "x"),
          "| ce", ce)
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-28s %s\n", lab, r))
}
try1("weights(wt * 2)", frm(bf(y | weights(wt * 2) ~ x), data = d, family = gaussian()))
try1("rate(time * 2)", frm(bf(yc | rate(time * 2) ~ x), data = d, family = poisson()))
try1("trials(n + 1)", frm(bf(yb | trials(n + 1) ~ x), data = d, family = binomial()))
try1("se(s / 2)", frm(bf(y | se(s / 2) ~ x), data = d, family = gaussian()))
try1("cens(c == 1)", frm(bf(y | cens(c == 1) ~ x), data = d, family = gaussian()))
try1("cens(ifelse)", frm(bf(y | cens(ifelse(c == 1, 'right', 'none')) ~ x),
                         data = d, family = gaussian()))
try1("trunc(lb = lb - 1)", frm(bf(y | trunc(lb = lb - 1) ~ x), gaussian(),
                                data = d))
try1("trunc(lb = min(y) - 1)", frm(bf(y | trunc(lb = min(y) - 1) ~ x),
                                    data = d, family = gaussian()))
try1("subset(x > 0)", frm(bf(y | subset(x > 0) ~ x), data = d, family = gaussian()))
try1("weights(wt)", frm(bf(y | weights(wt) ~ x), data = d, family = gaussian()))
