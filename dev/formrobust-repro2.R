# Item 2 probe: offset() in conditional_effects() and emmeans(), in
# several positions. Seed 21. FORMROBUST_LIB="" for the before arm.
LIB <- Sys.getenv("FORMROBUST_LIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(emmeans))
cat("frmtmb from", find.package("frmtmb"), "\n")
set.seed(21)
n <- 80
d <- data.frame(x = rnorm(n), time = runif(n, 1, 3), f = gl(2, 40),
                z = rnorm(n))
d$yc <- rpois(n, exp(0.2 + 0.3 * d$x) * d$time)
d$y <- rnorm(n, 0.5 * d$x + d$z, exp(0.1 * d$x) * d$time)
run <- function(lab, expr) {
  f <- tryCatch(expr, error = function(e) e)
  if (inherits(f, "error")) {
    cat(sprintf("%-26s FIT ERROR: %s\n", lab, conditionMessage(f)))
    return(invisible())
  }
  ce <- tryCatch({
    o <- conditional_effects(f, effects = "x")
    paste("ok", format(signif(o[[1]]$estimate__[1], 8)))
  }, error = function(e) paste("ERR", conditionMessage(e)))
  em <- tryCatch({
    o <- summary(emmeans(f, ~ f))
    paste("ok", format(signif(o$emmean[1], 8)))
  }, error = function(e) paste("ERR", conditionMessage(e)))
  cat(sprintf("%-26s CE: %s\n%-26s EMM: %s\n", lab, ce, "", em))
}
run("offset(log(time))",
    frm(bf(yc ~ x + f + offset(log(time))), data = d, family = poisson()))
run("offset(time)",
    frm(bf(yc ~ x + f + offset(time)), data = d, family = poisson()))
run("offset(log(time) + 0.1)",
    frm(bf(yc ~ x + f + offset(log(time) + 0.1)), data = d,
        family = poisson()))
run("sigma ~ offset(log(time))",
    frm(bf(y ~ x + f, sigma ~ x + offset(log(time))), data = d))
run("nl a ~ offset(z)",
    frm(bf(y ~ a + b * x, a ~ f + offset(z), b ~ 1, nl = TRUE), data = d))
run("no offset",
    frm(bf(yc ~ x + f), data = d, family = poisson()))
