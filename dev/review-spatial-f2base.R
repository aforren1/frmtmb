# Reviewer: the same formula fits on a given library (arg 1), for the logLik and names.
.libPaths(c(strsplit(commandArgs(TRUE)[1], ";")[[1]], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", as.character(packageVersion("frmtmb")), dirname(find.package("frmtmb")), "\n")
r <- function(lab, e) cat(sprintf("%-30s -> %s\n", lab, tryCatch({
  f <- suppressWarnings(suppressMessages(e))
  paste("logLik", format(as.numeric(logLik(f)), digits = 10), "b:", paste(names(f$estimates$beta), collapse = ","),
        "| fixef:", paste(rownames(fixef(f)), collapse = ","), "| variables:", paste(head(variables(f), 2), collapse = ","))
}, error = function(e) paste("ERR", substr(conditionMessage(e), 1, 160)))))
set.seed(3)
d <- data.frame(g = gl(15, 6), t = factor(rep(1:3, 30)), x = runif(90), z = runif(90, 1, 2))
d$y <- d$x + rnorm(15)[d$g] + rnorm(90)
r("y ~ x + exp(0 + t | g)", frm(y ~ x + exp(0 + t | g), data = d))
r("y ~ exp(x) + exp(0 + t | g)", frm(y ~ exp(x) + exp(0 + t | g), data = d))
r("y ~ exp(x) + (1 | g)", frm(y ~ exp(x) + (1 | g), data = d))
f <- frm(y ~ exp(x) + (1 | g), data = d)
p <- predict(f, newdata = data.frame(x = c(0, 1), g = "1"))
cat("predict exp(x) newdata:", format(p[, 1], digits = 8), "\n")
