## y ~ mo(k) + mok, where the data have a column named mok: brms's
## spelling of the mo() coefficient and of the plain column coincide.
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-brmsnames-lib", base)
} else base)
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(8)
n <- 200
d <- data.frame(k = sample(0:3, n, TRUE), mok = rnorm(n))
d$y <- 1 + 0.5 * d$k + 0.3 * d$mok + rnorm(n)
bp <- tryCatch(as.data.frame(q(brms::default_prior(y ~ mo(k) + mok,
                                                   data = d))),
               error = function(e) conditionMessage(e))
cat("brms default_prior:\n"); print(if (is.data.frame(bp))
  bp[, c("prior", "class", "coef")] else bp)
q(library(frmtmb))
cat("arm", arm, "\n")
fp <- tryCatch(as.data.frame(q(default_prior(y ~ mo(k) + mok, data = d))),
               error = function(e) conditionMessage(e))
print(if (is.data.frame(fp)) fp[, c("prior", "class", "coef")] else fp)
fit <- tryCatch(q(frm(y ~ mo(k) + mok, data = d,
                      prior = set_prior("normal(3, 0.001)", coef = "mok"))),
                error = function(e) conditionMessage(e))
if (is.character(fit)) cat("fit:", fit, "\n") else print(round(fixef(fit)[, 1:2], 3))
