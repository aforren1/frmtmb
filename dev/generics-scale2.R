# Item 2.5d, round 2: the rows the contract page was missing, and the
# two conditions the lognormal identity was stated without.
#
# Review found three holes and each is measured here rather than
# reasoned about:
#   * conditional_effects() is absent from the table and is on the
#     RESPONSE scale, the opposite of predict()'s default.
#   * the identity fitted() == exp(mu + sigma^2/2) is conditional on a
#     CONSTANT sigma and on NO truncation.
#   * posterior_summary(), pp_check(), bayes_R2() and hypothesis() are
#     stated nowhere.
#
#   Rscript dev/generics-scale2.R [LIB]
av <- commandArgs(trailingOnly = TRUE)
LIB <- if (length(av)) av[1] else "C:/Users/adf44/source/r/generics-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))

set.seed(2026)
n <- 400
dd <- data.frame(x = rnorm(n), g = factor(rep(1:20, 20)))
eta <- 8 + 0.4 * dd$x + rnorm(20, 0, 0.3)[dd$g]
dd$y <- exp(rnorm(n, eta, 0.4))
fit <- frm(bf(y ~ x + (1 | g)) + lognormal(), data = dd)

cat("```\n")
cat("== the rows the contract was missing, dev/generics-scale2.R ==\n")
cat("same 400-row lognormal fit, seed 2026\n\n")

cat("-- conditional_effects(), the row that matters most --\n")
ce <- conditional_effects(fit)
e1 <- ce[[1]]
cat(sprintf("effect plotted                  %s\n", names(ce)[1]))
cat(sprintf("estimate__ range                %.2f to %.2f\n",
            min(e1$estimate__), max(e1$estimate__)))
cat(sprintf("predict() on the same fit       %.4f\n", predict(fit)[1]))
cat(sprintf("fitted() range                  %.2f to %.2f\n",
            min(fitted(fit)), max(fitted(fit))))
cat("so conditional_effects() is RESPONSE, the opposite of predict()\n")

cat("\n-- the identity, and the two conditions it is stated without --\n")
mu <- predict(fit, type = "link")
sg <- sigma(fit)
cat(sprintf("constant sigma: identical(fitted, exp(mu+s^2/2))  %s\n",
            identical(fitted(fit), exp(mu + sg^2 / 2))))

# condition 1: a distributional sigma
f2 <- frm(bf(y ~ x + (1 | g), sigma ~ x) + lognormal(), data = dd)
s2 <- tryCatch(suppressWarnings(sigma(f2)), error = function(e) NA)
cat(sprintf("distributional sigma: sigma(fit) is              %s\n",
            paste(s2, collapse = ",")))
mu2 <- predict(f2, type = "link")
sv <- predict(f2, dpar = "sigma", type = "response")
cat(sprintf("  identity with sigma():   %s\n",
            if (all(is.na(exp(mu2 + s2^2 / 2)))) "NA, the formula fails"
            else "ok"))
cat(sprintf("  identity with predict(dpar = 'sigma'): identical %s\n",
            identical(fitted(f2), exp(mu2 + sv^2 / 2))))

# condition 2: truncation
dd3 <- dd[dd$y > 2000, ]
f3 <- frm(bf(y | trunc(lb = 2000) ~ x + (1 | g)) + lognormal(),
          data = dd3)
mu3 <- predict(f3, type = "link")
s3 <- sigma(f3)
naive <- exp(mu3 + s3^2 / 2)
rel <- max(abs(fitted(f3) - naive)) / max(abs(fitted(f3)))
cat(sprintf("truncated y | trunc(lb = 2000): max relative error of\n"))
cat(sprintf("  exp(mu + sigma^2/2) against fitted()            %.4f\n",
            rel))
cat(sprintf("  fitted()[1] %.2f against the naive formula %.2f\n",
            fitted(f3)[1], naive[1]))

cat("\n-- the four methods stated nowhere --\n")
ps <- tryCatch(posterior_summary(fit), error = function(e) conditionMessage(e))
cat(sprintf("posterior_summary(fit)          %s\n",
            if (is.character(ps)) paste("ERR:", substr(ps, 1, 46))
            else "matrix"))
cat(sprintf("posterior_summary(matrix)       %s\n",
            paste(colnames(posterior_summary(cbind(a = rnorm(50)))),
                  collapse = ",")))
for (g in c("pp_check", "bayes_R2")) {
  m <- tryCatch({ do.call(g, list(fit)); "returned" },
                error = function(e) paste("refuses:",
                                          substr(conditionMessage(e), 1, 40)))
  cat(sprintf("%-31s %s\n", paste0(g, "(fit)"), m))
}
h <- hypothesis(fit, "x = 0")
hd <- if (is.data.frame(h)) h else h$hypothesis
cat(sprintf("hypothesis() returns %s with columns %s
",
            paste(class(h), collapse = "/"),
            paste(names(hd), collapse = ",")))
est <- hd[[grep("stimate|^Est", names(hd))[1]]]
cat(sprintf("hypothesis(fit, 'x = 0') Estimate %.6f  vs fixef %.6f\n",
            est[1], fixef(fit)$mu[["x"]]))
cat("```\n")
