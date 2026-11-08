# What df = NA under a prior does to the readers of logLik() (fix round,
# review M2): print, AIC spellings, anova, lmtest::lrtest, frm_allfit,
# frm_multiple pooling, and frmtmb.sample's reads. Seed 201.
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
try_show <- function(label, expr) {
  r <- tryCatch(expr, error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-28s ERROR <%s> %s\n", label, class(r)[1L],
                substr(conditionMessage(r), 1, 110)))
  } else {
    cat(sprintf("%-28s %s\n", label,
                paste(utils::capture.output(print(r)), collapse = " | ")))
  }
}
set.seed(201)
d <- data.frame(x = rnorm(120), z = rnorm(120),
                g = factor(rep(1:12, each = 10)))
d$y <- 1 + 0.5 * d$x + rnorm(12)[d$g] + rnorm(120)
pr <- set_prior("normal(0, 0.1)", class = "b")
fl <- lm(y ~ x, data = d)
f0 <- frm(y ~ x + (1 | g), data = d)
fp <- frm(y ~ x + (1 | g), data = d, prior = pr)
fp2 <- frm(y ~ x + z + (1 | g), data = d, prior = pr)
try_show("print(logLik(fp))", logLik(fp))
try_show("attr df_why", attr(logLik(fp), "df_why"))
try_show("AIC(fp)", AIC(fp))
try_show("AIC(lm, fp)", AIC(fl, fp))
try_show("BIC(lm, fp)", BIC(fl, fp))
try_show("AIC(logLik(fp))", AIC(logLik(fp)))
try_show("BIC(logLik(fp))", BIC(logLik(fp)))
try_show("anova(fp, fp2)", anova(fp, fp2))
try_show("anova(f0, fp)", anova(f0, fp))
if (requireNamespace("lmtest", quietly = TRUE)) {
  try_show("lmtest::lrtest(fp, fp2)", lmtest::lrtest(fp, fp2))
} else cat("lmtest not installed\n")
af <- frm_allfit(fp, optimizers = list(nlminb = "nlminb", optim = "optim"))
try_show("print(frm_allfit(fp))", af)
imps <- lapply(1:3, function(i) {
  dd <- d
  dd$x <- dd$x + rnorm(nrow(dd), 0, 0.05)
  dd
})
m1 <- frm_multiple(bf(y ~ x + z + (1 | g)), data = imps, prior = pr)
m0 <- frm_multiple(bf(y ~ x + (1 | g)), data = imps, prior = pr)
try_show("summary(m1) coef rows", nrow(summary(m1)$coefficients %||%
                                         summary(m1)$fixed))
try_show("anova(m1, m0) D1", anova(m1, m0, method = "D1"))
try_show("anova(m1, m0) D3", anova(m1, m0))
try_show("logLik_penalized(m1)", logLik_penalized(m1))
try_show("df.residual(fp)", df.residual(fp))
try_show("summary(fp) Method line",
         grep("Method", utils::capture.output(print(summary(fp))),
              value = TRUE))
if (requireNamespace("frmtmb.sample", quietly = TRUE)) {
  ns <- asNamespace("frmtmb.sample")
  hits <- Filter(function(nm) {
    f <- get(nm, ns)
    is.function(f) && any(grepl("logLik(", deparse(f), fixed = TRUE))
  }, ls(ns, all.names = TRUE))
  cat("frmtmb.sample functions whose body calls logLik(:",
      if (length(hits)) paste(hits, collapse = ", ") else "none", "\n")
}
cat("done\n")
