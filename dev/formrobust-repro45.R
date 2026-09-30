# Items 4 and 5: update() of a nonlinear fit with a new body, and
# bernoulli on a two-valued response that is not 0/1. Seed 41.
# FORMROBUST_LIB="" for the before arm.
LIB <- Sys.getenv("FORMROBUST_LIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
tr <- function(lab, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-34s %s\n", lab, paste(format(r, digits = 10),
                                       collapse = " ")))
}
set.seed(41)
n <- 80
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$y <- 2 + 1.5 * d$x + rnorm(n, 0, 0.5)
f0 <- frm(bf(y ~ a + b * x, a ~ 1 + z, b ~ 1, nl = TRUE), data = d)
up <- tryCatch(update(f0, formula. = bf(y ~ a + b * x + 0, nl = TRUE)),
               error = function(e) conditionMessage(e))
tr("update nl, new body keeps a,b", if (is.character(up)) stop(up) else
  c(identical(logLik(up)[1], logLik(f0)[1]), formula(up)$pforms$a))
up2 <- tryCatch(update(f0, formula. = bf(y ~ a * exp(b * x), nl = TRUE)),
                error = function(e) conditionMessage(e))
tr("update nl, exp body", if (is.character(up2)) stop(up2) else
  names(fixef(up2)[, 1]))
tr("bf(y ~ a + b, nl = TRUE) alone", class(bf(y ~ a + b, nl = TRUE))[1])
tr("frm(bf(y ~ a + b, nl = TRUE))",
   frm(bf(y ~ a + b * x, nl = TRUE), data = d))
flin <- frm(bf(y ~ x, sigma ~ z), data = d)
tr("update linear keeps sigma ~ z",
   rownames(fixef(update(flin, y ~ x + z))))

# item 5
set.seed(5)
n <- 120
d5 <- data.frame(x = rnorm(n))
p <- plogis(-0.3 + 1.2 * d5$x)
yy <- rbinom(n, 1, p)
d5$y01 <- yy
d5$ym <- ifelse(yy == 1, -1, -2)            # -2 codes 0, -1 codes 1
d5$yf <- factor(ifelse(yy == 1, "yes", "no"))
d5$yl <- yy == 1
d5$y12 <- yy + 1
d5$yc <- ifelse(yy == 1, "yes", "no")
ref <- frm(bf(y01 ~ x), data = d5, family = bernoulli())
for (v in c("ym", "yf", "yl", "y12", "yc")) {
  f <- tryCatch(frm(bf(as.formula(paste(v, "~ x"))), data = d5,
                    family = bernoulli()),
                error = function(e) conditionMessage(e))
  tr(paste("bernoulli", v, "logLik identical"),
     if (is.character(f)) stop(f) else
       identical(logLik(f)[1], logLik(ref)[1]))
}
fm <- frm(bf(ym ~ x), data = d5, family = bernoulli())
nd <- d5[d5$ym == -1, ][1:5, ]
tr("residuals(newdata) on -1 rows",
   isTRUE(all.equal(residuals(fm, newdata = nd)[, "Estimate"],
                    residuals(ref, newdata = nd)[, "Estimate"],
                    tolerance = 0)))
tr("influence() on -1/-2", {
  i1 <- influence(fm); i0 <- influence(ref)
  identical(i1$fixef, i0$fixef) || isTRUE(all.equal(i1, i0, tolerance = 0))
})
tr("simulate codes", sort(unique(as.vector(unlist(simulate(fm, nsim = 2,
                                                             seed = 1))))))
tr("three values refused",
   frm(bf(y ~ x), data = data.frame(x = 1:6, y = c(0, 1, 2, 0, 1, 2)),
       family = bernoulli()))
