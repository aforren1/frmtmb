# Constructions for the candidates the spot check raised. Each one is
# run to the point where the divergence is OBSERVABLE, because "no test
# reaches it" is not evidence and neither is a mismatched regexp.
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "\n\n")

show <- function(label, expr) {
  cat("---", label, "---\n")
  out <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  print(out)
  cat("\n")
}

cat("### D1: family$link is a list, not a string\n")
show("str(student()$link)", utils::str(student()$link))
show("class(gaussian()$link)", class(frmtmb::gaussian()$link))
show("names(student())", names(student()))
show("student()$link_shape", student()$link_shape)

cat("### D2: an unsupported link is accepted by the constructor\n")
show("bernoulli('sqrt')$link", bernoulli("sqrt")$link)
show("exponential('cloglog')$link", exponential("cloglog")$link)
show("Beta('1/mu^2')$link", Beta("1/mu^2")$link)
# and does the fit use it?
set.seed(1)
d <- data.frame(y = rbinom(60, 1, 0.4), x = rnorm(60))
show("frm(y ~ x, bernoulli('sqrt')) coefficients",
     {
       f <- frm(y ~ x, d, family = bernoulli("sqrt"))
       list(coef = fixef(f), link = f$family$link)
     })
show("frm(y ~ x, bernoulli()) coefficients",
     {
       f <- frm(y ~ x, d, family = bernoulli())
       fixef(f)
     })

cat("### D3: bf() is not idempotent\n")
show("bf(bf(y ~ x, sigma ~ 1))", {
  form <- bf(y ~ x, sigma ~ 1)
  identical(form, bf(form))
})

cat("### D4: an UNORDERED factor response is accepted by an ordinal family\n")
show("cratio on factor(c('low','mid','high'))", {
  dd <- data.frame(y = factor(rep(c("low", "mid", "high"), 20),
                              levels = c("low", "mid", "high")),
                   x = rnorm(60))
  dd$y_unordered <- factor(as.character(dd$y))  # levels alphabetical
  fr <- frm(y_unordered ~ x, dd, family = cratio(), dry_run = "frame")
  list(class = class(dd$y_unordered),
       levels = levels(dd$y_unordered),
       response_codes = sort(unique(as.vector(fr$y))))
})

cat("### D5: duplicated group-level terms on one factor\n")
show("(1|g) + (x|g)", {
  dd <- data.frame(y = rnorm(50), x = rnorm(50), g = rep(1:5, 10))
  fr <- frm(y ~ x + (1 | g) + (x | g), dd, dry_run = "frame")
  "accepted"
})

cat("### D6: family given as c(family, link)\n")
show("c('weibull','sqrt')", frm_family(c("weibull", "sqrt")))

cat("### D7: mixture components with incompatible support\n")
show("mixture(poisson(), categorical())",
     class(mixture(poisson(), categorical())))
show("mixture(lognormal, exgaussian, poisson())",
     class(mixture(lognormal, exgaussian, poisson())))
