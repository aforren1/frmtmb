# Readers of a fit whose missing rows left the objective. Seed 11.
.libPaths(c("C:/Users/adf44/source/r/wt-refuse-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(11)
n <- 120
d <- data.frame(x = rnorm(n))
d$y <- rgamma(n, 3, 3 / exp(0.5 + 0.4 * d$x))
d$y[1:10] <- NA
f <- frm(bf(y | mi() ~ x), family = stats::Gamma(link = "log"), data = d)
g <- frm(bf(y ~ x), family = stats::Gamma(link = "log"),
         data = d[!is.na(d$y), ])
show <- function(lab, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 120)))
  cat(sprintf("%-20s %s\n", lab, paste(format(r, digits = 5), collapse = " ")))
}
show("logLik gap", as.numeric(logLik(f)) - as.numeric(logLik(g)))
show("nobs", c(nobs(f), nobs(g)))
show("attr nobs", attr(logLik(f), "nobs"))
show("cooks at missing", head(cooks.distance(f), 3))
show("cooks observed gap", max(abs(cooks.distance(f)[-(1:10)] - cooks.distance(g))))
show("resid at missing", head(residuals(f)[, 1], 3))
show("fitted at missing", head(fitted(f)[, 1], 2))
show("predict at missing", head(predict(f)[, 1], 2))
show("sim len", length(simulate(f, nsim = 1)[[1]]))
show("confint", confint(f)[1, ])
show("confint cc", confint(g)[1, ])
b <- tryCatch(frm_bootstrap(f, nsim = 4, seed = 1), error = function(e) e)
show("bootstrap", if (inherits(b, "error")) conditionMessage(b) else class(b)[1])
show("update", as.numeric(logLik(update(f))) - as.numeric(logLik(f)))
q <- tryCatch(frm(bf(y | mi() ~ x + (1 | g)), family = stats::Gamma(link = "log"),
                  data = transform(d, g = factor(rep(1:10, 12))), quadrature = TRUE),
              error = function(e) e)
show("quadrature", if (inherits(q, "error")) conditionMessage(q) else as.numeric(logLik(q)))
show("par_template", names(par_template(f)))
pp <- frm_multiple
show("summary", class(summary(f))[1])
show("dry_run frame", length(frm(bf(y | mi() ~ x), family = stats::Gamma(link = "log"),
                                 data = d, dry_run = "frame")$par_template$miss))
