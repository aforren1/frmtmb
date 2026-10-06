# brms 2.23.0's Stan code, data and default priors for ordinal mixtures
# and for hurdle_cumulative() with thres(gr = ) and cs(). Read-only
# survey; nothing compiles. Output: dev/ordmix-log-brms-code.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(1)
n <- 60
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
d$y <- sample(1:4, n, TRUE)
d$yh <- sample(0:4, n, TRUE)
show <- function(label, f, fam, data = d, prior = NULL, full = TRUE) {
  cat("\n################ ", label, " ################\n", sep = "")
  code <- tryCatch(stancode(f, data = data, family = fam, prior = prior),
                   error = function(e) e)
  if (inherits(code, "error")) {
    cat("ERROR:", conditionMessage(code), "\n")
    return(invisible())
  }
  cat(as.character(code), sep = "\n")
  dp <- tryCatch(default_prior(f, data = data, family = fam),
                 error = function(e) e)
  cat("-- default_prior --\n")
  if (inherits(dp, "error")) cat("ERROR:", conditionMessage(dp), "\n") else
    print(as.data.frame(dp)[, c("prior", "class", "coef", "group", "resp",
                                "dpar", "lb", "ub", "source")])
  sd <- standata(f, data = data, family = fam, prior = prior)
  cat("-- standata names --\n")
  print(names(sd))
  for (k in grep("nthres|disc|Jthres|ngrthres|con_theta", names(sd),
                 value = TRUE)) {
    cat(k, ":", head(sd[[k]], 6), "\n")
  }
}
mc <- mixture(cumulative, cumulative)
show("mix cum cum", bf(y ~ x), mc)
show("mix cum cum order none", bf(y ~ x), mixture(cumulative, cumulative,
                                                  order = "none"))
show("mix cum cum order mu", bf(y ~ x), mixture(cumulative, cumulative,
                                                order = "mu"))
show("mix probit sratio", bf(y ~ x), mixture(cumulative("probit"), sratio))
show("mix cum cum disc1 ~ x", bf(y ~ x, disc1 ~ x), mc)
show("mix cum cum disc1 ~ 0 + x", bf(y ~ x, disc1 ~ 0 + x), mc)
show("mix cum cum theta1 ~ z", bf(y ~ x, theta1 ~ z), mc)
show("mix cum cum mu2 ~ z", bf(y ~ x, mu2 ~ z), mc)
show("mix cum cum sum_to_zero", bf(y ~ x),
     mixture(cumulative(threshold = "sum_to_zero"),
             cumulative(threshold = "sum_to_zero")))
show("mix cum cum equidistant", bf(y ~ x),
     mixture(cumulative(threshold = "equidistant"),
             cumulative(threshold = "equidistant")))
show("mix cum flexible + equidistant", bf(y ~ x),
     mixture(cumulative(), cumulative(threshold = "equidistant")))
show("mix cum logit + cum probit", bf(y ~ x),
     mixture(cumulative(), cumulative("probit")))
show("mix cum cum thres(gr)", bf(y | thres(gr = g) ~ x), mc)
show("mix cum cum cs", bf(y ~ cs(x)), mc)
show("mix sratio acat cs", bf(y ~ cs(x)), mixture(sratio, acat))
show("mix cum cum thres(x = 4)", bf(y | thres(4) ~ x), mc)
show("mix cum cum (1|g)", bf(y ~ x + (1 | g)), mc)
show("mix 3 cum", bf(y ~ x), mixture(cumulative, cumulative, cumulative))
show("mix cum gaussian", bf(y ~ x), mixture(cumulative, gaussian))
show("mix hurdle_cumulative x2", bf(yh ~ x),
     mixture(hurdle_cumulative, hurdle_cumulative))
show("mix cum + hurdle_cum", bf(y ~ x),
     mixture(cumulative, hurdle_cumulative))

## hurdle_cumulative extensions
hc <- hurdle_cumulative()
show("hurdle thres(gr) hu ~ z", bf(yh | thres(gr = g) ~ x, hu ~ z), hc)
show("hurdle thres(gr) probit", bf(yh | thres(gr = g) ~ x),
     hurdle_cumulative("probit"))
show("hurdle cs(x)", bf(yh ~ cs(x)), hc)
show("hurdle cs(x) probit", bf(yh ~ cs(x)), hurdle_cumulative("probit"))
show("hurdle thres(gr) disc ~ z", bf(yh | thres(gr = g) ~ x, disc ~ 0 + z),
     hurdle_cumulative("probit"))
show("hurdle thres(gr) cs", bf(yh | thres(gr = g) ~ cs(x)), hc)
show("hurdle thres(gr) equidistant", bf(yh | thres(gr = g) ~ x),
     hurdle_cumulative(threshold = "equidistant"))
show("hurdle cs + disc", bf(yh ~ cs(x), disc ~ 0 + z),
     hurdle_cumulative("probit"))
