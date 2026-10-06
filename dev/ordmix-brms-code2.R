# brms 2.23.0 on more ordinal mixture combinations: order = "mu" (fixed
# thresholds) with each structure and family pair, and the refusals.
# Read-only survey. Output: dev/ordmix-log-brms-code2.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
set.seed(1)
n <- 60
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
d$y <- sample(1:4, n, TRUE)
d$yh <- sample(0:4, n, TRUE)
blocks <- function(code) {
  lines <- strsplit(as.character(code), "\n")[[1]]
  st <- grep("^(functions|data|transformed data|parameters|transformed parameters|model|generated quantities) \\{", lines)
  out <- character()
  for (i in seq_along(st)) {
    nm <- sub(" \\{.*", "", lines[st[i]])
    en <- if (i < length(st)) st[i + 1] - 1 else length(lines)
    if (nm %in% c("parameters", "transformed parameters", "model",
                  "generated quantities")) out <- c(out, lines[st[i]:en])
  }
  out <- out[!grepl("^ *//", out)]
  out
}
show <- function(label, f, fam, data = d) {
  cat("\n################ ", label, " ################\n", sep = "")
  code <- tryCatch(stancode(f, data = data, family = fam),
                   error = function(e) e)
  if (inherits(code, "error")) {
    cat("ERROR:", conditionMessage(code), "\n")
    return(invisible())
  }
  cat(blocks(code), sep = "\n")
  dp <- tryCatch(default_prior(f, data = data, family = fam),
                 error = function(e) e)
  cat("-- default_prior --\n")
  if (inherits(dp, "error")) cat("ERROR:", conditionMessage(dp), "\n") else
    print(as.data.frame(dp)[, c("prior", "class", "coef", "group", "resp",
                                "dpar", "lb", "ub", "source")])
}
show("order mu, cumulative + sratio", bf(y ~ x),
     mixture(cumulative, sratio, order = "mu"))
show("order mu, sratio + acat", bf(y ~ x),
     mixture(sratio, acat, order = "mu"))
show("order TRUE, cumulative x2", bf(y ~ x),
     mixture(cumulative, cumulative, order = TRUE))
show("order mu, sum_to_zero", bf(y ~ x),
     mixture(cumulative(threshold = "sum_to_zero"),
             cumulative(threshold = "sum_to_zero"), order = "mu"))
show("order mu, equidistant", bf(y ~ x),
     mixture(cumulative(threshold = "equidistant"),
             cumulative(threshold = "equidistant"), order = "mu"))
show("order mu, thres(gr)", bf(y | thres(gr = g) ~ x),
     mixture(cumulative, cumulative, order = "mu"))
show("order mu, disc1 ~ z, mu2 ~ z", bf(y ~ x, disc1 ~ 0 + z, mu2 ~ z),
     mixture(cumulative, cumulative, order = "mu"))
show("order mu, cs", bf(y ~ cs(x)), mixture(sratio, sratio, order = "mu"))
show("order mu, y ~ 1", bf(y ~ 1), mixture(cumulative, cumulative,
                                          order = "mu"))
show("order none, y ~ 1", bf(y ~ 1), mixture(cumulative, cumulative))
show("cumulative + poisson", bf(y ~ x), mixture(cumulative, poisson))
show("cumulative + categorical", bf(y ~ x), mixture(cumulative,
                                                   categorical))
show("order mu, mixed flexible + sum_to_zero", bf(y ~ x),
     mixture(cumulative, cumulative(threshold = "sum_to_zero"),
             order = "mu"))
show("order none, sum_to_zero + flexible", bf(y ~ x),
     mixture(cumulative(threshold = "sum_to_zero"), cumulative))
show("thres(gr) + equidistant mixture", bf(y | thres(gr = g) ~ x),
     mixture(cumulative(threshold = "equidistant"), cumulative))
show("theta1 ~ z order mu", bf(y ~ x, theta1 ~ z),
     mixture(cumulative, cumulative, order = "mu"))
show("disc1 ~ x, probit", bf(y ~ x, disc1 ~ 0 + x),
     mixture(cumulative("probit"), cumulative("probit")))
show("nl mu1", bf(y ~ 1, nlf(mu1 ~ a * x), a ~ 1, mu2 ~ x),
     mixture(cumulative, cumulative))
show("hurdle x2 order mu", bf(yh ~ x),
     mixture(hurdle_cumulative, hurdle_cumulative, order = "mu"),
     data = d)
show("hurdle x2, hu1 ~ z", bf(yh ~ x, hu1 ~ z),
     mixture(hurdle_cumulative, hurdle_cumulative), data = d)
show("mixture cum cum mvbf", mvbf(bf(y ~ x), bf(yh ~ x)) ,
     mixture(cumulative, cumulative))
