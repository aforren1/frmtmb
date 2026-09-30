# brms 2.23.0's Stan code, data and default priors for the ordinal
# families with disc and with each threshold option. Read-only survey;
# nothing compiles. Output: dev/ordinal-log-brms-code.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
set.seed(1)
n <- 60
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)),
                h = factor(sample(c("p", "q", "r"), n, TRUE)))
d$y <- sample(1:4, n, TRUE)
d$y2 <- sample(1:3, n, TRUE)
d$yh <- sample(0:4, n, TRUE)
blocks <- function(code, keep = c("parameters", "transformed parameters",
                                  "model", "generated quantities")) {
  code <- as.character(code)
  lines <- strsplit(code, "\n")[[1]]
  st <- grep("^(data|transformed data|parameters|transformed parameters|model|generated quantities) \\{", lines)
  out <- character()
  for (i in seq_along(st)) {
    nm <- sub(" \\{.*", "", lines[st[i]])
    en <- if (i < length(st)) st[i + 1] - 1 else length(lines)
    if (nm %in% keep) out <- c(out, lines[st[i]:en])
  }
  out
}
show <- function(label, f, fam, data = d, prior = NULL, full = FALSE) {
  cat("\n################ ", label, " ################\n", sep = "")
  code <- tryCatch(stancode(f, data = data, family = fam, prior = prior),
                   error = function(e) e)
  if (inherits(code, "error")) {
    cat("ERROR:", conditionMessage(code), "\n")
    return(invisible())
  }
  cat(if (full) as.character(code) else blocks(code), sep = "\n")
  dp <- default_prior(f, data = data, family = fam)
  cat("-- default_prior --\n")
  print(as.data.frame(dp)[, c("prior", "class", "coef", "group", "resp",
                              "dpar", "lb", "ub", "source")])
  sd <- standata(f, data = data, family = fam, prior = prior)
  cat("-- standata names --\n")
  print(names(sd))
  for (k in grep("nthres|disc|Jthres|ngrthres|nmthres", names(sd),
                 value = TRUE)) {
    cat(k, ":", head(sd[[k]], 6), "\n")
  }
}
show("cumulative flexible", bf(y ~ x), cumulative())
show("cumulative disc ~ z", bf(y ~ x, disc ~ z), cumulative())
show("cumulative disc ~ 0 + z", bf(y ~ x, disc ~ 0 + z), cumulative())
show("cumulative disc ~ 1 (intercept only)", bf(y ~ x, disc ~ 1),
     cumulative())
show("cumulative fixed disc = 2", bf(y ~ x, disc = 2), cumulative())
show("cumulative disc ~ z, probit", bf(y ~ x, disc ~ z),
     cumulative("probit"))
show("cumulative disc ~ z, link_disc identity", bf(y ~ x, disc ~ z),
     cumulative(link_disc = "identity"))
show("cumulative equidistant", bf(y ~ x),
     cumulative(threshold = "equidistant"), full = TRUE)
show("cumulative sum_to_zero", bf(y ~ x),
     cumulative(threshold = "sum_to_zero"), full = TRUE)
show("sratio equidistant cse(z)", bf(y ~ x + cse(z)),
     sratio(threshold = "equidistant"))
show("sratio sum_to_zero cs(z)", bf(y ~ x + cs(z)),
     sratio(threshold = "sum_to_zero"))
show("cratio equidistant disc ~ z", bf(y ~ x, disc ~ z),
     cratio(threshold = "equidistant"))
show("acat sum_to_zero disc ~ 0 + z", bf(y ~ x, disc ~ 0 + z),
     acat(threshold = "sum_to_zero"))
show("acat equidistant", bf(y ~ x), acat(threshold = "equidistant"))
show("cumulative thres(gr) flexible", bf(y | thres(gr = g) ~ x),
     cumulative())
show("cumulative thres(gr) equidistant", bf(y | thres(gr = g) ~ x),
     cumulative(threshold = "equidistant"), full = TRUE)
show("cumulative thres(gr) sum_to_zero", bf(y | thres(gr = g) ~ x),
     cumulative(threshold = "sum_to_zero"))
show("sratio thres(gr) equidistant disc ~ z",
     bf(y | thres(gr = g) ~ x, disc ~ z), sratio(threshold = "equidistant"))
show("cumulative thres(gr) + cs", bf(y | thres(gr = g) ~ x + cs(z)),
     sratio())
show("hurdle_cumulative equidistant", bf(yh ~ x),
     hurdle_cumulative(threshold = "equidistant"))
show("hurdle_cumulative sum_to_zero", bf(yh ~ x),
     hurdle_cumulative(threshold = "sum_to_zero"))
show("hurdle_cumulative equidistant disc ~ 0 + z", bf(yh ~ x, disc ~ 0 + z),
     hurdle_cumulative(threshold = "equidistant"))
show("cumulative thres(x = 5) equidistant", bf(y | thres(5) ~ x),
     cumulative(threshold = "equidistant"))
show("cumulative 2 categories equidistant", bf(y2b ~ x),
     cumulative(threshold = "equidistant"),
     data = transform(d, y2b = ifelse(y > 2, 2L, 1L)))
show("cumulative 2 categories sum_to_zero", bf(y2b ~ x),
     cumulative(threshold = "sum_to_zero"),
     data = transform(d, y2b = ifelse(y > 2, 2L, 1L)))
show("mv two ordinal equidistant",
     bf(y ~ x, family = cumulative(threshold = "equidistant")) +
       bf(y2 ~ x, disc ~ z, family = sratio(threshold = "equidistant")) +
       set_rescor(FALSE), NULL)
show("mixture cumulative cumulative", bf(y ~ x),
     mixture(cumulative, cumulative))
show("mixture cumulative equidistant", bf(y ~ x),
     mixture(cumulative(threshold = "equidistant"),
             cumulative(threshold = "equidistant")))
show("mixture cumulative disc1 ~ z", bf(y ~ x, disc1 ~ z),
     mixture(cumulative, cumulative))
show("cumulative disc ~ (1|g) center", bf(y ~ x, disc ~ 0 + z + (1 | g)),
     cumulative())
show("cumulative center = FALSE sum_to_zero", bf(y ~ x, center = FALSE),
     cumulative(threshold = "sum_to_zero"))
show("cumulative mo(), sum_to_zero",
     bf(y ~ x + mo(hm)), cumulative(threshold = "sum_to_zero"),
     data = transform(d, hm = as.integer(h)))
cat("\n== validate_prior: delta and Intercept under equidistant ==\n")
print(validate_prior(c(prior(normal(0, 1), class = delta),
                       prior(normal(0, 3), class = Intercept)),
                     bf(y ~ x), data = d,
                     family = cumulative(threshold = "equidistant")))
cat("\n== coef-level Intercept prior under equidistant ==\n")
print(tryCatch(validate_prior(prior(normal(0, 3), class = Intercept,
                                    coef = 1),
                              bf(y ~ x), data = d,
                              family = cumulative(threshold = "equidistant")),
               error = function(e) conditionMessage(e)))
cat("\n== delta prior on flexible ==\n")
print(tryCatch(validate_prior(prior(normal(0, 1), class = delta),
                              bf(y ~ x), data = d, family = cumulative()),
               error = function(e) conditionMessage(e)))
cat("\n== delta prior grouped ==\n")
print(tryCatch(validate_prior(prior(normal(0, 1), class = delta, group = a),
                              bf(y | thres(gr = g) ~ x), data = d,
                              family = cumulative(threshold = "equidistant")),
               error = function(e) conditionMessage(e)))
cat("\n== threshold argument refusals ==\n")
print(tryCatch(cumulative(threshold = "foo"), error = function(e)
  conditionMessage(e)))
print(tryCatch(cumulative(threshold = "equi"), error = function(e)
  conditionMessage(e)))
print(tryCatch(cumulative(link_disc = "logit"), error = function(e)
  conditionMessage(e)))
print(tryCatch(stancode(bf(y ~ x, disc = -1), data = d,
                        family = cumulative()),
               error = function(e) conditionMessage(e)))
cat("\n== brmsfamily fields ==\n")
str(unclass(cumulative(threshold = "equidistant")))
str(unclass(sratio(link_disc = "identity")))
