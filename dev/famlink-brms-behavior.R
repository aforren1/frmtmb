# What brms 2.23.0 itself does on the constructions this lane changes.
# Run: Rscript dev/famlink-brms-behavior.R > dev/famlink-brms-behavior-log.txt
.libPaths(c("C:/Users/adf44/source/r/pinlib",
             "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n\n")

show <- function(label, expr) {
  msgs <- character(0)
  warns <- character(0)
  res <- withCallingHandlers(
    tryCatch({ force(expr); "OK" },
             error = function(e) paste0("ERROR [", paste(class(e), collapse = "/"),
                                        "]: ", conditionMessage(e))),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  cat("==", label, "\n  result:", res, "\n")
  for (m in msgs) cat("  message:", trimws(m), "\n")
  for (w in warns) cat("  warning:", trimws(w), "\n")
}

set.seed(1)
d <- data.frame(x = rnorm(20))
d$yb <- rbinom(20, 1, 0.5)
d$y3 <- sample(1:3, 20, TRUE)
d$y2 <- sample(1:2, 20, TRUE)
d$fac <- factor(sample(c("high", "low", "mid"), 20, TRUE))
d$ord <- factor(d$fac, ordered = TRUE)
d$fac2 <- factor(sample(c("a", "b"), 20, TRUE))
d$one <- 1
d$two <- 2

# item 11: an unordered factor response to an ordinal family
for (f in c("cumulative", "sratio", "cratio", "acat")) {
  show(paste(f, "unordered factor"),
       standata(fac ~ x, d, family = f))
  show(paste(f, "ordered factor"), standata(ord ~ x, d, family = f))
}
show("categorical unordered factor", standata(fac ~ x, d, family = "categorical"))

# item 15: the bernoulli suggestion
show("binomial trials(1)", standata(yb | trials(1) ~ x, d, family = binomial()))
show("binomial trials(one) column", standata(yb | trials(one) ~ x, d, family = binomial()))
show("binomial trials(2)", standata(yb | trials(two) ~ x, d, family = binomial()))
show("binomial no trials", standata(yb ~ x, d, family = binomial()))
show("beta_binomial trials(1)", standata(yb | trials(1) ~ x, d, family = beta_binomial()))
show("zero_inflated_binomial trials(1)",
     standata(yb | trials(1) ~ x, d, family = zero_inflated_binomial()))
for (f in c("cumulative", "sratio", "cratio", "acat")) {
  show(paste(f, "2 categories integer"), standata(y2 ~ x, d, family = f))
  show(paste(f, "3 categories integer"), standata(y3 ~ x, d, family = f))
}
show("cumulative 2-level ordered factor",
     standata(o2 ~ x, transform(d, o2 = factor(fac2, ordered = TRUE)),
              family = "cumulative"))
show("categorical 2 levels", standata(fac2 ~ x, d, family = "categorical"))
show("categorical 3 levels", standata(fac ~ x, d, family = "categorical"))

# item 16 and section 7 contract 10
show("family = c('weibull','log')", standata(exp(x) ~ 1, d, family = c("weibull", "log")))
show("validate_family c('weibull','sqrt')", brms:::validate_family(c("weibull", "sqrt")))
show("validate_family c('categorical','probit')", brms:::validate_family(c("categorical", "probit")))
show("validate_family c('poisson')", print(brms:::validate_family("poisson")$link))
show("validate_family 3 long", brms:::validate_family(c("poisson", "log", "x")))
show("validate_family c('zi_poisson', 'identity')",
     print(brms:::validate_family(c("zi_poisson", "identity"))$link))
show("brmsfamily('Gamma')", print(brmsfamily("Gamma")$family))
show("brmsfamily('normal')", print(brmsfamily("normal")$family))
show("brmsfamily('hu_poisson')", print(brmsfamily("hu_poisson")$family))
show("brmsfamily('poisson', link_sigma = 'identity')",
     print(names(brmsfamily("poisson", link_sigma = "identity"))))
show("brmsfamily('poisson', 'sqrt')$link", print(brmsfamily("poisson", "sqrt")$link))
show("brmsfamily(c('poisson','sqrt'))", print(brmsfamily(c("poisson", "sqrt"))$link))
lk <- "probit"
show("bernoulli(lk) with lk <- 'probit'", print(bernoulli(lk)$link))
show("bernoulli(link = NULL)", print(bernoulli(link = NULL)$link))
show("bernoulli(link = NA)", print(bernoulli(link = NA)$link))
show("student(link_sigma = identity) unquoted secondary", student(link_sigma = identity))
show("weibull(link_shape = 'logit')", weibull(link_shape = "logit"))

# item 5 and contract 9
show("mixture(poisson(), categorical())", mixture(poisson(), categorical()))
show("mixture(lognormal, exgaussian, poisson())", mixture(lognormal, exgaussian, poisson()))
show("mixture(poisson, 'cumulative')", mixture(poisson, "cumulative"))
show("mixture(gaussian, hurdle_gamma)", mixture(gaussian, hurdle_gamma))
show("mixture(poisson, bernoulli)", mixture(poisson, bernoulli))
show("mixture(gaussian, Beta)", mixture(gaussian, Beta))
show("mixture(poisson, multinomial)", mixture(poisson, multinomial))
show("mixture(poisson, cox)", mixture(poisson, cox))

# item 1: the fields
str(unclass(student()), max.level = 1)
str(unclass(cumulative("probit")), max.level = 1)
str(unclass(categorical()), max.level = 1)
str(unclass(mixture(gaussian, student)), max.level = 1)
str(unclass(von_mises()), max.level = 1)
print(weibull())
print(student(), links = TRUE)
