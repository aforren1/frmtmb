print(brms:::link)
cat("=====inv_link=====\n")
print(brms:::inv_link)
cat("=====probit approx stan=====\n")
sc <- brms::stancode(brms::bf(y ~ x), data = data.frame(y = rbinom(20,1,.5), x = rnorm(20)),
                     family = brms::bernoulli(link = "probit_approx"))
cat(as.character(sc))
