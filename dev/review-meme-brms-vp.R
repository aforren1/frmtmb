# Review of lane meme: brms 2.23.0's validate_prior() fill on coef rows,
# class b against meanme/sdme, and the corme/cor group rows. No compile.
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
source("C:/Users/adf44/source/r/frmtmb-wt-meme/dev/review-meme-shapes.R")
d <- rv_data()
f <- bf(y ~ me(x, sx) + me(z, sz) + w + (1 + w | g))
pr <- c(set_prior("normal(0, 1)", class = "b"),
        set_prior("normal(0, 2)", class = "meanme"),
        set_prior("exponential(1)", class = "sdme"),
        set_prior("student_t(3, 0, 1)", class = "sd"),
        set_prior("lkj(2)", class = "cor"))
vp <- as.data.frame(validate_prior(pr, f, data = d))
print(vp[, c("prior", "class", "coef", "group", "lb", "ub", "source")])
cat("\n-- default_prior with (1 + w | g) and correlated me():\n")
print(as.data.frame(default_prior(f, data = d))[, c("prior", "class", "coef",
                                                     "group", "source")])
cat("\n-- flat corme: stancode statements\n")
sc <- stancode(bf(y ~ me(x, sx) + me(z, sz) + me(w, sw)), data = d,
               prior = c(set_prior("", class = "corme"),
                         set_prior("", class = "Intercept"),
                         set_prior("", class = "sigma")))
cat(grep("lprior|target", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
cat("\n-- brms: meanme with coef and bound\n")
r <- tryCatch(validate_prior(set_prior("normal(0, 1)", class = "meanme",
                                       coef = "mex", lb = 0),
                             bf(y ~ me(x, sx) + me(z, sz)), data = d),
              error = function(e) conditionMessage(e))
print(r)
r <- tryCatch(validate_prior(set_prior("normal(0, 1)", class = "b",
                                       coef = "w", lb = 0),
                             bf(y ~ me(x, sx) + w), data = d),
              error = function(e) conditionMessage(e))
cat("-- brms: b with coef and bound\n")
print(r)
r <- tryCatch(validate_prior(set_prior("normal(0, 1)", class = "sd",
                                       group = "g", coef = "w", lb = 0),
                             f, data = d),
              error = function(e) conditionMessage(e))
cat("-- brms: sd with coef and bound\n")
print(r)
