.libPaths(c("C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
p <- c(set_prior("student_t(3, 0, 2.5)", class = "sigma"),
       set_prior("normal(0, 1)", class = "Intercept", dpar = "sigma"),
       set_prior("gamma(2, 1)", class = "shape"), set_prior("normal(0,5)", class = "b"))
str(unclass(p)[[1]]); print(as.data.frame(p))
str(frmtmb:::parse_prior_dist("gamma(2, 1)"))
