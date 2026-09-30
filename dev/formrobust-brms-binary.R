# brms 2.23.0: which families are binary, and standata's Y coding for
# two-valued responses of several types.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
fams <- c("bernoulli", "binomial", "categorical", "cumulative", "gaussian",
          "zero_inflated_binomial", "beta_binomial")
for (f in fams) cat(f, brms:::is_binary(brms:::validate_family(f)), "\n")
show <- function(lab, y) {
  r <- tryCatch(as.vector(standata(y ~ 1, data = data.frame(y = y),
                                   family = bernoulli())$Y),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-22s %s\n", lab, paste(r, collapse = " ")))
}
show("-1,-2", c(-1, -2, -1, -2))
show("0,1", c(0, 1, 1, 0))
show("1,2", c(1, 2, 2, 1))
show("all 1", c(1, 1, 1))
show("all 0", c(0, 0, 0))
show("all 5", c(5, 5, 5))
show("TRUE/FALSE", c(TRUE, FALSE, TRUE))
show("factor b,a lv a,b", factor(c("b", "a", "b")))
show("factor lv z,a", factor(c("a", "z", "a"), levels = c("z", "a")))
show("character yes/no", c("yes", "no", "yes"))
show("three values", c(0, 1, 2))
show("factor one level", factor(c("a", "a")))
show("factor 3 lv 2 used", factor(c("a", "b"), levels = c("a", "b", "c")))
show("with NA", c(-1, NA, -2))
show("all TRUE", c(TRUE, TRUE, TRUE))
show("all FALSE", c(FALSE, FALSE))
