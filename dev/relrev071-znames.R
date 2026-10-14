.libPaths(c("C:/Users/adf44/source/r/rellib-r9", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
set.seed(11); n <- 100
dd <- data.frame(x = rnorm(n), g = factor(rep(1:10, length.out = n)), f = factor(sample(c("a","b","c"), n, TRUE)))
Y <- matrix(rgamma(3 * n, 2), n); Y <- Y / rowSums(Y); colnames(Y) <- c("a","b","c"); dd$Y <- Y
s <- brms::standata(brms::bf(Y ~ x + (1 + f || g)), data = dd, family = brms::dirichlet())
print(grep("^(Z|J|N|M)_", names(s), value = TRUE))
suppressMessages(library(frmtmb))
sys.source("C:/Users/adf44/source/r/frmtmb-wt-release/tests/testthat/helper-brms.R", envir = globalenv())
s2 <- brms_standata(brms::bf(Y ~ x + (1 + f || g)), data = dd, family = brms::dirichlet())
print(grep("^Z_", names(s2), value = TRUE)); print(class(s2))
nm <- grep(paste0("^Z_", 1, "_"), names(s2), value = TRUE)
print(nm); print(sub(paste0("^Z_", 1, "_(.*_)?([0-9]+)$"), "\2", nm))
