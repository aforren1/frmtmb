# brms 2.23.0: the columns of (1 + f || gr(g, cov = A)) and of
# (1 + x || gr(g, cov = A)), from standata and the prior table
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
ng <- 12
A <- 0.5 ^ abs(outer(1:ng, 1:ng, "-"))
dimnames(A) <- list(paste0("g", 1:ng), paste0("g", 1:ng))
d8 <- data.frame(g = factor(rep(paste0("g", 1:ng), each = 10),
                            levels = paste0("g", 1:ng)),
                 f = factor(rep(c("a", "b", "c"), 40)), x = rnorm(120))
d8$y <- rnorm(120)
for (fo in list(y ~ 1 + (1 + f || gr(g, cov = A)),
                y ~ 1 + (1 + x || gr(g, cov = A)),
                y ~ 1 + (0 + f || gr(g, cov = A)))) {
  cat("\n==", deparse(fo), "\n")
  sd <- standata(fo, data = d8, data2 = list(A = A))
  print(grep("^Z_|^M_|^NC_", names(sd), value = TRUE))
  print(unique(as.data.frame(sd[grep("^Z_", names(sd))])))
  dp <- default_prior(fo, data = d8, data2 = list(A = A))
  print(dp[dp$class == "sd", c("class", "coef", "group")])
}
