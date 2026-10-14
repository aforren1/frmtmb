# Fix round 2: does the location-constant refusal reach lca()? Its
# locations are the class predictors theta1 ... theta<K-1>.
#   Rscript dev/famsimplex-lca.R lane|base  -> dev/famsimplex-lca-<arm>.log
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (identical(arm, "lane")) "C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.latent)})
cat("arm:", arm, " lib:", find.package("frmtmb"), "\n")
set.seed(1)
n <- 300
cl <- rbinom(n, 1, 0.4) + 1
pr <- rbind(c(0.85, 0.80, 0.75, 0.90), c(0.15, 0.20, 0.25, 0.10),
            c(0.5, 0.5, 0.5, 0.5))
cl3 <- sample(1:3, n, TRUE)
Y <- matrix(0L, n, 4)
for (j in 1:4) Y[, j] <- 1L + rbinom(n, 1, pr[cl3, j])
dd <- data.frame(x = rnorm(n), o = 0.3)
dd$Y <- Y
run <- function(lab, expr) {
  r <- tryCatch(suppressMessages(suppressWarnings(expr)),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(lab, ":", if (inherits(r, "frmtmb_fit")) {
    sprintf("logLik %.10f; par %s", as.numeric(logLik(r)),
            paste(sprintf("%.6f", r$opt$par), collapse = " "))
  } else substr(r, 1, 150), "\n")
}
fam <- lca(K = 3)
cat("primary dpars:", fam$primary_dpars, "\n")
run("lca(3) free", frm(bf(Y ~ 1), family = fam, data = dd))
run("lca(3) theta1 = 0", frm(bf(Y ~ 1, theta1 = 0), family = fam, data = dd))
run("lca(3) theta2 = 0.5", frm(bf(Y ~ 1, theta2 = 0.5), family = fam,
                               data = dd))
run("lca(3) theta2 ~ 0 + offset(o)",
    frm(bf(Y ~ 1, theta2 ~ 0 + offset(o)), family = fam, data = dd))
