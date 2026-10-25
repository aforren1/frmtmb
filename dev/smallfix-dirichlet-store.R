# How set_prior() stores the two spellings of a distributional
# parameter's intercept, and brms 2.23.0 on class = "sigma" with sigma
# not modeled (validate_prior(), stancode(), brm() with the compile).
#   Rscript dev/smallfix-dirichlet-store.R <lib;lib>
# Log: dev/smallfix-out/dirichlet-store-<arm>.txt
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(a[1], ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
p1 <- set_prior("normal(0, 1)", class = "sigma")
p2 <- set_prior("normal(0, 1)", class = "Intercept", dpar = "sigma")
cat("stored identically (normal):", identical(unclass(p1), unclass(p2)), "\n")
str(unclass(p1)[[1]][c("class", "dpar", "coef")])
set.seed(1)
d <- data.frame(x = rnorm(40))
d$y <- 1 + d$x + rnorm(40)
bp <- brms::set_prior("dirichlet(2)", class = "sigma")
step <- function(lab, expr) {
  r <- tryCatch({ force(expr); "accepted" }, error = function(e) {
    paste("REFUSED:", substr(gsub("\\s+", " ", conditionMessage(e)), 1,
                             200))
  })
  cat(sprintf("  %-22s %s\n", lab, r))
}
cat("brms, class = \"sigma\", y ~ x (sigma not modeled)\n")
step("validate_prior", brms::validate_prior(bp, brms::bf(y ~ x), data = d))
step("stancode", sc <- brms::stancode(brms::bf(y ~ x), data = d, prior = bp))
cat(grep("dirichlet", sc, value = TRUE), "\n")
step("brm (compile)", suppressWarnings(suppressMessages(
  brms::brm(brms::bf(y ~ x), data = d, prior = bp, chains = 0, refresh = 0,
            silent = 2, backend = "rstan"))))
cat("frmtmb, class = \"sigma\", y ~ x\n")
step("set_prior", fp <- set_prior("dirichlet(2)", class = "sigma"))
step("validate_prior", validate_prior(fp, bf(y ~ x), data = d))
step("frm", frm(bf(y ~ x), family = gaussian(), data = d, prior = fp))
