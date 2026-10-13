# brms 2.23.0 on simo prior rows: which spellings it takes, and the
# concentration vector each puts in the Stan data. No compile.
# Usage: Rscript dev/moparity-brms-simo.R > dev/moparity-log/brms-simo.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
set.seed(1)
n <- 100
dd <- data.frame(x = sample(0:3, n, TRUE), z = rnorm(n),
                 g = factor(sample(1:5, n, TRUE)))
dd$y <- rnorm(n)
try_prior <- function(lab, fo, pr) {
  cat("\n==", lab, "\n")
  r <- tryCatch({
    sd <- standata(fo, data = dd, prior = pr)
    for (nm in grep("^con_simo", names(sd), value = TRUE)) {
      cat(nm, ":", sd[[nm]], "\n")
    }
    sc <- stancode(fo, data = dd, prior = pr)
    cat(grep("dirichlet", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
  }, error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
}
fo <- bf(y ~ mo(x) * z)
try_prior("coef mox1 dirichlet(c(5,1,1))", fo,
          prior(dirichlet(c(5, 1, 1)), class = simo, coef = mox1))
try_prior("coef mox1 dirichlet(2)", fo,
          prior(dirichlet(2), class = simo, coef = mox1))
try_prior("class-wide dirichlet(2)", fo, prior(dirichlet(2), class = simo))
try_prior("wrong length c(1,1)", fo,
          prior(dirichlet(c(1, 1)), class = simo, coef = mox1))
try_prior("normal on simo", fo, prior(normal(0, 1), class = simo,
                                      coef = mox1))
try_prior("coef mox:z1", fo,
          prior(dirichlet(c(1, 2, 3)), class = simo, coef = "mox:z1"))
try_prior("sigma dpar", bf(y ~ 1, sigma ~ mo(x)),
          prior(dirichlet(c(1, 2, 3)), class = simo, coef = mox1,
                dpar = sigma))
try_prior("lb on simo", fo,
          set_prior("dirichlet(2)", class = "simo", coef = "mox1", lb = 0))
cat("\n== default_prior rows\n")
dp <- default_prior(bf(y ~ mo(x) * z + (mo(x) | g), sigma ~ mo(x)),
                    data = dd)
print(dp[dp$class == "simo", c("prior", "class", "coef", "dpar", "source")])
