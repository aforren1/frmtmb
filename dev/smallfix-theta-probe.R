# brms 2.23.0's class "theta" in a mixture, run (review finding 9):
# what default_prior() lists and what a dirichlet() row on it becomes
# in the Stan code; and frmtmb's refusal of the same row.
#   Rscript dev/smallfix-theta-probe.R <lib;lib>
# Log: dev/smallfix-out/theta-probe-<arm>.txt
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(a[1], ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "| brms",
    format(packageVersion("brms")), "\n")
set.seed(2)
d <- data.frame(y = c(rnorm(60), rnorm(60, 4)))
mix <- brms::mixture(gaussian, gaussian)
dp <- brms::default_prior(brms::bf(y ~ 1), data = d, family = mix)
print(as.data.frame(dp)[as.data.frame(dp)$class == "theta",
                        c("prior", "class", "dpar")])
bp <- brms::set_prior("dirichlet(2, 2)", class = "theta")
sc <- brms::stancode(brms::bf(y ~ 1), data = d, family = mix, prior = bp)
cat(grep("theta|simplex", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
fp <- set_prior("dirichlet(2)", class = "theta")
e <- tryCatch(validate_prior(fp, bf(y ~ 1), data = d,
                             family = frmtmb::mixture(gaussian(),
                                                      gaussian())),
              error = function(e) conditionMessage(e))
cat("frmtmb validate_prior:", gsub("\\s+", " ", e), "\n")
