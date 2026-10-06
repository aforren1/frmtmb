# brms's default_prior() rows for gp() terms against frmtmb's
# get_prior(route = "sample"), and the fit-route table, on the designs
# of dev/gpby-brms-explore.R.
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
set.seed(1)
n <- 40
dd <- data.frame(x = runif(n, 0, 5), z = runif(n, 0, 3),
                 f = factor(sample(c("a", "b", "c"), n, TRUE)),
                 w = runif(n, 0.5, 2))
dd$y <- sin(dd$x) + rnorm(n, 0, 0.3)
forms <- c("y ~ gp(x)", "y ~ gp(x, by = f)", "y ~ gp(x, by = w)",
           "y ~ gp(x, z, by = f)", "y ~ gp(x, z, by = f, iso = FALSE)",
           "y ~ gp(x, by = f, k = 8)", "y ~ gp(x, by = f, cmc = FALSE)",
           "y ~ gp(x, by = f, gr = FALSE)", "y ~ gp(x, z)")
nbad <- 0L
nrow_cmp <- 0L
for (fm in forms) {
  b <- as.data.frame(brms::default_prior(brms::bf(as.formula(fm)),
                                         data = dd))
  b <- b[b$class %in% c("sdgp", "lscale"), c("prior", "class", "coef")]
  g <- as.data.frame(get_prior(frmtmb::bf(as.formula(fm)), data = dd,
                               route = "sample"))
  g <- g[g$class %in% c("sdgp", "lscale"), c("prior", "class", "coef")]
  key <- function(t) paste(t$class, t$coef)
  m <- merge(b, g, by = c("class", "coef"), all = TRUE,
             suffixes = c(".brms", ".frmtmb"))
  # brms writes "(flat)" for the class-wide lscale row; so does frmtmb
  flat <- function(p) ifelse(is.na(p) | !nzchar(p), "(flat)", p)
  m$same <- flat(m$prior.brms) == flat(m$prior.frmtmb)
  nrow_cmp <- nrow_cmp + nrow(m)
  nbad <- nbad + sum(!m$same | is.na(m$same))
  cat("==", fm, "\n")
  print(m, row.names = FALSE)
}
cat(sprintf("ROWS %d compared, %d differ\n", nrow_cmp, nbad))
fit <- frm(frmtmb::bf(y ~ gp(x, by = f)), data = dd)
print(get_prior(fit))
