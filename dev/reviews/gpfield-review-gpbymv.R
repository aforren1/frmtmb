# Review of lane gpfield: which warnings test-gp-by.R's mvbf fit (fm) and
# sigma fit (fs) give on the lane, against the `lost` set the test allows
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs(commandArgs(TRUE)[1])
set.seed(5)
n <- 60
d <- data.frame(x = round(stats::runif(n, 0, 6), 1),
                f = factor(rep(c("a", "b", "c"), length.out = n)),
                w = stats::runif(n, 0.5, 2))
d$y <- 0.5 + ifelse(d$f == "a", sin(d$x),
                    ifelse(d$f == "b", cos(d$x), 0.2 * d$x)) +
  stats::rnorm(n, 0, 0.3)
set.seed(31)
d$y2 <- 1 + cos(d$x) * (d$f == "a") + stats::rnorm(nrow(d), 0, 0.4)
d$ysig <- 0.5 + sin(d$x) +
  stats::rnorm(nrow(d), 0, exp(-1 + 0.3 * (d$f == "b") * sin(d$x)))
q <- gf_fit(bf(ysig ~ x, sigma ~ gp(x, by = f, k = 6)), data = d)
cat("FS warn:", unique(sub(":.*", "", q$warn)), "|", substr(q$warn, 1, 60), "\n")
q <- gf_fit(mvbf(bf(y ~ gp(x, by = f, k = 8)), bf(y2 ~ gp(x, k = 6))) +
              set_rescor(FALSE), data = d, family = gaussian())
cat("FM warn:", unique(sub(":.*", "", q$warn)), "|", substr(q$warn, 1, 60), "\n")
