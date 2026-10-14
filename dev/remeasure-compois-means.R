# The means behind R6's upstream entry: brms 2.23.0's Stan-side
# com_poisson pmf (proportional to mu^y / (y!)^nu, mu the rate) and its
# R-side brms:::dcom_poisson() (proportional to (mu^y / y!)^shape), both
# at mu = 3 and nu = shape = 2, summed exactly over 0:200.
#
#   Rscript dev/remeasure-compois-means.R
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
y <- 0:200
mean_of <- function(lp) {
  p <- exp(lp - max(lp))
  p <- p / sum(p)
  sum(y * p)
}
stan_lp <- y * log(3) - 2 * lgamma(y + 1)
r_lp <- 2 * (y * log(3) - lgamma(y + 1))
cat(sprintf("Stan side (rate 3, nu 2): mean %.6f\n", mean_of(stan_lp)))
cat(sprintf("R side (mu 3, shape 2), normalized exactly: mean %.6f\n",
            mean_of(r_lp)))
br <- exp(brms:::dcom_poisson(y, 3, 2, log = TRUE))
cat(sprintf("brms:::dcom_poisson(mu 3, shape 2): mass %.6f, mean %.6f\n",
            sum(br), sum(y * br) / sum(br)))
