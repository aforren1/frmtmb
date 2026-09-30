# Reviewer of lane defects, recheck: the data shared by
# dev/defects-rev-mi-brms.R and dev/defects-rev-mi-frmtmb.R.
mi_data <- function() {
  set.seed(21)
  n <- 80
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  d$y <- 1 + d$x + rnorm(n)
  d$ymi <- ifelse(runif(n) < 0.2, NA, d$y)
  d$xmi <- ifelse(runif(n) < 0.2, NA, d$x)
  d$y2 <- 0.5 * d$y + d$z + rnorm(n)
  # mi(sd = ): measured with known error, and some rows unmeasured
  d$sdy <- 0.3
  d$ymeas <- d$y + rnorm(n, 0, 0.3)
  d$ymeas[c(3, 9, 27, 41, 66)] <- NA
  d
}
mi_models <- list(
  uni = quote(bf(ymi | mi() ~ x)),
  mv = quote(bf(y ~ mi(xmi) + z) + bf(xmi | mi() ~ z) + set_rescor(FALSE)),
  mv2 = quote(bf(ymi | mi() ~ x) + bf(y2 ~ z) + set_rescor(FALSE)),
  sdy = quote(bf(ymeas | mi(sdy) ~ x))
)
