# Reviewer of lane defects, recheck of B2: data and models shared by the
# brms and frmtmb arms (seed 5, n = 150, the construction of
# dev/defects-rev-probe5.R plus a second ordinal response).
mvcs_data <- function() {
  set.seed(5)
  n <- 150
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  d$yo <- cut(d$x + rnorm(n), c(-Inf, -0.5, 0.5, Inf), labels = FALSE)
  d$yg <- d$z + rnorm(n)
  d$yo2 <- cut(d$z + 0.5 * d$x + rnorm(n), c(-Inf, -1, 0, 1, Inf),
               labels = FALSE)
  d
}
mvcs_models <- list(
  one_cs = quote(bf(yg ~ x) + bf(yo ~ x + cs(z), family = sratio()) +
                   set_rescor(FALSE)),
  two_cs = quote(bf(yo ~ x + cs(z), family = sratio()) +
                   bf(yo2 ~ z + cs(x), family = acat()) + set_rescor(FALSE)),
  no_ord = quote(bf(yg ~ x) + bf(yo ~ x, family = sratio()) +
                   set_rescor(FALSE))
)
