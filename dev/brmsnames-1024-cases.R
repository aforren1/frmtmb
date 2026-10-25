## The formula terms whose coefficient names are compared with brms
## 2.23.0, and the data they are fitted to (seed 1024). Sourced by
## dev/brmsnames-1024-brms.R, dev/brmsnames-1024-frm.R and
## dev/brmsnames-1024-table.R; it attaches no package.
##
## Each case is list(id, form, family, data). `form` is a function of
## the bf()/mvbf() of the package that runs it, so the same text builds
## brms's and frmtmb's formula. `family` is a string, read by each
## package's own family constructor.

brmsnames_data <- function() {
  set.seed(1024)
  n <- 200
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  xp = stats::runif(n, 0.5, 2),
                  k = sample(1:4, n, replace = TRUE),
                  f = factor(sample(c("a", "b", "c"), n, replace = TRUE)),
                  g = factor(rep(1:10, length.out = n)),
                  e = stats::runif(n, 1, 3))
  d$y <- 1 + 0.4 * exp(d$x / 2) + 0.3 * d$z + 0.2 * d$k +
    stats::rnorm(10, 0, 0.3)[d$g] + stats::rnorm(n, 0, exp(0.1 * d$z))
  d$y2 <- 0.5 * d$x + stats::rnorm(n)
  d$xo <- d$x + stats::rnorm(n, 0, 0.3)
  d$sdx <- 0.3
  d$xm <- d$z + stats::rnorm(n, 0, 0.5)
  d$xm[sample(n, 20)] <- NA
  d$yo <- cut(d$x + stats::rlogis(n), c(-Inf, -1, 0, 1, Inf),
              labels = FALSE)
  d$yc <- stats::rpois(n, d$e * exp(0.2 + 0.3 * d$x))
  d$yn <- 2 * exp(0.3 * d$x) + stats::rnorm(n, 0, 0.2)
  d
}

brmsnames_cases <- function(bf, mvbf, set_rescor) {
  g <- "gaussian"
  list(
    list(id = "exp", form = bf(y ~ exp(x)), family = g),
    list(id = "log", form = bf(y ~ log(xp)), family = g),
    list(id = "sqrt", form = bf(y ~ sqrt(xp)), family = g),
    list(id = "Isq", form = bf(y ~ I(x^2)), family = g),
    list(id = "poly", form = bf(y ~ poly(x, 2)), family = g),
    list(id = "polyraw", form = bf(y ~ poly(x, 2, raw = TRUE)), family = g),
    list(id = "scale", form = bf(y ~ scale(x)), family = g),
    list(id = "factor", form = bf(y ~ factor(k)), family = g),
    list(id = "exp:z", form = bf(y ~ exp(x):z), family = g),
    list(id = "exp*f", form = bf(y ~ exp(x) * f), family = g),
    list(id = "log:factor", form = bf(y ~ log(xp):factor(k)), family = g),
    list(id = "poly*z", form = bf(y ~ poly(x, 2) * z), family = g),
    list(id = "I(exp)", form = bf(y ~ I(exp(x))), family = g),
    list(id = "exp+log+sqrt", form = bf(y ~ exp(x) + log(xp) + sqrt(xp)),
         family = g),
    list(id = "exp(x/2)", form = bf(y ~ exp(x / 2)), family = g),
    list(id = "s", form = bf(y ~ s(x)), family = g),
    list(id = "s+exp", form = bf(y ~ s(xp) + exp(z)), family = g),
    list(id = "s(exp)", form = bf(y ~ s(exp(x))), family = g),
    list(id = "t2", form = bf(y ~ t2(x, z)), family = g),
    list(id = "mo", form = bf(y ~ mo(k)), family = g),
    list(id = "mo:z", form = bf(y ~ mo(k):z), family = g),
    list(id = "mo+exp", form = bf(y ~ mo(k) + exp(x)), family = g),
    list(id = "me", form = bf(y ~ me(xo, sdx)), family = g),
    list(id = "me+exp", form = bf(y ~ me(xo, sdx) + exp(z)), family = g),
    list(id = "mi", form = bf(y ~ mi(xm) + exp(x)) +
           bf(xm | mi() ~ z) + set_rescor(FALSE), family = g),
    list(id = "cs", form = bf(yo ~ cs(x) + z), family = "acat"),
    list(id = "cs(exp)", form = bf(yo ~ cs(exp(x)) + z), family = "acat"),
    list(id = "offset", form = bf(yc ~ x + offset(log(e))),
         family = "poisson"),
    list(id = "offset+exp", form = bf(yc ~ exp(x) + offset(log(e))),
         family = "poisson"),
    list(id = "nl", form = bf(yn ~ a * exp(b * x), a ~ 1, b ~ 1,
                              nl = TRUE), family = g),
    list(id = "nl a~exp", form = bf(y ~ a + b * x, a ~ exp(z), b ~ 1,
                                    nl = TRUE), family = g),
    list(id = "sigma~exp", form = bf(y ~ x, sigma ~ exp(z)), family = g),
    list(id = "re exp", form = bf(y ~ exp(x) + (1 + exp(x) | g)),
         family = g),
    list(id = "re log", form = bf(y ~ z + (1 + log(xp) | g)), family = g),
    list(id = "mv exp", form = mvbf(bf(y ~ exp(x)), bf(y2 ~ log(xp)),
                                    rescor = FALSE), family = g),
    list(id = "re exp ||", form = bf(y ~ exp(x) + (1 + exp(x) || g)),
         family = g),
    list(id = "re mo", form = bf(y ~ mo(k) + (1 + mo(k) | g)), family = g),
    list(id = "offset(exp)", form = bf(y ~ x + offset(exp(z / 4))),
         family = g),
    list(id = "nl a~exp sigma~z", form = bf(y ~ a + b * x, a ~ exp(z),
                                            b ~ 1 + z, sigma ~ z,
                                            nl = TRUE), family = g),
    list(id = "mv nl", form = mvbf(bf(yn ~ a * exp(b * x),
                                      a ~ 1 + exp(z) + (1 | g), b ~ 1,
                                      sigma ~ z, nl = TRUE),
                                   bf(y ~ x + exp(z)), rescor = FALSE),
         family = g)
  )
}

## brms's prior spelling of one default_prior() row, as a set_prior()
## argument list: the class, coef, group, dpar, nlpar and resp it names.
brmsnames_prior_args <- function(row) {
  dist <- if (row$class %in% c("sd", "sds", "sigma")) {
    "student_t(3, 0, 5)"
  } else if (identical(row$class, "simo")) {
    "dirichlet(1)"
  } else {
    "normal(0, 5)"
  }
  a <- list(prior = dist, class = row$class, coef = row$coef)
  for (nm in c("group", "dpar", "nlpar", "resp")) {
    if (nzchar(row[[nm]])) a[[nm]] <- row[[nm]]
  }
  a
}
