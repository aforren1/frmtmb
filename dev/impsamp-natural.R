# Does the refusal fire on its own, after the slopes' step is scale-free?
# Very small draw counts on the designs most likely to rise: tiny groups, big
# variance, heavy tails, separation, and a correlated block. Prints each fit
# that fires or caps. Log dev/impsamp-out/natural.txt
.libPaths(c("C:/Users/adf44/source/r/wt-impsamp-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
sq <- function(e) suppressWarnings(suppressMessages(e))
mk <- list(
  bern_sd4 = function(s) {
    set.seed(s); g <- factor(rep(1:40, each = 3L)); x <- rnorm(120)
    list(d = data.frame(y = rbinom(120, 1, plogis(0.5 * x + rnorm(40, 0, 4)[g])),
                        x = x, g = g), f = bf(y ~ x + (1 | g)), fam = bernoulli())
  },
  bern_slope_2 = function(s) {
    set.seed(s); g <- factor(rep(1:30, each = 2L)); x <- rnorm(60)
    list(d = data.frame(y = rbinom(60, 1, 0.4), x = x, g = g),
         f = bf(y ~ x + (x | g)), fam = bernoulli())
  },
  bern_slope_sep = function(s) {
    set.seed(s); g <- factor(rep(1:20, each = 4L)); x <- rnorm(80)
    list(d = data.frame(y = as.integer(x + rnorm(20, 0, 2)[g] + rnorm(80, 0, .1) > 0),
                        x = x, g = g), f = bf(y ~ x + (x | g)), fam = bernoulli())
  },
  pois_tiny = function(s) {
    set.seed(s); g <- factor(rep(1:30, each = 2L)); x <- rnorm(60)
    list(d = data.frame(y = rpois(60, exp(-1 + rnorm(30, 0, 2)[g])), x = x, g = g),
         f = bf(y ~ x + (1 | g)), fam = poisson())
  },
  student_out = function(s) {
    set.seed(s); g <- factor(rep(1:15, each = 3L)); x <- rnorm(45)
    y <- 0.3 + 0.5 * x + rnorm(15, 0, .7)[g] + 0.3 * rt(45, 2)
    y[seq(1, 45, by = 9)] <- y[seq(1, 45, by = 9)] + 6
    list(d = data.frame(y = y, x = x, g = g), f = bf(y ~ x + (1 | g), nu = 3),
         fam = student())
  },
  gamma_two = function(s) {
    set.seed(s); g <- factor(rep(1:20, each = 3L)); x <- rnorm(60)
    list(d = data.frame(y = rgamma(60, 2, 2 / exp(0.2 * x + rnorm(20, 0, 1.5)[g] +
                                                    rnorm(20, 0, 1)[g] * x)),
                        x = x, g = g), f = bf(y ~ x + (1 | g) + (0 + x | g)),
         fam = Gamma(link = "log"))
  },
  nb_mu_sigma = function(s) {
    set.seed(s); g <- factor(rep(1:20, each = 5L)); x <- rnorm(100)
    list(d = data.frame(y = rnbinom(100, mu = exp(0.5 + rnorm(20, 0, 1)[g]),
                                    size = exp(rnorm(20, 0, 1)[g])), x = x, g = g),
         f = bf(y ~ x + (1 | g), shape ~ 1 + (1 | g)), fam = negbinomial())
  })
n <- c(fires = 0L, capped = 0L, fits = 0L, other = 0L)
for (nm in names(mk)) for (s in 1:12) for (N in c(2L, 4L, 8L, 16L, 64L)) {
  m <- mk[[nm]](s)
  r <- tryCatch(sq(frm(m$f, family = m$fam, data = m$d, importance = N)),
                error = function(e) e)
  n[["fits"]] <- n[["fits"]] + 1L
  if (inherits(r, "error")) {
    msg <- conditionMessage(r)
    if (grepl("did not settle", msg, fixed = TRUE)) {
      n[["fires"]] <- n[["fires"]] + 1L
      cat("FIRES", nm, s, N, ":", substr(msg, 1, 400), "\n")
    } else {
      n[["other"]] <- n[["other"]] + 1L
      cat("OTHER", nm, s, N, ":", substr(msg, 1, 160), "\n")
    }
  } else if (isTRUE(r$importance$capped)) {
    n[["capped"]] <- n[["capped"]] + 1L
    cat("CAPPED", nm, s, N, "moves",
        paste(format(r$importance$moves, digits = 3), collapse = ","),
        "ess", format(r$importance$ess_min, digits = 3), "\n")
  }
}
print(n)
