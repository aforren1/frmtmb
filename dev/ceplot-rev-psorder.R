# Reviewer check (lane ceplot): posterior_samples(pars = ) column order
# against brms 2.23.0 on six model shapes (names only: brms fitted at
# fixed_param, one chain of one iteration), and how many warnings
# parnames() gives on a fit with brms loaded.
#   Rscript dev/ceplot-rev-psorder.R > dev/ceplot-rev-log/psorder.txt
# Data seed 11.
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(brms)
  library(frmtmb)
  library(frmtmb.sample)
})
source("C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-rev-shapes.R")
qb <- function(e) suppressMessages(suppressWarnings(e))
set.seed(11)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:8, 25)),
                o = sample(1:4, n, TRUE))
d$y <- rnorm(n, 1 + 0.5 * d$x + rnorm(8)[d$g], exp(0.1 * d$z))
d$y2 <- rnorm(n, 0.4 * d$z)
d$yp <- 2 * exp(0.3 * d$x) + 0.2 * d$z + rnorm(n, 0, 0.3)
d$yc <- rpois(n, exp(0.2 + 0.3 * d$x))
d$yo <- factor(cut(d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf)), ordered = TRUE)
cases <- list(
  sigma_re = list(frmtmb::bf(y ~ x + (1 | g), sigma ~ z), stats::gaussian(),
                  brms::bf(y ~ x + (1 | g), sigma ~ z), stats::gaussian()),
  mv = list(frmtmb::bf(y ~ x) + frmtmb::bf(y2 ~ z), stats::gaussian(),
            brms::bf(y ~ x) + brms::bf(y2 ~ z) + brms::set_rescor(FALSE),
            stats::gaussian()),
  cum = list(frmtmb::bf(yo ~ x + z), frmtmb::cumulative(),
             brms::bf(yo ~ x + z), brms::cumulative()),
  nl = list(frmtmb::bf(yp ~ a * exp(b * x), a ~ 1 + z, b ~ 1, nl = TRUE),
            stats::gaussian(),
            brms::bf(yp ~ a * exp(b * x), a ~ 1 + z, b ~ 1, nl = TRUE),
            stats::gaussian()),
  nl_sigma = list(frmtmb::bf(yp ~ a * exp(b * x), a ~ 1 + z, b ~ 1, sigma ~ x,
                            nl = TRUE), stats::gaussian(),
                  brms::bf(yp ~ a * exp(b * x), a ~ 1 + z, b ~ 1, sigma ~ x,
                           nl = TRUE), stats::gaussian()),
  zi = list(frmtmb::bf(yc ~ x, zi ~ z), frmtmb::zero_inflated_poisson(),
            brms::bf(yc ~ x, zi ~ z), brms::zero_inflated_poisson()),
  mo = list(frmtmb::bf(y ~ mo(o) + x), stats::gaussian(),
            brms::bf(y ~ mo(o) + x), stats::gaussian()))
pats <- list("^b_", "Intercept", "sigma|zi", c("^sd_", "^b_"), "x",
             c("^b_x$", "^b_Intercept$"))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fit <- frm(cs[[1]], family = cs[[2]], data = d)
  ds <- hand(fit, 4)
  b <- qb(brm(cs[[3]], family = cs[[4]], data = d, algorithm = "fixed_param",
              chains = 1, iter = 1, warmup = 0, refresh = 0, seed = 1,
              silent = 2, init = 0))
  cat("==", nm, "\n")
  for (p in pats) {
    a <- names(qb(posterior_samples(b, pars = p)))
    f <- names(qb(posterior_samples(ds, pars = p)))
    # brms has lprior/lp__ and frmtmb theta_ columns; compare on the
    # names both packages have
    common <- intersect(a, f)
    ok <- identical(a[a %in% common], f[f %in% common])
    cat(sprintf("  pars %-22s same order on %2d common names: %s%s\n",
                paste(p, collapse = ","), length(common), ok,
                if (ok) "" else paste0("\n     brms:   ",
                                       paste(a[a %in% common], collapse = " "),
                                       "\n     frmtmb: ",
                                       paste(f[f %in% common], collapse = " "))))
    only <- setdiff(union(a, f), common)
    if (length(only)) cat("     in one package only:", only, "\n")
  }
}
fit <- frm(frmtmb::bf(y ~ x + (1 | g)), family = stats::gaussian(), data = d)
nw <- 0L
withCallingHandlers(invisible(parnames(fit)), warning = function(w) {
  nw <<- nw + 1L
  invokeRestart("muffleWarning")
})
nb <- 0L
withCallingHandlers(invisible(parnames(b)), warning = function(w) {
  nb <<- nb + 1L
  invokeRestart("muffleWarning")
})
cat("parnames() warnings with brms attached: frmtmb fit", nw, "| brmsfit", nb,
    "\n")
cat("done\n")
