# Review of lane priors2, item 1: fits with no bound at a zero-density
# end, bitwise against 0.82.0. One line per fit: objective, logLik and
# every outer parameter in hex.
# Rscript priors2-review-bitwise.R <lane|base> <out.tsv>
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-priors2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(7)
n <- 240
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:20, 12)))
d$y <- 1 + 0.5 * d$x + rnorm(20)[d$g] * 0.6 + rnorm(n, sd = 0.8)
d$cnt <- rpois(n, exp(0.2 + 0.4 * d$x + rnorm(20, sd = 0.4)[d$g]))
k <- rbinom(n, 1, 0.4)
d$ym <- ifelse(k == 1, rnorm(n, -1 + 0.3 * d$x), rnorm(n, 2.5, 0.7))
d$ym3 <- ifelse(k == 1, rnorm(n, -2), ifelse(runif(n) < 0.5,
                                             rnorm(n, 1.5), rnorm(n, 5)))
d$o <- cut(d$y + rlogis(n), c(-Inf, 0, 1, 2, Inf), labels = FALSE)
d$o <- factor(d$o, ordered = TRUE)
d$t <- rep(1:12, each = 20)
d$ynl <- 3 * exp(-0.4 * abs(d$x) * 2) + rnorm(n, sd = 0.2)
d$xa <- abs(d$x) * 2
cases <- list(
  gauss = function() frm(y ~ x + z, data = d),
  glmm = function() frm(y ~ x + (1 + x | g), data = d),
  pois_glmm = function() frm(cnt ~ x + (1 | g), data = d,
                             family = poisson()),
  nb = function() frm(cnt ~ x, data = d, family = negbinomial()),
  mix2 = function() frm(ym ~ x, data = d,
                        family = mixture(gaussian(), gaussian())),
  mix3 = function() frm(ym3 ~ 1, data = d,
                        family = mixture(gaussian(), gaussian(),
                                         gaussian())),
  ordinal = function() frm(o ~ x + z, data = d, family = cumulative()),
  nonlinear = function() frm(bf(ynl ~ a * exp(-b * xa), a ~ 1, b ~ 1,
                                nl = TRUE), data = d,
                             prior = set_prior("normal(1, 2)", nlpar = "a") +
                               set_prior("normal(0.5, 1)", nlpar = "b")),
  dpar = function() frm(bf(y ~ x, sigma ~ z), data = d),
  ar = function() frm(y ~ x + ar(t, g), data = d),
  prior_b_sd = function() frm(y ~ x + (1 | g), data = d,
                              prior = set_prior("normal(0, 1)", class = "b") +
                                set_prior("student_t(3, 0, 2)",
                                          class = "sd")),
  # bounds that are not at a zero-density end, through a link
  sigma_lb = function() frm(y ~ x, data = d,
                            prior = set_prior("", class = "sigma",
                                              lb = 0.3, ub = 3)),
  sd_ub = function() frm(y ~ x + (1 | g), data = d,
                         prior = set_prior("", class = "sd", ub = 3)),
  b_box = function() frm(y ~ x, data = d,
                         prior = set_prior("normal(0, 1)", class = "b",
                                           lb = -2, ub = 2)),
  mix2_theta_flat = function() {
    f <- frm(ym ~ x, data = d, family = mixture(gaussian(), gaussian()))
    frm(ym ~ x, data = d, family = mixture(gaussian(), gaussian()),
        prior = default_prior(f))
  }
)
hx <- function(v) paste(sprintf("%a", as.numeric(v)), collapse = ",")
out <- character(0)
for (nm in names(cases)) {
  r <- tryCatch(suppressWarnings(suppressMessages(cases[[nm]]())),
                error = function(e) e)
  line <- if (inherits(r, "error")) {
    paste(nm, "ERROR", gsub("[\t\n]", " ", conditionMessage(r)), sep = "\t")
  } else {
    paste(nm, r$opt$convergence, hx(r$opt$objective),
          hx(as.numeric(logLik(r))), hx(r$opt$par), sep = "\t")
  }
  cat(substr(line, 1, 120), "\n")
  out <- c(out, line)
}
writeLines(out, args[2])
