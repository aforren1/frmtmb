# Review of lane alapse, check 2: every reader of the covariance on a
# refused fit; zi's block against glm(); a fit with random effects is
# not refused; check_se ignore/stop.
#   Rscript dev/reviews/alapse-review-readers.R <base|lane>
lib_arg <- commandArgs(TRUE)[1]
source("dev/reviews/alapse-review-common.R")
set.seed(5)
d <- data.frame(x = stats::rnorm(100), z = stats::rnorm(100),
                g = factor(sample(letters[1:3], 100, TRUE)))
d$y <- 1 + 0.5 * d$x + stats::rnorm(100)

r <- cap(frm(y ~ x + g, family = asym_laplace(), data = d))
say("frm", r)
f <- r$v
cat("lost:", paste(names(ns$sdr_of(f)$se_lost), ns$sdr_of(f)$se_lost,
                   sep = "=", collapse = " "), "\n")
show <- function(lab, expr) {
  x <- cap(expr)
  say(lab, x)
  if (!inherits(x$v, "error")) {
    out <- utils::capture.output(print(x$v))
    cat(paste0("  | ", head(out, 25)), sep = "\n")
  }
  invisible(x$v)
}
show("summary", summary(f))
show("vcov", vcov(f))
show("confint", confint(f))
show("confint-profile", confint(f, parm = "x", method = "profile"))
show("fixef", fixef(f))
show("predict-se", predict(f, newdata = d[1:3, ], se.fit = TRUE))
show("fitted", head(fitted(f), 3))
show("hypothesis", hypothesis(f, "x = 0"))
show("diagnose", diagnose(f))
if (requireNamespace("emmeans", quietly = TRUE)) {
  show("emmeans", emmeans::emmeans(f, ~ g))
}
if (requireNamespace("marginaleffects", quietly = TRUE)) {
  show("avg_slopes", marginaleffects::avg_slopes(f, variables = "x"))
  show("predictions", head(marginaleffects::predictions(f), 3))
}
if (requireNamespace("frmtmb.sample", quietly = TRUE)) {
  # the start covariance frm_sample() reads (sample.R: sqrt(diag(sdr_of(fit)$cov.fixed)))
  cat("frm_sample start se:", format(sqrt(diag(ns$sdr_of(f)$cov.fixed))), "\n")
  s <- cap(frmtmb.sample::frm_sample(f, chains = 1, iter = 300, warmup = 150,
                                     refresh = 0, seed = 1))
  say("frm_sample", s)
  if (!inherits(s$v, "error")) {
    show("frm_sample summary", summary(s$v)$fixed)
  }
}

# zi: its block against glm on the zero indicator, with a factor
set.seed(9)
dz <- data.frame(x = stats::rnorm(300),
                 h = factor(sample(c("u", "v", "w"), 300, TRUE)))
dz$y <- ifelse(stats::runif(300) < stats::plogis(-1.5 + dz$x + (dz$h == "v")),
               0, brms::rasym_laplace(300, 1 + 0.5 * dz$x, 1, 0.5))
rz <- cap(frm(bf(y ~ x, zi ~ x + h), family = zero_inflated_asym_laplace(),
              data = dz))
say("frm zi", rz)
fz <- rz$v
se <- sqrt(diag(ns$sdr_of(fz)$cov.fixed))
nm <- ns$outer_par_names(fz)
zi <- grepl("^zi_", nm)
g <- stats::glm(I(y == 0) ~ x + h, family = stats::binomial(), data = dz)
cat("zi est frm:", format(fz$opt$par[zi]), "\n")
cat("zi est glm:", format(coef(g)), "\n")
cat("zi se ratio frm/glm:", format(se[zi] / sqrt(diag(vcov(g)))), "\n")
cat("lost:", paste(names(ns$sdr_of(fz)$se_lost), collapse = " "), "\n")

# random effects, not refused (outside mu, censored, log link)
d$gg <- factor(rep(1:10, 10))
for (k in 1:3) {
  fo <- switch(k, bf(y ~ x, sigma ~ (1 | gg)), NULL, NULL)
  rr <- if (k == 1) {
    cap(frm(fo, family = asym_laplace(), data = d))
  } else if (k == 2) {
    d$cy <- ifelse(d$y > 2, "right", "none")
    d$yc <- pmin(d$y, 2)
    cap(frm(yc | cens(cy) ~ x + (1 | gg), family = asym_laplace(), data = d))
  } else {
    d$yp <- exp(d$y / 2)
    cap(frm(yp ~ x + (1 | gg), family = asym_laplace(link = "log"), data = d))
  }
  cat("== RE case", k, "\n")
  say("frm", rr)
  if (!inherits(rr$v, "error")) {
    cat("  ald_cert null:", is.null(rr$v$opt$ald_cert), " lost:",
        paste(names(ns$sdr_of(rr$v)$se_lost), ns$sdr_of(rr$v)$se_lost,
              sep = "=", collapse = " "), "\n")
  }
}

# check_se semantics
for (act in c("ignore", "stop", "warn")) {
  rc <- cap(frm(y ~ x, family = asym_laplace(), data = d,
                control = frmtmb_control(check_se = act)))
  cat("== check_se", act, "\n")
  say("frm", rc)
  if (!inherits(rc$v, "error")) {
    show("summary", summary(rc$v)$fixed)
  }
}
cat("DONE\n")
