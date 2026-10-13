# Reviewer of lane optrap (re-check copy of dev/optrap-rev-scale.R, output to dev/optrap-rev2-log/), claim 1: families, links and designs the
# lane's scale study does not cover. One case per process; every fit at
# every scale on one build. A fit's log likelihood is unscaled to the
# response at scale 1 (Jacobian `jac`), and an external fit at the same
# scale is recorded where one fits the same model class.
#   Rscript dev/optrap-rev-scale.R <lib or ""> <case> <tag> [seeds]
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
case <- a[2]
tag <- a[3]
seeds <- if (length(a) > 3) eval(parse(text = a[4])) else 1:20
suppressMessages({library(frmtmb); library(lme4); library(glmmTMB)})
cat("lib", dirname(find.package("frmtmb")), "case", case, "\n")
# OPTRAP_REV_GU=full: the convergence check's trip-wire in the full
# par_units instead of grad_units() (dev/optrap-rev-gu.R)
if (identical(Sys.getenv("OPTRAP_REV_GU"), "full")) {
  utils::assignInNamespace("grad_units", function(fit) {
    fit$par_units %||% 1
  }, "frmtmb")
  cat("grad_units patched: full par_units\n")
}
q <- function(expr) suppressMessages(suppressWarnings(expr))
scales <- 10^c(-6, -3, 0, 3, 6)
if (case == "bigmean") scales <- c(0, 1e3, 1e6, 1e9)
if (case == "mv") scales <- 1:5
mv_k <- list(c(1, 1), c(1e6, 1e-3), c(1e-3, 1e6), c(1e3, 1e3), c(1e6, 1e6))
gen <- function(s) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  b <- rnorm(20)[d$g]
  switch(case,
    gauss_log = {
      d$y <- exp(1.5 + 0.3 * d$x + 0.3 * b) + rnorm(100, 0, 0.5)
    },
    lognormal = {
      d$y <- exp(0.5 + 0.5 * d$x + 0.5 * b + rnorm(100, 0, 0.5))
    },
    gamma = {
      mu <- exp(1 + 0.3 * d$x + 0.3 * b)
      d$y <- rgamma(100, shape = 3, rate = 3 / mu)
    },
    student_sigx = {
      d$y <- 1 + 0.5 * d$x + 0.5 * b + exp(0.3 * d$x) * rt(100, 5)
    },
    gauss_sigx = {
      d$y <- 1 + 0.5 * d$x + 0.5 * b + exp(0.3 * d$x) * rnorm(100)
    },
    skewn = {
      d$y <- 1 + 0.5 * d$x + 0.5 * b + (rexp(100) - 1)
    },
    weights = {
      d$y <- 1 + 0.5 * d$x + 0.5 * b + rnorm(100)
      d$w <- round(runif(100, 0.5, 2), 2)
    },
    offset = {
      d$o <- 2 * rnorm(100)
      d$y <- 1 + 0.5 * d$x + d$o + 0.5 * b + rnorm(100)
    },
    mixture = {
      z <- rbinom(100, 1, 0.4)
      d$y <- 0.5 * d$x + ifelse(z == 1, rnorm(100, 4, 1), rnorm(100, 0, 1))
    },
    bigmean = {
      d$y <- 1 + 0.5 * d$x + 0.5 * b + rnorm(100)
    },
    mv = {
      d$y1 <- 1 + 0.5 * d$x + 0.5 * b + rnorm(100)
      d$y2 <- 2 - 0.3 * d$x + 0.5 * rnorm(20)[d$g] + rnorm(100)
    },
    exgauss = {
      d$y <- 1 + 0.5 * d$x + 0.5 * b + rnorm(100) + rexp(100)
    },
    asym_laplace = {
      d$y <- 1 + 0.5 * d$x + 0.5 * b + rexp(100) - rexp(100)
    },
    huber = {
      d$y <- 1 + 0.5 * d$x + 0.5 * b + rt(100, 3)
    },
    shifted_lognormal = {
      d$y <- 2 + exp(0.3 * d$x + 0.3 * b + rnorm(100, 0, 0.5))
    },
    hurdle_lognormal = {
      d$y <- ifelse(rbinom(100, 1, 0.3) == 1, 0,
                    exp(0.5 + 0.5 * d$x + 0.5 * b + rnorm(100, 0, 0.5)))
    })
  d
}
fitit <- function(d) {
  switch(case,
    gauss_log = frm(bf(y ~ x + (1 | g)), family = gaussian(link = "log"),
                    data = d),
    lognormal = frm(bf(y ~ x + (1 | g)), family = frmtmb::lognormal(), data = d),
    gamma = frm(bf(y ~ x + (1 | g)), family = Gamma(link = "log"),
                data = d),
    student_sigx = frm(bf(y ~ x + (1 | g), sigma ~ x), family = student(),
                       data = d),
    gauss_sigx = frm(bf(y ~ x + (1 | g), sigma ~ x), family = gaussian(),
                     data = d),
    skewn = frm(bf(y ~ x + (1 | g)), family = skew_normal(), data = d),
    weights = frm(bf(y | weights(w) ~ x + (1 | g)), family = gaussian(),
                  data = d),
    offset = frm(bf(y ~ x + offset(o) + (1 | g)), family = gaussian(),
                 data = d),
    mixture = frm(bf(y ~ x), family = mixture(gaussian, gaussian),
                  data = d),
    bigmean = frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d),
    mv = frm(mvbf(bf(y1 ~ x + (1 | g)), bf(y2 ~ x + (1 | g))) +
               set_rescor(FALSE), data = d),
    exgauss = frm(bf(y ~ x + (1 | g)), family = exgaussian(), data = d),
    asym_laplace = frm(bf(y ~ x + (1 | g)), family = asym_laplace(),
                       data = d),
    huber = frm(bf(y ~ x + (1 | g)), family = huber(), data = d),
    shifted_lognormal = frm(bf(y ~ x + (1 | g)),
                            family = shifted_lognormal(), data = d),
    hurdle_lognormal = frm(bf(y ~ x + (1 | g)),
                           family = hurdle_lognormal(), data = d))
}
ext <- function(d) {
  ll <- function(m) as.numeric(logLik(m))
  r <- tryCatch(switch(case,
    gauss_log = ll(glmmTMB(y ~ x + (1 | g), family = gaussian(link = "log"),
                           data = d)),
    lognormal = ll(lmer(log(y) ~ x + (1 | g), data = d, REML = FALSE)) -
      sum(log(d$y)),
    gamma = ll(glmmTMB(y ~ x + (1 | g), family = Gamma(link = "log"),
                       data = d)),
    student_sigx = ll(glmmTMB(y ~ x + (1 | g), dispformula = ~x,
                              family = t_family(), data = d)),
    gauss_sigx = ll(glmmTMB(y ~ x + (1 | g), dispformula = ~x,
                            family = gaussian(), data = d)),
    skewn = ll(glmmTMB(y ~ x + (1 | g), family = skewnormal(), data = d)),
    offset = ll(lmer(y ~ x + offset(o) + (1 | g), data = d, REML = FALSE)),
    bigmean = ll(lmer(y ~ x + (1 | g), data = d, REML = FALSE)),
    mv = ll(lmer(y1 ~ x + (1 | g), data = d, REML = FALSE)) +
      ll(lmer(y2 ~ x + (1 | g), data = d, REML = FALSE)),
    NA_real_), error = function(e) NA_real_)
  q(r)
}
X <- list()
for (s in seeds) {
  d0 <- gen(s)
  for (sc in scales) {
    d <- d0
    jac <- 0
    if (case == "mv") {
      k <- mv_k[[sc]]
      d$y1 <- d$y1 * k[1]
      d$y2 <- d$y2 * k[2]
      jac <- nrow(d) * sum(log(k))
    } else if (case == "bigmean") {
      d$y <- d$y + sc
    } else {
      d$y <- d$y * sc
      if (case == "offset") d$o <- d$o * sc
      jac <- if (case == "weights") sum(d$w) * log(sc) else
        if (case == "hurdle_lognormal") sum(d$y > 0) * log(sc) else
          nrow(d) * log(sc)
    }
    w <- character()
    t0 <- proc.time()[[3]]
    f <- tryCatch(withCallingHandlers(suppressMessages(fitit(d)),
      warning = function(x) {
        w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")
      }), error = function(e) e)
    el <- proc.time()[[3]] - t0
    lx <- ext(d)
    if (inherits(f, "error")) {
      row <- data.frame(case = case, seed = s, scale = sc, code = NA,
                        ll = NA, llu = NA, ext = lx + jac, nwarn = NA,
                        units = NA, secs = el,
                        warn = paste("ERROR", substr(conditionMessage(f),
                                                     1, 80)))
    } else {
      ll <- as.numeric(logLik(f))
      u <- tryCatch(frmtmb:::predictor_units(f$frame),
                    error = function(e) NULL)
      row <- data.frame(case = case, seed = s, scale = sc,
                        code = f$opt$convergence, ll = ll, llu = ll + jac,
                        ext = lx + jac, nwarn = length(w),
                        units = if (is.null(u) || !length(u)) "1" else
                          paste(signif(u, 3), collapse = ","),
                        secs = el,
                        warn = gsub("[\t\n]", " ",
                                    paste(unique(substr(w, 1, 70)),
                                          collapse = " || ")))
    }
    X[[length(X) + 1L]] <- row
  }
}
X <- do.call(rbind, X)
dir.create("dev/optrap-rev2-log", showWarnings = FALSE)
write.table(X, sprintf("dev/optrap-rev2-log/scale-%s-%s.tsv", case, tag),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("done", nrow(X), "\n")
