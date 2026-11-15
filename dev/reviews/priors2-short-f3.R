# Short check of lane priors2's fix round, item 3 (F3): bound_inside()
# and internal_bound() end, land inside, and refuse what they must.
# Rscript priors2-short-f3.R <lane|base>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-priors2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
ib <- frmtmb:::internal_bound
if (arm == "lane") {
  links <- list(logit = make.link("logit"), probit = make.link("probit"),
                cloglog = make.link("cloglog"), log = make.link("log"))
  ncall <- 0L
  count_link <- function(lk) {
    f <- lk$linkinv
    lk$linkinv <- function(x) {
      ncall <<- ncall + 1L
      f(x)
    }
    lk
  }
  run <- function(name, vs) {
    lk <- count_link(links[[name]])
    pm <- list(scale = "link", link = lk)
    out <- 0L
    maxcall <- 0L
    t0 <- proc.time()[["elapsed"]]
    for (side in c("lb", "ub")) {
      for (v in vs) {
        ncall <<- 0L
        b <- ib(v, pm, side)
        maxcall <- max(maxcall, ncall)
        y <- links[[name]]$linkinv(b)
        ok <- if (side == "lb") y >= v else y <= v
        if (!isTRUE(ok)) out <- out + 1L
      }
    }
    el <- proc.time()[["elapsed"]] - t0
    cat(sprintf("%-8s %-26s ends %5d  outside %d  max linkinv calls %d  %.3f ms/bound\n",
                name, attr(vs, "lab"), 2L * length(vs), out, maxcall,
                1000 * el / (2 * length(vs))))
  }
  lab <- function(v, s) structure(v, lab = s)
  set.seed(3)
  run("logit", lab(1 - 10^-seq(2, 5, by = 0.01), "1 - 10^-seq(2, 5, .01)"))
  run("logit", lab(10^-seq(2, 12, by = 0.01), "10^-seq(2, 12, .01)"))
  # 2^-1070 is refused through make.link("logit"): priors2-short-f3b.R
  run("logit", lab(c(0.99981802991413904, 1 - 1e-15),
                   "review ub, 1 - 1e-15"))
  run("logit", lab(runif(5000), "runif(5000)"))
  run("probit", lab(runif(5000), "runif(5000)"))
  run("probit", lab(1 - 10^-seq(2, 15, by = 0.01), "1 - 10^-seq(2, 15, .01)"))
  run("cloglog", lab(runif(5000), "runif(5000)"))
  run("log", lab(exp(runif(5000, -700, 700)), "exp(U(-700, 700))"))
  run("log", lab(1 + (-2000:2000) * .Machine$double.eps,
                 "1 + k eps, k in -2000:2000"))
  # the sd scale goes through exp()
  pm_sd <- list(scale = "sd")
  v <- exp(runif(5000, -50, 50))
  bad <- sum(exp(vapply(v, ib, 0, pm = pm_sd, which = "lb")) < v) +
    sum(exp(vapply(v, ib, 0, pm = pm_sd, which = "ub")) > v)
  cat("sd scale, 10000 ends, outside", bad, "\n")
  # the helper refuses by name where no value maps inside
  r <- tryCatch(frmtmb:::bound_inside(0, 2, stats::plogis, "lb"),
                error = conditionMessage)
  cat("bound_inside(0, 2, plogis, lb):", r, "\n")
}
# what internal_bound() does with the edge values
pl <- list(scale = "link", link = make.link("logit"))
pg <- list(scale = "link", link = make.link("log"))
edge <- function(lab, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR", conditionMessage(e)))
  cat(sprintf("%-28s %s\n", lab, format(r, digits = 17)))
}
edge("logit ub = 1", ib(1, pl, "ub"))
edge("logit lb = 0", ib(0, pl, "lb"))
edge("logit lb = 1", ib(1, pl, "lb"))
edge("logit ub = 0", ib(0, pl, "ub"))
edge("logit lb = 2", ib(2, pl, "lb"))
edge("logit ub = 2", ib(2, pl, "ub"))
edge("logit lb = -1", ib(-1, pl, "lb"))
edge("logit ub = -1", ib(-1, pl, "ub"))
edge("log lb = 0", ib(0, pg, "lb"))
edge("log ub = 0", ib(0, pg, "ub"))
edge("log lb = -1", ib(-1, pg, "lb"))
edge("log lb = Inf", ib(Inf, pg, "lb"))
# through frm(): the user's path
set.seed(1)
dz <- data.frame(z = rnorm(120))
dz$cnt <- rpois(120, 2) * rbinom(120, 1, 0.5)
dz$y <- 1 + 0.3 * dz$z + rnorm(120)
fit_line <- function(lab, expr) {
  r <- tryCatch(suppressWarnings(suppressMessages(expr)),
                error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-34s ERROR %s\n", lab, substr(conditionMessage(r), 1, 140)))
  } else {
    cat(sprintf("%-34s code %d objective %.15g estimate %s\n", lab,
                r$opt$convergence, r$opt$objective,
                paste(format(r$opt$par, digits = 10), collapse = " ")))
  }
}
fit_line("zi lb = 2 (logit)", frm(cnt ~ z, data = dz,
                                  family = zero_inflated_poisson(),
                                  prior = set_prior("", class = "zi", lb = 2)))
fit_line("zi lb = 1 (logit)", frm(cnt ~ z, data = dz,
                                  family = zero_inflated_poisson(),
                                  prior = set_prior("", class = "zi", lb = 1)))
fit_line("zi ub = 1 (logit)", frm(cnt ~ z, data = dz,
                                  family = zero_inflated_poisson(),
                                  prior = set_prior("", class = "zi", ub = 1)))
fit_line("zi ub = 0.3 (logit)", frm(cnt ~ z, data = dz,
                                    family = zero_inflated_poisson(),
                                    prior = set_prior("", class = "zi",
                                                      ub = 0.3)))
fit_line("zi none", frm(cnt ~ z, data = dz, family = zero_inflated_poisson()))
fit_line("sigma lb = 0 (log)", frm(y ~ z, data = dz,
                                   prior = set_prior("", class = "sigma",
                                                     lb = 0)))
fit_line("sigma none", frm(y ~ z, data = dz))
fit_line("sigma ub = 0 (log)", frm(y ~ z, data = dz,
                                   prior = set_prior("", class = "sigma",
                                                     ub = 0)))
fit_line("shape lb = 0 (log, negbinomial)",
         frm(cnt ~ z, data = dz, family = negbinomial(),
             prior = set_prior("", class = "shape", lb = 0)))
