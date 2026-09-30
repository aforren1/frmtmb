# Reviewer: no-regression. The same calls, none using the new options,
# on the base build and on the lane; each arm saves its outputs to
# dev/postfit2-rev-log/noreg-<arm>.rds, and `compare` diffs them.
#   Rscript dev/postfit2-rev-noreg.R base
#   Rscript dev/postfit2-rev-noreg.R lane
#   Rscript dev/postfit2-rev-noreg.R compare
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "lane"
wt <- "C:/Users/adf44/source/r/frmtmb-wt-postfit2"
out_dir <- file.path(wt, "dev/postfit2-rev-log")
libs <- c("C:/Users/adf44/source/r/wt-postfit2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
say <- function(...) cat(sprintf(...), "\n", sep = "")

if (arm == "compare") {
  .libPaths(libs[-1])
  a <- readRDS(file.path(out_dir, "noreg-base.rds"))
  b <- readRDS(file.path(out_dir, "noreg-lane.rds"))
  strip <- function(x) {
    if (is.data.frame(x)) {
      x <- as.list(x)
      return(lapply(x, function(v) { attributes(v) <- attributes(v)[c("levels", "class")]; v }))
    }
    if (is.list(x)) return(lapply(unclass(x), strip))
    x
  }
  attr_names <- function(x) {
    if (is.data.frame(x)) return(sort(names(attributes(x))))
    if (is.list(x)) return(lapply(unclass(x), attr_names))
    NULL
  }
  relgap <- function(p, q) {
    p <- unlist(p); q <- unlist(q)
    if (!is.numeric(p) || !is.numeric(q) || length(p) != length(q)) return(NA)
    ok <- is.finite(p) & is.finite(q)
    max(abs(p[ok] - q[ok])) / max(abs(p[ok]), 1e-300)
  }
  say("cases: base %d, lane %d", length(a), length(b))
  for (nm in union(names(a), names(b))) {
    x <- a[[nm]]; y <- b[[nm]]
    if (inherits(x, "error") || inherits(y, "error") ||
        is.character(x) || is.character(y)) {
      say("%-34s base: %s | lane: %s", nm,
          if (is.character(x)) substr(x, 1, 90) else "value",
          if (is.character(y)) substr(y, 1, 90) else "value")
      next
    }
    same_vals <- identical(strip(x), strip(y))
    an <- identical(attr_names(x), attr_names(y))
    msg <- if (same_vals) "values identical" else {
      g <- tryCatch(relgap(strip(x), strip(y)), error = function(e) NA)
      sprintf("VALUES DIFFER, max rel gap %s", format(g, digits = 3))
    }
    extra <- ""
    if (!an) {
      na <- unlist(attr_names(x)); nb <- unlist(attr_names(y))
      extra <- sprintf("; attributes added: %s; removed: %s",
                       paste(setdiff(unique(nb), unique(na)), collapse = ","),
                       paste(setdiff(unique(na), unique(nb)), collapse = ","))
    }
    say("%-34s %s%s", nm, msg, extra)
  }
  quit(save = "no")
}

if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
say("ARM %s from %s", arm, find.package("frmtmb"))
res <- list()
keep <- function(nm, expr) {
  r <- tryCatch(suppressWarnings(suppressMessages(expr)),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  if (inherits(r, "frmtmb_conditional_effects") && !is.null(attr(r, "boot"))) {
    attr(r, "boot") <- NULL
  }
  res[[nm]] <<- r
}
hand_draws <- function(fit, n = 40, sd = 0.03, seed = 1) {
  tpl <- fit$frame$par_template
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  set.seed(seed)
  M <- matrix(rep(est, each = n) + rnorm(n * length(est), 0, sd), n,
              dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  structure(list(stanfit = NULL, draws = M, fit = fit),
            class = "frmtmb_draws")
}

set.seed(31)
d <- data.frame(x = rnorm(120), z = runif(120),
                f = factor(sample(c("a", "b", "c"), 120, TRUE)),
                g = factor(rep(1:12, 10)))
d$y <- rnorm(120, 1 + 0.5 * d$x + 0.4 * (d$f == "b") +
               rnorm(12, 0, 0.7)[d$g] + sin(3 * d$z), 0.8)
d$cnt <- rpois(120, exp(0.3 + 0.3 * d$x))
m1 <- frm(bf(y ~ x + f + (1 | g)), family = gaussian(), data = d)
keep("m1 default", conditional_effects(m1))
keep("m1 x:f", conditional_effects(m1, "x:f"))
keep("m1 cond list", conditional_effects(m1, "x", conditions = list(f = "b")))
cd <- data.frame(f = c("a", "c"))
keep("m1 cond df", conditional_effects(m1, "x", conditions = cd))
keep("m1 int_cond", conditional_effects(m1, "x:f",
                                        int_conditions = list(x = c(-1, 1))))
keep("m1 re NULL wald", conditional_effects(m1, "x", re_formula = NULL))
keep("m1 re NULL boot", conditional_effects(m1, "x", re_formula = NULL,
                                            band = "boot", boot = 15,
                                            seed = 1, resolution = 5))
keep("m1 re NA boot", conditional_effects(m1, "x", band = "boot", boot = 15,
                                          seed = 1, resolution = 5))
keep("m1 predict", conditional_effects(m1, "x", method = "predict",
                                       seed = 1, resolution = 5))
keep("m1 profile", conditional_effects(m1, "x", band = "profile",
                                       resolution = 5))
keep("m1 surface=FALSE x:z?", conditional_effects(m1, "x:f",
                                                  surface = FALSE))
m2 <- frm(bf(y ~ s(z) + x), family = gaussian(), data = d)
keep("m2 default", conditional_effects(m2))
keep("m2 boot", conditional_effects(m2, "z", band = "boot", boot = 15,
                                    seed = 1, resolution = 5))
keep("m2 predict nd", predict(m2, newdata = d[1:10, ]))
keep("m2 fitted", fitted(m2))
m3 <- frm(bf(y ~ t2(x, z) + s(z, by = f) + f), family = gaussian(),
          data = d)
keep("m3 default", conditional_effects(m3))
keep("m3 predict nd", predict(m3, newdata = d[1:10, ]))
keep("m3 x:z", conditional_effects(m3, "x:z"))
m4 <- frm(bf(cnt ~ x + (1 | g)), family = poisson(), data = d)
keep("m4 default", conditional_effects(m4))
keep("m4 re NULL", conditional_effects(m4, "x", re_formula = NULL))
d$ord <- cut(d$y, c(-Inf, 0, 1, 2, Inf), ordered_result = TRUE)
m5 <- frm(bf(ord ~ x), family = cumulative(), data = d)
keep("m5 default", conditional_effects(m5))
keep("m5 cats_mean", conditional_effects(m5, categorical = FALSE))
m6 <- frm(bf(y ~ x, sigma ~ x), family = gaussian(), data = d)
keep("m6 sigma", conditional_effects(m6, dpar = "sigma"))
d$y2 <- rnorm(120, d$x)
m7 <- frm(mvbf(bf(y ~ x), bf(y2 ~ x)), family = gaussian(), data = d)
keep("m7 mv", conditional_effects(m7))
m8 <- frm(bf(cnt ~ x), family = zero_inflated_poisson(), data = d)
keep("m8 zip", conditional_effects(m8))
keep("m8 zip boot", conditional_effects(m8, "x", band = "boot", boot = 10,
                                        seed = 1, resolution = 5))
## draws, hand-built
ds1 <- hand_draws(m1)
keep("d1 default", conditional_effects(ds1))
keep("d1 re NULL", conditional_effects(ds1, "x", re_formula = NULL,
                                       seed = 3))
keep("d1 x:f cond", conditional_effects(ds1, "x:f",
                                        conditions = data.frame(g = NA)))
keep("d1 robust F", conditional_effects(ds1, "x", robust = FALSE))
ds2 <- hand_draws(m2)
keep("d2 default", conditional_effects(ds2))
ds5 <- hand_draws(m5)
keep("d5 default", conditional_effects(ds5))
## update and smooth prediction (the frame change)
keep("update m1", coef(update(m1, . ~ . - f)))
keep("update m3 data", coef(update(m3, data = d[1:100, ])))
saveRDS(res, file.path(out_dir, paste0("noreg-", arm, ".rds")))
say("saved %d cases", length(res))
