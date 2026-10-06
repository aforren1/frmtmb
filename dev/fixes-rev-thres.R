# Reviewer of lane fixes, claim 4a: core ord_thres_linpred() at frmtmb's
# estimates against brms 2.23.0's posterior_linpred(incl_thres = TRUE)
# at the same values (Fixed_param draws, brms_fixed_fit()), in sample
# and on newdata, for every shape the lane names and three it did not.
#   Rscript dev/fixes-rev-thres.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = normalizePath("dev/fixes-rev-stan-cache"))
suppressMessages({library(testthat); library(frmtmb)})
cat("LIB", find.package("frmtmb"), "\n")
henv <- new.env(parent = asNamespace("frmtmb"))
for (h in c("helper-brms.R", "helper-brms-methods.R")) {
  sys.source(file.path("tests/testthat", h), envir = henv)
}
brms_fixed_fit <- henv$brms_fixed_fit
set.seed(8080)
n <- 220
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
u <- stats::rlogis(n) / exp(0.3 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
d$yg <- ifelse(d$g == "b", pmin(d$y, 4L), d$y)
d$yh <- ifelse(runif(n) < 0.25, 0L, d$y)
nd <- data.frame(x = c(-1, 0, 2), z = c(0.5, -1, 0),
                 g = factor(c("b", "a", "b"), levels = c("a", "b")),
                 y = 1L, yg = 1L)
shapes <- list(
  cum = list(bf(y ~ x), cumulative(), brms::bf(y ~ x), brms::cumulative()),
  cum_disc = list(bf(y ~ x, disc ~ 0 + z), cumulative(),
                  brms::bf(y ~ x, disc ~ 0 + z), brms::cumulative()),
  cum_eq_disc = list(bf(y ~ x, disc ~ 0 + z),
                     cumulative(threshold = "equidistant"),
                     brms::bf(y ~ x, disc ~ 0 + z),
                     brms::cumulative(threshold = "equidistant")),
  sr_cs = list(bf(y ~ x + cs(z)), sratio(), brms::bf(y ~ x + cs(z)),
               brms::sratio()),
  acat_cs = list(bf(y ~ x + cs(z)), acat(), brms::bf(y ~ x + cs(z)),
                 brms::acat()),
  cr_cs_disc = list(bf(y ~ x + cs(z), disc ~ 0 + x), cratio(),
                    brms::bf(y ~ x + cs(z), disc ~ 0 + x), brms::cratio()),
  acat_eq = list(bf(y ~ x), acat(threshold = "equidistant"),
                 brms::bf(y ~ x), brms::acat(threshold = "equidistant")),
  sr_stz = list(bf(y ~ x), sratio(threshold = "sum_to_zero"),
                brms::bf(y ~ x), brms::sratio(threshold = "sum_to_zero")),
  cum_gr = list(bf(yg | thres(gr = g) ~ x), cumulative(),
                brms::bf(yg | thres(gr = g) ~ x), brms::cumulative()),
  sr_gr_cs = list(bf(yg | thres(gr = g) ~ x + cs(z)), sratio(),
                  brms::bf(yg | thres(gr = g) ~ x + cs(z)), brms::sratio()),
  cum_gr_disc = list(bf(yg | thres(gr = g) ~ x, disc ~ 0 + z), cumulative(),
                     brms::bf(yg | thres(gr = g) ~ x, disc ~ 0 + z),
                     brms::cumulative()),
  hurdle = list(bf(yh ~ x), hurdle_cumulative(), brms::bf(yh ~ x),
                brms::hurdle_cumulative()))
cmp <- function(a, b) {
  if (!identical(dim(a), dim(b))) return(paste("DIM", paste(dim(a),
                                                            collapse = "x"),
                                               "vs", paste(dim(b),
                                                           collapse = "x")))
  na_same <- identical(is.na(a), is.na(b))
  m <- max(abs(a - b) / pmax(abs(b), 1), na.rm = TRUE)
  sprintf("NA pattern same %s, max rel diff %.3g", na_same, m)
}
for (nm in names(shapes)) {
  s <- shapes[[nm]]
  cat("==", nm, "\n")
  fit <- tryCatch(suppressWarnings(frm(s[[1]], family = s[[2]], data = d)),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    cat("  frm ERROR:", conditionMessage(fit), "\n"); next
  }
  bb <- tryCatch(brms_fixed_fit(s[[3]], s[[4]], d, fit, ndraws = 3),
                 error = function(e) e)
  if (inherits(bb, "error")) {
    cat("  brms fixed fit ERROR:", conditionMessage(bb), "\n"); next
  }
  pl <- tryCatch(brms::posterior_linpred(bb, incl_thres = TRUE),
                 error = function(e) e)
  if (inherits(pl, "error")) {
    cat("  brms ERROR:", conditionMessage(pl), "\n"); next
  }
  cat("  brms dim:", dim(pl), " dimnames[[3]]:", dimnames(pl)[[3]],
      " dimnames[[1]]:", dimnames(pl)[[1]], "\n")
  if (nm == "hurdle") {
    cat("  brms draw 1 rows 1-2:\n"); print(pl[1, 1:2, ])
    hu <- brms::posterior_epred(bb, dpar = "hu")[1, 1:2]
    cat("  brms hu rows 1-2:", hu, "\n")
    r <- tryCatch(frmtmb:::ord_thres_linpred(fit),
                  error = function(e) conditionMessage(e))
    cat("  frm:", if (is.character(r)) r else "returned", "\n")
    next
  }
  a <- frmtmb:::ord_thres_linpred(fit)
  cat("  in sample:", cmp(a, pl[1, , ]), "\n")
  pn <- tryCatch(brms::posterior_linpred(bb, incl_thres = TRUE,
                                         newdata = nd)[1, , , drop = TRUE],
                 error = function(e) conditionMessage(e))
  an <- tryCatch(frmtmb:::ord_thres_linpred(fit, newdata = nd),
                 error = function(e) conditionMessage(e))
  if (is.character(pn) || is.character(an)) {
    cat("  newdata: brms", if (is.character(pn)) pn else "ok", "| frm",
        if (is.character(an)) an else "ok", "\n")
  } else {
    cat("  newdata:", cmp(an, pn), "\n")
    print(rbind(frm = an[1, ], brms = pn[1, ]))
  }
}
