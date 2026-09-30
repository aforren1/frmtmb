# Reviewer, claim 1 supplement: the paths the first battery did not reach
# (the expected-category curve, categorical = FALSE; ord_cat_moments();
# the ordinal response/pearson residual values that plot() uses;
# pp_check() with a newdata response; influence() refits through
# thres_pin_recode()). Same data as dev/fams2-rev-ord.R.
#   Rscript dev/fams2-rev-ord2.R base|lane
arm <- commandArgs(TRUE)[1]
base_lib <- c("C:/Users/adf44/source/r/rellib-r3",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lane") c("C:/Users/adf44/source/r/wt-fams2-lib",
                               base_lib) else base_lib)
suppressPackageStartupMessages(library(frmtmb))
out_dir <- "C:/Users/adf44/source/r/frmtmb-wt-fams2/dev/fams2-rev-out"
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(cw) {
      w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning")
    })
  if (inherits(v, "ggplot")) v <- v$data
  if (inherits(v, "frmtmb_conditional_effects")) v <- unclass(v)
  list(value = v, warnings = w)
}
set.seed(20260929)
n <- 400; ng <- 20
g <- factor(rep(seq_len(ng), length.out = n))
x <- rnorm(n)
fc <- factor(sample(c("a", "b", "c"), n, TRUE))
grp <- factor(sample(c("u", "v"), n, TRUE))
ue <- rnorm(ng, 0, 0.7)[g]
lat <- 0.9 * x + ue + 0.5 * (fc == "b") - 0.4 * (fc == "c") + rlogis(n)
y <- as.integer(cut(lat, c(-Inf, -1, 0.2, 1.3, Inf)))
lat2 <- -0.6 * x + rlogis(n)
y2 <- as.integer(cut(lat2, c(-Inf, -0.5, 0.8, Inf)))
d <- data.frame(y = y, y2 = y2, x = x, g = g, fc = fc, grp = grp)
d$yf <- factor(c("lo", "mid", "hi", "top")[d$y],
               levels = c("lo", "mid", "hi", "top"), ordered = TRUE)
d$yu <- factor(c("lo", "mid", "hi", "top")[d$y],
               levels = c("lo", "never", "mid", "hi", "top", "nobody"),
               ordered = TRUE)
pr <- set_prior("normal(0, 3)", class = "Intercept")
fits <- list(
  cum_logit_re = frm(y ~ x + fc + (1 | g), data = d, family = cumulative()),
  cum_probit = frm(y ~ x, data = d, family = cumulative(link = "probit")),
  thresK = suppressWarnings(frm(y | thres(5) ~ x, data = d,
                                family = cratio(), prior = pr)),
  thres_gr = frm(y | thres(gr = grp) ~ x, data = d, family = cumulative()),
  sratio_cs = frm(y ~ x + cs(fc), data = d, family = sratio()),
  acat_cs = frm(y ~ x + cs(fc), data = d, family = acat()),
  cratio_ofactor = frm(yf ~ x + (1 | g), data = d, family = cratio()),
  cum_unused = frm(yu ~ x, data = d, family = cumulative())
)
res <- list()
for (nm in names(fits)) {
  f <- fits[[nm]]
  rs <- f$spec$responses[[1]]
  nd <- d[1:6, ]
  r <- list(
    ce_expected = cap(conditional_effects(f, categorical = FALSE)),
    ce_expected_re = cap(conditional_effects(f, categorical = FALSE,
                                             re_formula = NULL)),
    moments = cap(frmtmb:::ord_cat_moments(f, rs)),
    rv_pearson = cap(frmtmb:::residual_values(f, type = "pearson")),
    rv_response = cap(frmtmb:::residual_values(f, type = "response")),
    ppc_nd = cap({ set.seed(3); pp_check(f, type = "bars", ndraws = 4,
                                         newdata = nd) }),
    ppc_ecdf = cap({ set.seed(4); pp_check(f, type = "ecdf_overlay",
                                           ndraws = 4) }),
    predict_nd = cap({ set.seed(5); predict(f, newdata = nd) })
  )
  res[[nm]] <- r
}
# influence() refits on subsets, through thres_pin_recode()
d50 <- d[1:60, ]
d50$yu2 <- droplevels(d50$yu)
d50$yu2 <- factor(d50$yu, levels = levels(d$yu), ordered = TRUE)
fi <- frm(yu2 ~ x, data = d50, family = cumulative())
res$influence_unused <- cap(unclass(influence(fi)))
fi2 <- frm(yf ~ x, data = d50, family = sratio())
res$influence_ofactor <- cap(unclass(influence(fi2)))
mv <- frm(bf(y ~ x) + cumulative() + bf(y2 ~ x) + acat(), data = d)
res$mv_ce_expected <- cap(conditional_effects(mv, categorical = FALSE))
res$mv_ce_expected_y2 <- cap(conditional_effects(mv, resp = "y2",
                                                 categorical = FALSE))
# strip functions and environments
strip <- function(x) {
  if (is.function(x) || is.environment(x)) return("<fn/env>")
  if (is.list(x)) { a <- attributes(x); y <- lapply(x, strip)
    a$names <- NULL; attributes(y) <- c(list(names = names(x)), a); y } else x
}
saveRDS(strip(res), file.path(out_dir, paste0("ord2-", arm, ".rds")))
cat("saved\n")
