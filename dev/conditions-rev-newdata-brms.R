# Reviewer, lane wt-conditions (recheck 1): the newdata refusals case by
# case against brms 2.23.0's own validate_newdata() on an empty brmsfit
# (no Stan compile), and against frmtmb's predict() on base and lane.
#   Rscript dev/conditions-rev-newdata-brms.R base|lane   (seed 20260917)
arm <- commandArgs(TRUE)[1L]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/conditions-lib"
            else "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
set.seed(20260917)
n <- 60
d <- data.frame(y = rnorm(n), x = rnorm(n),
                g = factor(rep(1:10, each = 6)),
                f = factor(rep(c("a", "b", "c"), 20)),
                w = runif(n, 0.5, 1.5), cc = rbinom(n, 1, 0.2),
                k = rbinom(n, 10, 0.4), ntr = 10L)
nd <- d[1:6, ]
v <- list(
  base = nd,
  new_g = transform(nd, g = factor(c("99", "1", "2", "3", "4", "98"))),
  no_g = nd[, setdiff(names(nd), "g")],
  no_y = nd[, setdiff(names(nd), c("y", "k"))],
  char_f = transform(nd, f = as.character(f)),
  na_f = transform(nd, f = factor(c("a", NA, "b", "c", "a", "b"))),
  new_f = transform(nd, f = factor(c("a", "zz", "b", "c", "a", "b"))),
  no_x = nd[, setdiff(names(nd), "x")],
  no_w = nd[, setdiff(names(nd), "w")],
  no_cc = nd[, setdiff(names(nd), "cc")],
  no_ntr = nd[, setdiff(names(nd), "ntr")],
  sub_f = droplevels(nd[nd$f == "a", ]))
specs <- list(
  mixed = list(fo = y ~ x + f + (1 | g), fam = gaussian()),
  slope = list(fo = y ~ x + (f | g), fam = gaussian()),
  wt = list(fo = y | weights(w) ~ x + f, fam = gaussian()),
  cens = list(fo = y | cens(cc) ~ x + f, fam = gaussian()),
  binom = list(fo = k | trials(ntr) ~ x + f + (1 | g), fam = binomial()))
settings <- list(keep = list(re = NULL, anl = FALSE),
                 allow = list(re = NULL, anl = TRUE),
                 reNA = list(re = NA, anl = FALSE))
one <- function(expr) tryCatch({force(expr); "OK"}, error = function(e)
  paste0("ERR: ", substr(gsub("[[:space:]]+", " ", conditionMessage(e)), 1, 90)))
rows <- list()
for (s in names(specs)) {
  fit <- suppressWarnings(frm(specs[[s]]$fo, data = d, family = specs[[s]]$fam))
  bf_ <- suppressWarnings(suppressMessages(brms::brm(specs[[s]]$fo, data = d,
    family = specs[[s]]$fam, empty = TRUE)))
  for (vn in names(v)) for (st in names(settings)) {
    re <- settings[[st]]$re; anl <- settings[[st]]$anl
    fr <- one(predict(fit, newdata = v[[vn]], re_formula = re,
                      allow_new_levels = anl))
    br <- one(brms:::validate_newdata(v[[vn]], bf_, re_formula = re,
                                      allow_new_levels = anl))
    rows[[length(rows) + 1L]] <- data.frame(spec = s, newdata = vn,
      setting = st, frmtmb = fr, brms = br, stringsAsFactors = FALSE)
  }
}
out <- do.call(rbind, rows)
saveRDS(out, sprintf("C:/Users/adf44/source/r/frmtmb-wt-conditions/dev/conditions-rev-log/newdata-brms-%s.rds", arm))
cat("rows", nrow(out), "\n")
