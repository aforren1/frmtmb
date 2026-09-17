# Exploration only: where a frame keeps the response, the addition terms,
# the design matrices and the fixed dpars. Run: Rscript dev/brmsport-probe2.R
lib <- Sys.getenv("BRMSPORT_LIB", "C:/Users/adf44/source/r/brmsport-lib")
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
show <- function(label, expr) {
  cat("\n====", label, "\n")
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(r)) print(r) else str(r, max.level = 3, give.attr = FALSE,
                                        list.len = 30)
  invisible(r)
}
set.seed(2)
dat <- data.frame(y = rnorm(9), s = 1:9, w = 1:9, c1 = rep(-1:1, 3),
                  c2 = rep(c("left", "none", "right"), 3),
                  c3 = c(rep(c(TRUE, FALSE), 4), FALSE),
                  c4 = c(sample(-1:1, 5, TRUE), rep(2, 4)), t = 11:19)
fr <- show("cens c4 y+2", frm(y | cens(c4, y + 2) ~ 1, dat,
                               dry_run = "frame")$aterm_values)
show("trials", frm(s | trials(t) ~ 1, dat, family = binomial(),
                    dry_run = "frame")$aterm_values)
show("se", frm(y | se(s) ~ 1, dat, dry_run = "frame")$aterm_values)
f <- frm(y ~ 1 + poly(x, 3), data.frame(y = rnorm(10), x = rnorm(10)),
         dry_run = "frame")
show("linpred mu names", names(f$linpreds[[1]]))
show("linpred mu", f$linpreds[[1]][c("X", "dpar", "resp")])
show("ordinal", frm(y ~ x, data.frame(y = c(1:5, 1:4, 4), x = rnorm(10)),
                    family = cumulative(), dry_run = "frame")[
                      c("y", "y_levels", "par_template")])
show("fixed nu", frm(bf(y ~ 1, nu = 3), list(y = 1:10), family = student(),
                     dry_run = "frame")[c("par_template", "map",
                                          "betad_fixed_idx")])
show("mv", frm(bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE),
               data.frame(y1 = 1:10, y2 = 11:20, x = rep(0, 10)),
               dry_run = "frame")[c("y", "n_obs", "spec")])
show("rate", frm(y | rate(time) ~ x, data.frame(y = rpois(10, 1),
                 x = rnorm(10), time = 1:10), family = poisson(),
                 dry_run = "frame")$aterm_values)
show("re", frm(y ~ x + (1 | g), data.frame(y = rnorm(10), x = rnorm(10),
               g = rep(1:5, 2)), dry_run = "frame")$re_blocks)
show("dry_run values", formals(frm)$dry_run)
