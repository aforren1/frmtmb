# Reviewer probes (lane ceplot):
#  1. a nonlinear PARAMETER name as an effect (it is on the valid list)
#  2. the ordinal warning under method = "predict" (brms warns only for
#     posterior_epred)
#   Rscript dev/ceplot-rev-misc2.R > dev/ceplot-rev-log/misc2.txt
# Data seeds 11 and 2.
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
show <- function(label, expr) {
  w <- character(0)
  r <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }), error = function(e) paste("ERROR", conditionMessage(e)))
  cat(label, ":", if (is.character(r)) substr(r, 1, 200) else
    paste("answered, effects", paste(names(r), collapse = ",")), "\n")
  for (x in unique(w)) cat("   warning:", substr(gsub("\n", " ", x), 1, 200),
                           "\n")
}
set.seed(11)
dn <- data.frame(x = rnorm(200), z = rnorm(200))
dn$yp <- 2 * exp(0.3 * dn$x) + 0.2 * dn$z + rnorm(200, 0, 0.3)
fn <- frm(bf(yp ~ a * exp(b * x), a ~ 1 + z, b ~ 1, nl = TRUE),
          family = gaussian(), data = dn)
show("1 effects = c('a', 'x')", conditional_effects(fn, effects = c("a", "x"),
                                                    resolution = 3))
show("1 effects = 'b'", conditional_effects(fn, effects = "b",
                                            resolution = 3))
set.seed(2)
do <- data.frame(x = rnorm(200))
do$y <- factor(cut(do$x + rlogis(200), c(-Inf, -1, 0, 1, Inf)),
               ordered = TRUE)
fo <- frm(bf(y ~ x), family = cumulative(), data = do)
show("2 categorical = FALSE, method = 'predict'",
     conditional_effects(fo, categorical = FALSE, method = "predict",
                         resolution = 3))
show("2 categorical = FALSE, method = 'posterior_predict'",
     conditional_effects(fo, categorical = FALSE,
                         method = "posterior_predict", resolution = 3))
