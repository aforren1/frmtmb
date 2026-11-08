# The whole warning of the reviewer's gp + new-level call (review m5). Seed 8150.
# the row check's verdict on it. Data seed 8150.
source("C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/mfxse-lib.R")
ns <- asNamespace("frmtmb")
set.seed(8150)
mg <- 90
dg <- data.frame(x = runif(mg, 0, 6), g = factor(sample(1:6, mg, TRUE)))
dg$y <- sin(dg$x) + rnorm(6, 0, 0.5)[dg$g] + rnorm(mg, sd = 0.3)
fg <- quiet(frm(y ~ gp(x) + (1 | g), data = dg))
ndg <- data.frame(x = c(dg$x[1], 3.33, 9), g = factor(c("1", "new", "2")))
seen <- NULL
f0 <- get("get_predict.frmtmb_fit", envir = ns)
registerS3method("get_predict", "frmtmb_fit", function(model, newdata, ...) {
  if (is.null(seen)) seen <<- newdata
  f0(model, newdata, ...)
}, envir = asNamespace("marginaleffects"))
p <- quiet(predictions(fg, newdata = ndg, allow_new_levels = TRUE))
print(seen)
lp <- fg$frame$linpreds[[1]]
ed <- ns$lp_eta_design(fg, lp, seen, TRUE, TRUE)
cat("n", ed$n, "nrow", nrow(seen), "\n")
for (sp in ed$sm_parts) print(sp$extra_var)
# the whole message of the call, untruncated
p <- withCallingHandlers(predictions(fg, newdata = ndg, allow_new_levels = TRUE),
                         warning = function(w) {
                           cat("W:", conditionMessage(w), "\n")
                           invokeRestart("muffleWarning")
                         })
