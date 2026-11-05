# Short check of fix round 2, item 5: where a categorical me() contrast with a row of NA
# stops. Seed 7421, the recipe of review-predread-r2-na.R.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7421)
n <- 240
d <- data.frame(z = rnorm(n), lat = rnorm(n), sx = 0.4)
d$xo <- d$lat + rnorm(n, 0, 0.4)
eta <- 0.4 + 0.7 * d$lat
d$cat <- as.integer(cut(eta + rlogis(n), c(-Inf, 0, 1.5, Inf), labels = FALSE))
nd <- data.frame(xo = c(1, NA, 0, 1), sx = c(0.4, 0.4, NA, 0.4), z = c(0, 0, 0, 1))
fit <- quiet(frm(bf(cat ~ me(xo, sx) + z), family = categorical(), data = d))
r <- withCallingHandlers(tryCatch(suppressWarnings(frm_curve(fit, newdata = nd, contrast = transform(nd, z = z + 1), simultaneous = FALSE)$.se), error = function(e) conditionMessage(e)),
  error = function(e) {
    cs <- sys.calls()
    cat("CALLS:\n ", paste(tail(vapply(cs, function(cl) substr(paste(deparse(cl), collapse = " "), 1, 90), ""), 9), collapse = "\n  "), "\n")
  })
print(r)
print(body(frmtmb:::me_same_value)[[1]])
g <- frmtmb:::curve_sim_crit
print(names(formals(g)))
