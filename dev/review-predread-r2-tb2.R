source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7421)
n <- 240
d <- data.frame(z = rnorm(n), lat = rnorm(n), sx = 0.4)
d$xo <- d$lat + rnorm(n, 0, 0.4)
eta <- 0.4 + 0.7 * d$lat
d$cat <- as.integer(cut(eta + rlogis(n), c(-Inf, 0, 1.5, Inf), labels = FALSE))
nd <- data.frame(xo = c(1, NA, 0, 1), sx = c(0.4, 0.4, NA, 0.4), z = c(0, 0, 0, 1))
fit <- quiet(frm(bf(cat ~ me(xo, sx) + z), family = categorical(), data = d))
tryCatch(withCallingHandlers(suppressWarnings(frm_curve(fit, newdata = nd, contrast = transform(nd, z = z + 1), simultaneous = FALSE)$.se),
  error = function(e) {
    cs <- sys.calls()
    cat("CALLS:\n ", paste(tail(vapply(cs, function(cl) substr(paste(deparse(cl), collapse = " "), 1, 100), ""), 8), collapse = "\n  "), "\n")
  }), error = function(e) print(conditionMessage(e)))
