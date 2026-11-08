# Review of lane mfxse: why avg_slopes() at new levels gives no warning
# (follow.log): replay interop_rows_check()'s design rebuild on the newdata
# the armed get_predict() receives. Data seed 7411 (mre_data()).
lane <- "C:/Users/adf44/source/r/wt-mfxse-lib"
.libPaths(c(lane, "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(marginaleffects)})
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
quiet <- function(e) suppressWarnings(suppressMessages(e))
set.seed(7411)
J <- 12L; m <- 20L; n <- J * m
d <- data.frame(x = rnorm(n), f = factor(sample(c("a", "b"), n, TRUE)),
                g = factor(rep(seq_len(J), each = m)))
d$cnt <- rpois(n, exp(0.5 + 0.4 * d$x + rnorm(J, 0, 0.8)[d$g]))
fit <- quiet(frm(cnt ~ x + f + (1 | g), family = poisson(), data = d))
ndn <- data.frame(x = c(0, 1, 2), f = "a", g = factor(c("3", "new", "new2")))
seen <- NULL
gp0 <- ns$get_predict.frmtmb_fit
registerS3method("get_predict", "frmtmb_fit", function(model, newdata, type = "response", ...) {
  if (is.null(seen) && isTRUE(model$cache[["mfx_rows_armed"]])) seen <<- list(nd = newdata, dots = list(...))
  gp0(model, newdata, type, ...)
}, envir = asNamespace("marginaleffects"))
invisible(quiet(avg_slopes(fit, newdata = ndn, variables = "x", allow_new_levels = TRUE)))
cat("class of newdata:", class(seen$nd), "\n")
cat("columns:", names(seen$nd), "\n")
cat("dots:", names(seen$dots), "\n")
lp <- fit$frame$linpreds[[1]]
r <- tryCatch(ns$lp_eta_design(fit, lp, seen$nd, TRUE, isTRUE(seen$dots$allow_new_levels)),
              error = function(e) paste("ERROR:", conditionMessage(e)))
if (is.character(r)) cat(r, "\n") else cat("design rows", r$n, "\n")
r2 <- tryCatch(ns$lp_eta_design(fit, lp, as.data.frame(seen$nd), TRUE, TRUE),
               error = function(e) paste("ERROR:", conditionMessage(e)))
if (is.character(r2)) cat("as.data.frame:", r2, "\n") else cat("as.data.frame: design rows", r2$n, "\n")
# the check itself on that newdata, armed by hand
fit$cache[["mfx_rows_armed"]] <- TRUE
w <- tryCatch(ns$interop_rows_check(fit, seen$nd, seen$dots, seen$nd$rowid > 0),
              warning = function(w) paste("WARNS:", class(w)[1], substr(conditionMessage(w), 1, 70)))
cat("interop_rows_check on that newdata:", if (is.character(w)) w else "no warning", "\n")
# does marginaleffects muffle warnings raised inside the slope Jacobian?
ms <- grep("suppressWarnings", deparse(get("comparisons", asNamespace("marginaleffects"))), value = TRUE)
cat("suppressWarnings lines in comparisons():", length(ms), "\n"); cat(trimws(ms), sep = "\n")
# every warning of avg_slopes() at new levels, with its class
ws <- list()
invisible(withCallingHandlers(avg_slopes(fit, newdata = ndn, variables = "x", allow_new_levels = TRUE),
  warning = function(w) { ws[[length(ws) + 1]] <<- w; invokeRestart("muffleWarning") }))
for (w in ws) cat("avg_slopes warning:", paste(class(w), collapse = "/"), "|", substr(conditionMessage(w), 1, 80), "\n")
