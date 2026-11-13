# Short check of lane sampfront's fix round: value placement of subset rows.
#   Rscript dev/reviews/sampfront-short-place.R
.libPaths(c("C:/Users/adf44/source/r/wt-sampfront-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb.sample", find.package("frmtmb.sample"), "\n")
full_draws <- function(fit, n = 20L, seed = 1L) {
  tpl <- fit$frame[["par_template"]]
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  lab <- frmtmb::brms_par_labels(fit)
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, 0.05),
              n, dimnames = list(NULL, lab))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  structure(list(stanfit = NULL, draws = M, fit = fit),
            class = "frmtmb_draws")
}
q <- function(...) suppressWarnings(suppressMessages(frm(...)))
try_sub <- function(lab, x, nd, meant) {
  cat("\n-- ", lab, "\n", sep = "")
  ll <- suppressWarnings(log_lik(x))
  r <- tryCatch(suppressWarnings(log_lik(x, newdata = nd)),
                error = function(e) e)
  if (inherits(r, "error")) {
    m <- gsub("\\s+", " ", conditionMessage(r))
    cat("  REFUSED:", if (grepl("This newdata", m))
      sub(".*(This newdata)", "\\1", m) else substr(m, 1, 160), "\n")
    if (!is.null(meant)) cat("  (a correct subset was refused)\n")
    return(invisible())
  }
  if (is.null(meant)) {
    cat("  ACCEPTED, dim", dim(r), " (expected a refusal)\n"); return()
  }
  dd <- max(abs(unclass(r) - unclass(ll)[, meant, drop = FALSE]))
  cat("  ACCEPTED, dim", dim(r), "; max |col - fitted col of meant row| =",
      format(dd), if (dd > 0) " SILENT WRONG" else "", "\n")
}
set.seed(11)
d <- data.frame(g = factor(rep(1:6, each = 8)), x = rnorm(48),
                k = rep(1:8, 6))
d$day <- as.Date("2026-01-01") + d$k
d$tm <- as.POSIXct("2026-01-01", tz = "UTC") + 3600 * d$k
d$y <- 1 + 0.5 * d$x + 0.1 * d$k + rnorm(6)[d$g] * 0.7 + rnorm(48)
rn0 <- function(df) { rownames(df) <- NULL; df }
# duplicated fitted rows: rows 41:43 are copies of rows 1:3
dd <- d; dd[41:43, ] <- dd[1:3, ]
xd <- full_draws(q(bf(y ~ x + (1 | g)), data = dd))
try_sub("duplicate fitted rows, value-placed (names reset) rows 41:43",
        xd, rn0(dd[41:43, ]), 41:43)
try_sub("duplicate fitted rows, value-placed rows 1:3", xd, rn0(dd[1:3, ]),
        1:3)
# integer vs double
xk <- full_draws(q(bf(y ~ x + k + (1 | g)), data = d))
dk <- rn0(d[11:15, ]); dk$k <- as.double(dk$k)
try_sub("k fitted integer, newdata double, names reset", xk, dk, 11:15)
dk2 <- d[11:15, ]; dk2$k <- as.double(dk2$k)
try_sub("k fitted integer, newdata double, names kept", xk, dk2, 11:15)
# factor with another level set
dl <- rn0(d[11:15, ]); dl$g <- factor(as.character(dl$g), levels = c(9, 2, 1))
try_sub("g with another level set, names reset", xk, dl, 11:15)
# NA in a read variable
dn <- rn0(d[11:15, ]); dn$x[2] <- NA
try_sub("NA in x, names reset", xk, dn, NULL)
dz <- rn0(d[11:15, ]); dz$x[2] <- -0 * dz$x[2]
# Date and POSIXct read through a transform
xdt <- full_draws(q(bf(y ~ x + as.numeric(day - min(day)) + (1 | g)),
                    data = d))
try_sub("Date variable via as.numeric(day - min(day)), names reset", xdt,
        rn0(d[11:15, ]), 11:15)
xtm <- tryCatch(full_draws(q(bf(y ~ x + as.numeric(tm) + (1 | g)), data = d)),
                error = function(e) e)
if (!inherits(xtm, "error")) {
  try_sub("POSIXct via as.numeric(tm), names reset", xtm, rn0(d[11:15, ]),
          11:15)
  dt2 <- rn0(d[11:15, ]); dt2$tm <- as.POSIXct(format(dt2$tm), tz = "EST")
  try_sub("POSIXct, same clock time in another zone", xtm, dt2, NULL)
} else cat("POSIXct fit failed:", conditionMessage(xtm), "\n")
# a moved value after a name match must not fall through to a wrong row:
# row 12's x set to row 30's values in every variable
dm <- d[12, ]; dm[1, ] <- d[30, ]; rownames(dm) <- "12"
try_sub("row named 12 holding row 30's values", xk, dm, 30)
cat("\nDONE\n")
