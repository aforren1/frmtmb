# Release review of 0.82.0, probe 2: does the restart's message survive
# an escape that runs after it (escape_stationary(), mo_search(),
# sd_up_escape())? On the lane (nlstop) it was said at the restart.
# usage: probe2.R <lib>
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r21"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
rec <- new.env()
trace_esc <- function() {
  rec$restart <- FALSE
  rec$stat_gain <- NA
}
# record whether the restart won, without changing what it returns
orig <- get("prior_conflict_escape", ns)
unlockBinding("prior_conflict_escape", ns)
assign("prior_conflict_escape", function(...) {
  o <- orig(...)
  rec$restart <- !is.null(o$prior_escape)
  o
}, envir = ns)
run <- function(label, expr) {
  trace_esc()
  msgs <- character()
  fit <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) {
      msgs <<- c(msgs, class(m)[1])
      invokeRestart("muffleMessage")
    })
  if (inherits(fit, "error")) {
    cat(label, "ERROR", conditionMessage(fit), "\n")
    return(invisible())
  }
  cat(sprintf("%-28s restart_won=%s record=%s message=%s stat_gain=%s mo=%s sd_esc=%s ll=%.4f\n",
              label, rec$restart, !is.null(fit$opt$prior_escape),
              "frmtmb_prior_multistart" %in% msgs,
              format(fit$opt$stationary_escape[["gain"]] %||% NA, digits = 3),
              format(fit$opt$mo_search[["gain"]] %||% NA, digits = 3),
              format(fit$opt$sd_escape[["gain"]] %||% NA, digits = 3),
              -fit$opt$objective))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
pr <- c(set_prior("normal(0, 1000)", class = "b"),
        set_prior("normal(1e6, 1e6)", class = "Intercept"))
for (seed in 6201:6206) {
  set.seed(seed)
  d <- data.frame(x = rnorm(200))
  # skewed residuals: the skew_normal fit has a stationary point at alpha 0
  e <- rexp(200, 1 / 1e3) - 1e3
  d$y <- 1e6 * (1 + 0.05 * d$x) + e
  run(paste0("skew_normal seed ", seed),
      frm(y ~ x, data = d, family = skew_normal(), prior = pr))
  run(paste0("gaussian seed ", seed),
      frm(y ~ x, data = d, prior = pr))
  # a group term whose sd the stop leaves near 0 (sd_up_escape())
  d$g <- factor(rep(1:20, each = 10))
  d$y2 <- d$y + rnorm(20, 0, 30)[d$g]
  run(paste0("gaussian (1|g) seed ", seed),
      frm(y2 ~ x + (1 | g), data = d, prior = pr))
  # a monotonic term (mo_search())
  d$m <- sample(1:4, 200, TRUE)
  d$y3 <- d$y + 2e3 * c(0, 0.1, 0.9, 1)[d$m]
  run(paste0("gaussian mo(m) seed ", seed),
      frm(y3 ~ x + mo(m), data = d, prior = pr))
}
