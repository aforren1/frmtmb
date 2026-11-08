# Short check of the 0.80.0 release fixes, item 2: `0 + f + w1 * w2` at 1e5
# and 1e6 (columns span the constant), and two designs whose columns do
# not span it, `0 + w1 + w2` and `0 + w1 * w2`, against the base build's
# frame. The rank decision is read from the frame (dry_run = "frame").
#
#   Rscript dev/review-rel-noint2.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(strsplit(lib, ";")[[1]], "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "\n")
run <- function(expr) {
  cs <- character()
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(w) {
      cs <<- c(cs, paste0("W:", class(w)[1L]))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      cs <<- c(cs, paste0("M:", substr(conditionMessage(m), 1, 50)))
      invokeRestart("muffleMessage")
    }), error = function(e) {
      cs <<- c(cs, paste0("E:", substr(conditionMessage(e), 1, 80)))
      NULL
    })
  list(v = v, cs = unique(cs))
}
ll <- function(f) if (is.null(f)) NA else as.numeric(logLik(f))
for (seed in 1:3) for (m in c(1e3, 1e5, 1e6)) {
  set.seed(seed)
  n <- 300
  d <- data.frame(z1 = rnorm(n), z2 = rnorm(n),
                  f = factor(sample(letters[1:3], n, TRUE)))
  d$y <- c(a = 0, b = 1, c = -1)[as.character(d$f)] + 0.5 * d$z1 -
    0.4 * d$z2 + 0.6 * d$z1 * d$z2 + rnorm(n, 0, 0.5)
  d$w1 <- d$z1 + m
  d$w2 <- d$z2 + m
  # the no-span designs: the data are generated with no intercept so the
  # model on w is a fair one, and its "centered" twin is the same columns
  d$y0 <- 0.5 * d$w1 - 0.5 * d$w2 + 0.6 * d$z1 * d$z2 + rnorm(n, 0, 0.5)
  cases <- list(
    span_prod = list(y ~ 0 + f + w1 * w2, y ~ 0 + f + z1 * z2),
    nospan_add = list(y0 ~ 0 + w1 + w2, NULL),
    nospan_prod = list(y0 ~ 0 + w1 * w2, NULL))
  for (nm in names(cases)) {
    a <- run(frm(cases[[nm]][[1]], data = d))
    b <- if (!is.null(cases[[nm]][[2]])) run(frm(cases[[nm]][[2]], data = d))
    fr <- run(frm(cases[[nm]][[1]], data = d, dry_run = "frame"))
    kept <- if (!is.null(a$v)) length(a$v$estimates$beta) else NA
    cat(sprintf("seed %d m %.0e %-11s ll %.4f beta %s%s | %s\n", seed, m, nm,
                ll(a$v), kept,
                if (!is.null(b)) sprintf(" centered %.4f short %.4f",
                                         ll(b$v), ll(b$v) - ll(a$v)) else "",
                if (length(a$cs)) paste(a$cs, collapse = " ; ") else "-"))
  }
}
cat("NOINT2 DONE\n")
