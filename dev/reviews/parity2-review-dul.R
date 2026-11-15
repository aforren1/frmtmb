# Review of lane parity2, item 5337: drop_unused_levels = FALSE on an
# interaction or nested grouping term, brms 2.23.0 standata (no compile)
# against frmtmb (arm base|lane). Seed 3.
#   Rscript dev/reviews/parity2-review-dul.R base|lane
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "lane"
libs <- c("C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-parity2-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, " lib:", find.package("frmtmb"),
    as.character(packageVersion("frmtmb")), " brms",
    as.character(packageVersion("brms")), "\n")
set.seed(3)
n <- 120
du <- data.frame(g = factor(sample(1:6, n, TRUE), levels = 1:7),
                 h = factor(sample(1:4, n, TRUE)), x = rnorm(n))
du$y <- rnorm(n, 1 + 0.4 * du$x)
cat("observed g:h combos:", nrow(unique(du[c("g", "h")])), "\n")
for (f in c("y ~ x + (1 | g)", "y ~ x + (1 | g:h)", "y ~ x + (1 | g/h)",
            "y ~ x + (1 | h) + (1 | h:g)")) {
  fo <- stats::as.formula(f)
  sd <- brms::make_standata(fo, data = du, drop_unused_levels = FALSE)
  nn <- unlist(sd[grep("^N_[0-9]+$", names(sd))])
  sd0 <- brms::make_standata(fo, data = du)
  nn0 <- unlist(sd0[grep("^N_[0-9]+$", names(sd0))])
  fit <- tryCatch(suppressMessages(suppressWarnings(
    frm(fo, family = gaussian(), data = du, drop_unused_levels = FALSE))),
    error = function(e) e)
  fr <- if (inherits(fit, "error")) conditionMessage(fit) else
    paste(vapply(ranef(fit), function(r) NROW(r), 1L), collapse = " ")
  fit0 <- suppressMessages(suppressWarnings(
    frm(fo, family = gaussian(), data = du)))
  fr0 <- paste(vapply(ranef(fit0), function(r) NROW(r), 1L), collapse = " ")
  cat(sprintf("%-30s brms N (FALSE) %s (default) %s | frmtmb levels (FALSE) %s (default) %s\n",
              f, paste(nn, collapse = " "), paste(nn0, collapse = " "), fr,
              fr0))
  if (!inherits(fit, "error")) {
    cat("   logLik FALSE", as.numeric(logLik(fit)), " default",
        as.numeric(logLik(fit0)), "\n")
  }
}
