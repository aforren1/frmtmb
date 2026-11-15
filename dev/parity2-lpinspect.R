# Lane parity2: what a fit's linear predictors hold, per shape, to write
# brms's default-effects rule against (base build).
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-parity2/dev/parity2-data.R")
show <- function(f, fam = gaussian()) {
  fit <- suppressWarnings(suppressMessages(frm(f, family = fam, data = d)))
  for (nm in names(fit$frame$linpreds)) {
    lp <- fit$frame$linpreds[[nm]]
    cat("--", nm, " resp", lp$resp, " dpar", lp$dpar, "\n")
    if (!is.null(lp$terms)) {
      cat("  term.labels:", attr(lp$terms, "term.labels"), "\n")
    }
    cat("  slots:", paste(names(lp), collapse = " "), "\n")
    for (s in c("smooths", "mo", "mi", "me", "cs", "gps")) {
      for (e in lp[[s]] %||% list()) {
        cat("  ", s, ": label", e$label %||% "", " term",
            deparse1(e$sm$term %||% e$expr %||% e$var), " mult",
            deparse1(e$mult_expr), " by", deparse1(e$sm$by %||% ""), "\n")
      }
    }
    if (!is.null(lp$nl_body)) cat("  nl_body:", deparse1(lp$nl_body),
                                  " data_list:", names(lp$data_list),
                                  " nl_pars:", lp$nl_pars, "\n")
  }
  cat("rspec dpars:", names(fit$spec$responses[[1]]$dpars), " nlpars:",
      fit$spec$responses[[1]]$nlpars, "\n")
}
show(bf(y ~ exp(x):z))
show(bf(y ~ t2(x, z)))
show(bf(y ~ mo(k):z + s(w, by = x)))
show(bf(yo ~ cs(x) + z), sratio())
show(bf(y ~ x, sigma ~ exp(z)))
show(bf(y ~ a * exp(b * x), a ~ 1 + w, b ~ z, nl = TRUE))
show(bf(y ~ x + (1 + z | g) + offset(w) + me(x2, w) + gp(z, w)))
show(bf(y ~ 0 + Intercept + x))
