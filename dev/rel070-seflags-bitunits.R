# rel070 seflags: optrap's response unit and Newton polish on the 30 fits
# of dev/seflags-rev4-bitwise.R (same data and specs), rellib-r8.
.libPaths(c("C:/Users/adf44/source/r/rellib-r8", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb)); ns <- asNamespace("frmtmb")
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-release/dev/seflags-rev4-bitwise.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^res <- list", src) - 1)]))
for (s in 1:3) { d <- mk(s); for (nm in names(specs)) { sp <- specs[[nm]]
  f <- suppressWarnings(suppressMessages(frm(sp[[1]], family = sp[[2]], data = d)))
  u <- ns$predictor_units(f$frame)
  raw <- tryCatch({lp <- f$frame$linpreds[[1]]; y <- f$frame$y[[lp$resp]]; ns$predictor_resid_sd(lp, y, is.finite(y))}, error = function(e) NA)
  cat(sprintf("%-7s seed %d | units %s | raw spread %.4g | polish %s\n", nm, s,
    if (length(u)) paste(names(u), u) else "- (1)", raw,
    if (is.null(f$opt$polish)) "-" else paste(signif(f$opt$polish, 3), collapse = " ")))
}}
