# Probe for the guardrail-1 test: brms 2.23.0's formula vocabulary
# (families, addition terms, predictor specials) and how frmtmb on the
# library given in SPATIAL_LIBS answers each name with no extension
# attached. Prints one line per name.
libs <- strsplit(Sys.getenv("SPATIAL_LIBS"), ";", fixed = TRUE)[[1]]
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
bns <- asNamespace("brms")
fams <- sub("^[.]family_", "", ls(bns, all.names = TRUE,
                                   pattern = "^[.]family_"))
ats <- sub("^resp_", "", grep("^resp_", getNamespaceExports("brms"),
                              value = TRUE))
cat("BRMS families", length(fams), ":", paste(sort(fams), collapse = " "),
    "\n")
cat("BRMS aterms", length(ats), ":", paste(sort(ats), collapse = " "), "\n")
sp <- c("s", "t2", "gp", "mo", "me", "mi", "cs", "cse", "car", "sar", "fcor",
        "ar", "ma", "arma", "cosy", "unstr", "gr", "mm", "mmc")
cat("BRMS specials exported:",
    paste(sp[sp %in% getNamespaceExports("brms")], collapse = " "),
    "| not exported:",
    paste(sp[!sp %in% getNamespaceExports("brms")], collapse = " "), "\n")
for (f in sort(fams)) {
  r <- tryCatch({ brmsfamily(f); "OK" },
                error = function(e) paste("ERR", conditionMessage(e)))
  cat("FAM", f, substr(r, 1, 150), "\n")
}
d <- data.frame(y = rpois(20, 3), x = rnorm(20), n = 10, g = gl(4, 5),
                w = 1, ce = 0, t = 1:20)
for (a in sort(ats)) {
  fo <- stats::as.formula(paste0("y | ", a, "(w) ~ x"))
  r <- tryCatch({ frm(bf(fo), family = gaussian(), data = d,
                      dry_run = "spec"); "OK" },
                error = function(e) paste("ERR", conditionMessage(e)))
  cat("ATERM", a, substr(gsub("\n", " ", r), 1, 150), "\n")
}
