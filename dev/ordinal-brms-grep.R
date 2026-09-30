.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
ns <- asNamespace("brms")
nm <- ls(ns, all.names = TRUE)
body_of <- function(f) {
  x <- get(f, ns)
  if (!is.function(x)) return("")
  paste(deparse(x), collapse = "\n")
}
bodies <- vapply(nm, body_of, "")
for (pat in c("equidistant", "sum_to_zero", "\"disc\"", "disc\\b",
              "delta")) {
  hit <- nm[grepl(pat, bodies)]
  cat("== functions mentioning", pat, "==\n")
  print(hit)
}
cat("\n== lines mentioning disc ==\n")
for (f in nm[grepl("disc", bodies)]) {
  ln <- strsplit(bodies[[f]], "\n")[[1]]
  h <- grep("disc", ln, value = TRUE)
  h <- h[!grepl("discrete", h)]
  if (length(h)) {
    cat("--", f, "\n")
    cat(paste0("   ", trimws(h)), sep = "\n")
  }
}
cat("\n== lines mentioning delta / equidistant / sum_to_zero ==\n")
for (f in nm[grepl("delta|equidistant|sum_to_zero", bodies)]) {
  ln <- strsplit(bodies[[f]], "\n")[[1]]
  h <- grep("delta|equidistant|sum_to_zero", ln, value = TRUE)
  if (length(h)) {
    cat("--", f, "\n")
    cat(paste0("   ", trimws(h)), sep = "\n")
  }
}
