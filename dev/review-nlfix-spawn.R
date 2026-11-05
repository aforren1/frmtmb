# Item 10: can this session start a child R process? Every 10 s for 6 minutes, through system2 with
# stdout = TRUE, as the failing tests do. Its parent chain is stopped from outside after 2 minutes.
rs <- file.path(R.home("bin"), "Rscript.exe")
for (i in 1:36) {
  t0 <- Sys.time()
  out <- tryCatch(suppressWarnings(system2(rs, c("-e", shQuote("cat('CHILDOK')")), stdout = TRUE, stderr = TRUE)), error = function(e) paste("ERR", conditionMessage(e)))
  st <- attr(out, "status")
  cat(format(t0, "%H:%M:%S"), "out=[", paste(out, collapse = " "), "] status=", if (is.null(st)) "0" else st, "\n"); flush.console()
  Sys.sleep(10)
}
