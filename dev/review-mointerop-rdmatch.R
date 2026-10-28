.libPaths(c("C:/Users/adf44/source/r/wt-mointerop-lib", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
db <- tools::Rd_db("frmtmb")
for (t in c("frmtmb-marginaleffects.Rd", "frm.Rd")) {
  a <- paste(capture.output(tools::Rd2txt(db[[t]], options = list(underline_titles = FALSE))), collapse = "\n")
  b <- paste(capture.output(tools::Rd2txt(tools::parse_Rd(file.path("man", t)), options = list(underline_titles = FALSE))), collapse = "\n")
  cat(t, "installed == tree:", identical(a, b), "\n")
  if (!identical(a, b)) {
    al <- strsplit(a, "\n")[[1]]; bl <- strsplit(b, "\n")[[1]]
    cat("only installed:\n"); print(setdiff(al, bl)); cat("only tree:\n"); print(setdiff(bl, al))
  }
}
