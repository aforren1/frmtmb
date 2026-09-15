SP <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/genrev-out"
rd <- function(tag) {
  f <- file.path(SP, paste0(tag, ".txt"))
  if (!file.exists(f)) return(NULL)
  ln <- readLines(f, warn = FALSE)
  if (!any(ln == "GENREVDONE")) {
    cat("INCOMPLETE:", tag, "\n"); cat(tail(ln, 6), sep = "\n"); return(NULL)
  }
  rows <- grep("^ROW\t", ln, value = TRUE)
  p <- do.call(rbind, strsplit(sub("^ROW\t", "", rows), "\t"))
  d <- as.data.frame(p, stringsAsFactors = FALSE)
  names(d) <- c("gen", "genowner", "methowner", "via", "gs3")
  d
}
report <- function(tag, ctrl_tag, label) {
  d <- rd(tag); k <- rd(ctrl_tag)
  if (is.null(d) || is.null(k)) return(invisible())
  m <- merge(k, d, by = "gen", suffixes = c(".ctrl", ".t"))
  # denominator: generics the control can dispatch on the class
  den <- m$via.ctrl == "class"
  lost <- den & m$via.t != "class"
  wrong <- den & m$via.t == "class" &
    m$methowner.t != m$methowner.ctrl
  silent <- den & m$via.t == "default"
  g3den <- m$gs3.ctrl != "<none>"
  g3lost <- g3den & m$gs3.t == "<none>"
  cat(sprintf("%-28s  dispatch-damage %2d/%2d   getS3method-damage %2d/%2d\n",
              label, sum(lost), sum(den), sum(g3lost), sum(g3den)))
  if (sum(lost)) cat("    lost dispatch: ",
                     paste(sort(m$gen[lost]), collapse = ", "), "\n")
  if (sum(silent)) cat("    SILENTLY answered by a .default instead: ",
                       paste(sort(m$gen[silent]), collapse = ", "), "\n")
  if (sum(wrong)) cat("    WRONG method owner: ",
                      paste(sort(paste0(m$gen[wrong], "(", m$methowner.ctrl[wrong],
                                        "->", m$methowner.t[wrong], ")")),
                            collapse = ", "), "\n")
  if (sum(g3lost)) cat("    lost getS3method: ",
                       paste(sort(m$gen[g3lost]), collapse = ", "), "\n")
  invisible()
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  for (a in args) {
    p <- strsplit(a, "=")[[1]]
    report(p[1], p[2], p[3])
  }
}
