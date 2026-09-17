# Lane wt-conditions: how many source lines a head rename pushes past
# 80 columns, per added width, before choosing the helper names.
roots <- c("R", file.path(list.dirs("extensions", recursive = FALSE), "R"))
over <- c(`+1` = 0L, `+4` = 0L, `+0 now` = 0L)
for (r in roots) for (f in list.files(r, "[.][Rr]$", full.names = TRUE)) {
  src <- readLines(f, warn = FALSE)
  pd <- getParseData(parse(f, keep.source = TRUE))
  hit <- pd[pd$token == "SYMBOL_FUNCTION_CALL" &
              pd$text %in% c("stop", "warning", "message"), ]
  for (ln in unique(hit$line1)) {
    k <- sum(hit$line1 == ln)
    w <- nchar(src[ln], type = "width")
    over["+0 now"] <- over["+0 now"] + (w > 80)
    over["+1"] <- over["+1"] + (w + k > 80)
    over["+4"] <- over["+4"] + (w + 4L * k > 80)
  }
}
print(over)
