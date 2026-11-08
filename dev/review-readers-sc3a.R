# Review of lane readers: part A of dev/review-readers-sc3.R alone, rerun
# after its data were made positive (a lognormal part needs y > 0).
lines <- readLines("dev/review-readers-sc3.R")
stop_at <- grep("== B. hurdle_lognormal", lines, fixed = TRUE)[1] - 1L
eval(parse(text = lines[seq_len(stop_at)]))
