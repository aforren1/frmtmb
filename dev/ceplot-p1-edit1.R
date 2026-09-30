# Lane ceplot punch 1: one-off edit of R/predict.R for B1, keying the
# "old_levels" choice by grouping factor and carrying which rows are
# new in each design part.
f <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot/R/predict.R"
s <- paste(readLines(f), collapse = "\n")
rep1 <- function(s, old, new) {
  n <- lengths(regmatches(s, gregexpr(old, s, fixed = TRUE)))
  if (n != 1L) stop("expected one match of: ", old, " got ", n)
  sub(old, new, s, fixed = TRUE)
}
s <- rep1(s,
"      re_parts[[length(re_parts) + 1L]] <- list(bk = bk, comp = comp,
                                                mm = mm, j = j,
                                                new_key = gv)",
"      # is_new: the rows that read an unseen level of this block, which
      # a by-level's block does not share with the rows it skips
      re_parts[[length(re_parts) + 1L]] <- list(bk = bk, comp = comp,
                                                mm = mm, j = j,
                                                new_key = gv,
                                                is_new = is_new &
                                                  is.na(j))")
s <- rep1(s,
"    list(bk = bk, comp = comp,
         mm = mmk,
         j = new_level_pick_apply(pick, bk, as.character(gv[[k]]),
                                  iw$J[, k], is_new[, k]),
         new_key = as.character(gv[[k]]))",
"    jk <- new_level_pick_apply(pick, bk, as.character(gv[[k]]),
                               iw$J[, k], is_new[, k])
    list(bk = bk, comp = comp, mm = mmk, j = jk,
         new_key = as.character(gv[[k]]),
         is_new = is_new[, k] & is.na(jk))")
writeLines(strsplit(s, "\n", fixed = TRUE)[[1]], f)
cat("edited\n")
