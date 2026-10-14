p <- "extensions/frmtmb.sample/R/sample.R"
x <- paste(readLines(p), collapse = "\n")
old <- "gev_xi_blocked_text <- function(tm, of = \"\") {
  paste0(\"xi\", of, \" of gen_extreme_value() is sampled on its \",
         tm$link$name, \" scale, not through brms's tmp_xi, because of \",
         paste(tm$blocked, collapse = \" and \"), \". The support's edge \",
         \"is then a wall the sampler meets: on the data of \",
         \"dev/famuni-out/tmpxi.txt, 957 to 3388 of 4000 transitions \",
         \"diverged that way and none through tmp_xi. Check the \",
         \"divergences, or give xi a formula (xi ~ 1)\")
}"
new <- "gev_xi_blocked_text <- function(tm, of = \"\") {
  paste0(tm$dpar %||% \"xi\", of, \" of gen_extreme_value() is sampled \",
         \"on its \", tm$link$name, \" scale, not through brms's tmp_xi, \",
         \"because of \", paste(tm$blocked, collapse = \" and \"), \", and \",
         \"a prior on it is a density on xi itself. The support's edge is \",
         \"then a wall the sampler meets: on the data of \",
         \"dev/famuni-out/tmpxi.txt, 957 to 3388 of 4000 transitions \",
         \"diverged that way and none through tmp_xi. Check the divergences\")
}"
stopifnot(grepl(old, x, fixed = TRUE))
x <- sub(old, new, x, fixed = TRUE)
writeLines(strsplit(x, "\n")[[1]], p)
