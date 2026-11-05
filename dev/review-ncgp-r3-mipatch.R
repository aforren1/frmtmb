# Short check: the exact change to refuse's mi_retape(), applied to the trial merge only.
p <- "C:/Users/adf44/source/r/rev-ncgp-src/mtree3/R/frame.R"
s <- readLines(p)
i <- grep("fit$obj <- RTMB::MakeADFun(build_objective(fit$frame), fit$estimates,", s, fixed = TRUE)
stopifnot(length(i) == 1L)
stopifnot(grepl("random = random, map = fit$frame[[\"map\"]],", s[i + 1L], fixed = TRUE), grepl("silent = TRUE)", s[i + 2L], fixed = TRUE))
s[i] <- "  fit$obj <- gp_tape_scope(fit$frame, RTMB::MakeADFun("
s[i + 1L] <- "    build_objective(fit$frame), gp_tape_order(fit$frame, fit$estimates),"
s[i + 2L] <- "    random = random, map = fit$frame[[\"map\"]], silent = TRUE))"
writeLines(s, p)
cat("mi_retape patched\n")
