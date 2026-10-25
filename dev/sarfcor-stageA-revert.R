# Undo a one-off part-B edit script on the stage-A copy of a file, so
# that the stage-A tree carries part A alone. Usage:
#   Rscript dev/sarfcor-stageA-revert.R <edit-script> <file-in-repo>
args <- commandArgs(TRUE)
sp <- "C:/Users/adf44/AppData/Local/Temp/12/claude/c--Users-adf44-source-r-frmtmb/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad/sarfcor-stageA"
src <- readLines(args[1])
target <- file.path(sp, args[2])
src <- sub("^p <- .*$", paste0("p <- \"", target, "\""), src)
code <- paste(src, collapse = "\n")
env <- new.env()
# the edit's rep1(old, new), run backwards
env$rep1_rev <- TRUE
code <- sub("txt <<- sub(old, new, txt, fixed = TRUE)",
            "txt <<- sub(new, old, txt, fixed = TRUE)", code, fixed = TRUE)
code <- sub("n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))",
            "n <- lengths(regmatches(txt, gregexpr(new, txt, fixed = TRUE)))",
            code, fixed = TRUE)
eval(parse(text = code), envir = env)
cat("reverted", args[1], "on", target, "\n")
