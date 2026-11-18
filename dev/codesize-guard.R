# Q4 (second half) of the codesize lane: how much of core's code is guard
# and check logic (SE check, limit rules, escapes, starts, refusals,
# input checks) against modeling code. A name-and-file classification, so
# an estimate: a function counts as guard when its file is one of the
# fit-end check files or its name matches the pattern below.
#   Rscript dev/codesize-guard.R
# Reads dev/codesize-out/functions.tsv.
out <- "C:/Users/adf44/source/r/frmtmb-wt-codesize/dev/codesize-out"
funs <- read.delim(file.path(out, "functions.tsv"), stringsAsFactors = FALSE)
gfiles <- c("se-check.R", "se-reparam.R", "fit-end.R", "fit-open-ends.R",
            "fit-restart.R", "weights-unbounded.R", "ald-laplace.R",
            "mixture-starts.R", "start-ls.R", "gp-start.R")
pat <- paste0("check|guard|refus|limit|escape|start|restart|warn|^stop_|",
              "verdict|nonest|^se_|^jc_|validate|assert|verify|settle|",
              "headroom|probe|_ok$|^is_ok|reject|forbid|^ensure|safe|",
              "_message$|_msg$|said|unsettled|conflict|edge|nocurv|flat")
nm <- tolower(funs$name)
funs$kind <- ifelse(funs$file %in% gfiles, "guard (fit-end check files)",
             ifelse(grepl(pat, nm), "guard (by name, other files)",
                    "other"))
k <- aggregate(cbind(n = 1, code = funs$code, comment = funs$comment,
                     rox = funs$rox) ~ kind, funs, sum)
cat("BEGIN codesize-guard\n")
cat(sprintf("pattern: %s\n", pat))
cat("| kind | functions | code | other comments | roxygen |\n|---|---|---|---|---|\n")
cat(sprintf("| %s | %d | %d | %d | %d |\n", k$kind, k$n, k$code, k$comment,
            k$rox), sep = "")
cat(sprintf("guard share of function code: %.1f%%\n",
            100 * sum(k$code[k$kind != "other"]) / sum(k$code)))
g <- funs[funs$kind == "guard (by name, other files)", ]
pf <- aggregate(cbind(n = 1, code = g$code) ~ file, g, sum)
pf <- pf[order(-pf$code), ][1:15, ]
cat("Guard-named functions in other files, top 15 files:\n")
cat(sprintf("  %s: %d functions, %d code lines\n", pf$file, pf$n, pf$code),
    sep = "")
top <- funs[order(-funs$code), ][1:25, ]
cat("\nLargest 25 functions (code lines):\n")
cat(sprintf("  %s (%s:%d) %d code, %d comment, %d roxygen; %s\n", top$name,
            top$file, top$first, top$code, top$comment, top$rox, top$kind),
    sep = "")
cat("END codesize-guard\n")
