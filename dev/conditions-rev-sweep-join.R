# Reviewer, lane wt-conditions (round 2): join the base and lane sweeps
# made by dev/conditions-rev-runset.sh, one row per condition an
# expect_error(), expect_warning(), expect_message() or
# expect_condition() caught, and compare message, call and the try()
# text on the same test expression.
#   Rscript dev/conditions-rev-sweep-join.R <base dir> <lane dir>
av <- commandArgs(trailingOnly = TRUE)
rd <- function(dir) {
  fs <- list.files(dir, "[.]rds$", full.names = TRUE)
  do.call(rbind, lapply(fs, readRDS))
}
b <- rd(av[1L]); l <- rd(av[2L])
cat("rows base", nrow(b), "lane", nrow(l), "\n")
# keys are per file and per deparsed expression with an occurrence count
k <- function(x) paste(x$file, x$key)
m <- merge(transform(b, k = k(b)), transform(l, k = k(l)), by = "k",
           suffixes = c(".b", ".l"))
cat("joined", nrow(m), " only base", length(setdiff(k(b), k(l))),
    " only lane", length(setdiff(k(l), k(b))), "\n")
kind <- ifelse(grepl("^expect_error", m$key.b), "error",
        ifelse(grepl("^expect_warning", m$key.b), "warning",
        ifelse(grepl("^expect_message", m$key.b), "message", "condition")))
lane_frm <- grepl("frmtmb_(error|warning|message)", m$class.l)
same_msg <- m$message.b == m$message.l
same_call <- m$call.b == m$call.l
same_try <- m$trytext.b == m$trytext.l
cat("\nby kind, lane condition is frmtmb-classed:\n")
print(table(kind, lane_frm))
cat("\nclassed rows only (converted sites), identical fields:\n")
s <- lane_frm
print(data.frame(
  kind = names(table(kind[s])),
  n = as.vector(table(kind[s])),
  msg = as.vector(tapply(same_msg[s], kind[s], sum)),
  call = as.vector(tapply(same_call[s], kind[s], sum)),
  trytext = as.vector(tapply(same_try[s], kind[s], sum))))
cat("\ndistinct converted messages compared:",
    length(unique(m$message.l[s])), "\n")
cat("classed rows whose base class was not simple* (already classed):\n")
print(table(m$class.b[s & !grepl("^simple", m$class.b)]))
bad <- s & !(same_msg & same_call & same_try)
cat("\nclassed rows that differ:", sum(bad), "\n")
if (any(bad)) {
  options(width = 250)
  print(head(m[bad, c("k", "message.b", "message.l", "call.b", "call.l")],
             40), right = FALSE)
}
cat("\nunclassed rows that differ between arms:",
    sum(!s & !(same_msg & same_call)), "\n")
nc <- !s & !(same_msg & same_call)
if (any(nc)) print(head(m[nc, c("k", "class.b", "class.l", "message.b",
                                "message.l")], 20), right = FALSE)
cat("\nnon-null calls among classed rows (lane):",
    sum(m$call.l[s] != "NULL"), "\n")
print(head(unique(m$call.l[s & m$call.l != "NULL"]), 20))
