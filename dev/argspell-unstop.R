## Undo the dots guard where the method REFUSES unconditionally.
##
## A method whose whole body is one `stop()` cannot swallow anything:
## every call to it errors. Putting `frm_check_dots(...)` in front of the
## refusal only replaces a message that says what to do instead with one
## that says the argument is unknown, and it turns a legitimate
## positional argument into the wrong complaint:
## `conditional_effects(fm, "x")` on a frm_multiple() result reported
## "given 1 argument with no name" instead of "no pooled version".
##
## Run: Rscript dev/argspell-unstop.R
refusers <- c("stop", "fit_no_draws", "multiple_no_draws")
n <- 0L
for (fp in c(list.files("R", pattern = "[.]R$", full.names = TRUE),
             list.files("extensions/frmtmb.sample/R", pattern = "[.]R$",
                        full.names = TRUE))) {
  txt <- readLines(fp, warn = FALSE)
  ex <- parse(text = paste(txt, collapse = "\n"), keep.source = TRUE)
  drop <- integer()
  for (i in seq_along(ex)) {
    e <- ex[[i]]
    if (!is.call(e) || !identical(as.character(e[[1L]]), "<-")) next
    rhs <- e[[3L]]
    if (!is.call(rhs) || !identical(as.character(rhs[[1L]]), "function")) next
    b <- rhs[[3L]]
    if (!is.call(b) || !identical(as.character(b[[1L]]), "{")) next
    if (length(b) != 3L) next          # `{`, the guard, one more statement
    g <- b[[2L]]
    if (!is.call(g) ||
        !identical(as.character(g[[1L]]), "frm_check_dots")) next
    if (length(g) != 2L) next          # a bare frm_check_dots(...)
    s2 <- b[[3L]]
    if (!is.call(s2) || !as.character(s2[[1L]])[1L] %in% refusers) next
    srl <- attr(b, "srcref")
    drop <- c(drop, as.integer(srl[[2L]])[1L])
    cat("  unguarding", as.character(e[[2L]]), "in", basename(fp), "\n")
  }
  if (!length(drop)) next
  txt <- txt[-drop]
  con <- file(fp, open = "wb")
  writeLines(txt, con, sep = "\n")
  close(con)
  n <- n + length(drop)
}
cat("removed", n, "guards from unconditional refusals\n")
