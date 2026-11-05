# Proposed-rule mutants of fix round 2, in a source copy only.
root <- commandArgs(trailingOnly = TRUE)[1]; p <- file.path(root, "R/par-template.R")
s <- readChar(p, file.info(p)$size, useBytes = TRUE)
crlf <- grepl("\r\n", s, fixed = TRUE)
rp <- function(a, b) {
  if (crlf) { a <- gsub("\n", "\r\n", a, fixed = TRUE); b <- gsub("\n", "\r\n", b, fixed = TRUE) }
  if (!grepl(a, s, fixed = TRUE, useBytes = TRUE)) stop("missing: ", a)
  s <<- sub(a, b, s, fixed = TRUE, useBytes = TRUE)
}
M <- 'Sys.getenv("REV_MUT")'
# no stop on two agreeing restarts; a first agreement stops only an exact alias
rp("if (from_flat && same(ll, floor)) break", paste0("if (!(", M, ' %in% c("no_early", "sign4")) && from_flat && same(ll, floor)) break'))
rp("if (agrees && !(saddle_floor && same(ll, floor))) break", paste0("if (!(", M, ' %in% c("no_early", "alias_only", "sign4")) && agrees && !(saddle_floor && same(ll, floor))) break'))
rp("from_flat = !on_saddle)", paste0("from_flat = !on_saddle && (", M, ' != "alias_only" || any(vapply(held, function(w) grepl("are not identified", conditionMessage(w)), NA))))'))
# all-negative is one of the four patterns (the two-coefficient case: every sign pair)
rp("  signs <- list(rep(1, length(pars)))\n  for (i in which(owner == flip)) {\n    if (length(signs) >= 4L) break",
   paste0("  signs <- list(rep(1, length(pars)))\n  cap4 <- if (", M, ' == "sign4") 3L else 4L\n  for (i in which(owner == flip)) {\n    if (length(signs) >= cap4) break'))
rp("  lapply(signs, function(s) {\n    st <- start %||% list()",
   paste0("  if (", M, ' == "sign4" && length(which(owner == flip)) >= 2L) { neg <- rep(1, length(pars)); neg[owner == flip] <- -1; signs[[length(signs) + 1L]] <- neg }\n  lapply(signs, function(s) {\n    st <- start %||% list()'))
# a saddle beside a ridge: move every body
rp("if (!on_saddle && is.null(start) && flagged) {", paste0("if ((!on_saddle || ", M, ' == "union") && is.null(start) && flagged) {'))
rp("      if (length(bc$coefs)) {\n        bodies[[length(bodies) + 1L]] <-",
   paste0('      dup <- any(vapply(bodies, function(q) identical(q$lp[["dpar"]], lp[["dpar"]]) && identical(q$lp[["resp"]], lp[["resp"]]), NA))\n      if (length(bc$coefs) && !dup) {\n        bodies[[length(bodies) + 1L]] <-'))
writeBin(charToRaw(s), p); cat("mutated\n")
