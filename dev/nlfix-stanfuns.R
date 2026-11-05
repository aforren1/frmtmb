# Which functions do brms 2.23.0's own vignettes, Rd examples, R examples
# and test suite call in nonlinear bodies? Reads dev/brms-suite/brms, parses
# every bf(...)/nlf(...) call it can and lists the names in call position
# of the body of each nl formula. Writes nothing; prints name, count and
# where.
root <- "dev/brms-suite/brms"
files <- c(list.files(file.path(root, "vignettes"), "[.](Rmd|ltx)$",
                      full.names = TRUE),
           list.files(file.path(root, "man"), "[.]Rd$", full.names = TRUE),
           list.files(file.path(root, "R"), "[.]R$", full.names = TRUE),
           list.files(file.path(root, "tests/testthat"), "[.]R$",
                      full.names = TRUE))
base_ops <- c("+", "-", "*", "/", "^", "(", "~", "bf", "nlf", "lf", "|", "||",
              "brmsformula", "%*%", "<=", ">=", "<", ">", "==", "[", "[[",
              ":", "$", "c", "::", "=", "<-")
hits <- list()
walk <- function(e, f) {
  if (!is.call(e)) return(invisible())
  fn <- e[[1]]
  if (is.name(fn) && as.character(fn) %in% c("bf", "brmsformula", "nlf")) {
    # the first formula argument, or any formula naming an nlf parameter
    is_nl <- as.character(fn) == "nlf" || !is.null(e$nl)
    if (is_nl) {
      rhs <- NULL
      a1 <- if (length(e) >= 2L) e[[2]]
      if (is.call(a1) && identical(a1[[1]], as.name("~"))) rhs <- a1[[length(a1)]]
      if (!is.null(rhs)) {
        nms <- setdiff(setdiff(all.names(rhs), all.vars(rhs)), base_ops)
        for (n in nms) hits[[length(hits) + 1L]] <<- data.frame(
          name = n, file = sub(paste0(root, "/"), "", f, fixed = TRUE),
          body = substr(paste(deparse(rhs), collapse = " "), 1, 70))
      }
    }
  }
  for (i in seq_along(e)) {
    a <- tryCatch(e[[i]], error = function(x) NULL)
    if (!is.null(a) && is.call(a)) walk(a, f)
  }
  invisible()
}
chunks <- function(f) {
  ln <- readLines(f, warn = FALSE)
  if (grepl("[.]Rd$", f)) {
    # examples blocks: Rd markup removed enough to parse by call
    ln <- gsub("\\\\dontrun\\{|\\\\donttest\\{|\\\\%", "", ln)
  }
  if (grepl("[.]R$", f)) ln <- sub("^#' ?", "", ln)
  txt <- paste(ln, collapse = "\n")
  # every bf( or nlf( call: parse from its start, extending until it parses
  pos <- gregexpr("\\b(bf|nlf|brmsformula)\\(", txt)[[1]]
  out <- list()
  for (p in pos[pos > 0]) {
    tail_txt <- substr(txt, p, p + 800)
    # balanced parens
    ch <- strsplit(tail_txt, "")[[1]]
    d <- cumsum((ch == "(") - (ch == ")"))
    end <- which(d == 0L & cumsum(ch == "(") > 0L)[1L]
    if (is.na(end)) next
    e <- tryCatch(parse(text = substr(tail_txt, 1, end))[[1]],
                  error = function(x) NULL)
    if (!is.null(e)) out[[length(out) + 1L]] <- e
  }
  out
}
for (f in files) for (e in chunks(f)) walk(e, f)
h <- do.call(rbind, hits)
cat("files read:", length(files), "\n")
cat("calls found in nl bodies:\n")
tab <- sort(table(h$name), decreasing = TRUE)
print(tab)
cat("\nby name:\n")
for (n in names(tab)) {
  s <- h[h$name == n, ]
  s <- s[!duplicated(s$body), ]
  cat(sprintf("%-14s %s\n", n, paste(sprintf("%s: %s", s$file, s$body)[
    seq_len(min(3, nrow(s)))], collapse = " ;; ")))
}
# the R-side mirrors brms carries, so a body runs off the tape too
nh <- file.path(root, "R/numeric-helpers.R")
cat("\nbrms numeric-helpers.R defines:\n")
cat(paste(sub(" <- function.*", "", grep("^[a-zA-Z_.0-9]+ <- function",
                                         readLines(nh), value = TRUE)),
          collapse = ", "), "\n")
