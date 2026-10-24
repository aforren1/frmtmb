# Guardrail 3 of dev/core-boundary.md for car() and spde(), measured on
# the base commit 19b162b9 (frmtmb 0.73.1), read with `git show`.
#
# 1. Lines in R/: every top-level object of R/ whose name is car's or
#    spde's (car_*, spde_*, parse_car_call, parse_spde_call, the two
#    registry entries) with the roxygen and comment lines above it, plus
#    the lines elsewhere in R/ that name car or spde (branches, lists,
#    compat rows), counted per file.
# 2. Churn: commits touching those objects' files in the last five
#    rounds, and the commits whose diff touches a car/spde line.
#
# Run from the worktree root. Writes nothing; paste the output.
base <- "19b162b9"
files <- system2("git", c("ls-tree", "--name-only", base, "R/"),
                 stdout = TRUE)
read_base <- function(f) system2("git", c("show", paste0(base, ":", f)),
                                 stdout = TRUE)
owner <- function(nm) {
  if (grepl("^(car_|parse_car_call$|block_is_esicar$)", nm)) return("car")
  if (grepl("^(spde_|parse_spde_call$)", nm)) return("spde")
  NA_character_
}
tot <- c(car = 0L, spde = 0L)
other <- list()
for (f in files) {
  x <- read_base(f)
  ex <- parse(text = x, keep.source = TRUE)
  sr <- attr(ex, "srcref")
  prev_end <- 0L
  own_lines <- integer(0)
  for (i in seq_along(ex)) {
    e <- ex[[i]]
    nm <- NA_character_
    if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=")) {
      lhs <- e[[2]]
      nm <- if (is.name(lhs)) as.character(lhs) else deparse1(lhs)
    }
    o <- NA_character_
    if (!is.na(nm)) {
      o <- owner(nm)
      if (is.na(o) && grepl('covstruct_registry[[]["]car["]', nm)) o <- "car"
      if (is.na(o) && grepl('covstruct_registry[[]["]spde["]', nm)) o <- "spde"
    }
    s <- sr[[i]][1]; e2 <- sr[[i]][3]
    # the comment and roxygen lines directly above belong to the object
    k <- s - 1L
    while (k > prev_end && grepl("^\\s*#", x[k])) k <- k - 1L
    if (!is.na(o)) {
      tot[o] <- tot[o] + (e2 - k)
      own_lines <- c(own_lines, (k + 1L):e2)
    }
    prev_end <- e2
  }
  # lines outside those objects that name car or spde in code
  rest <- setdiff(seq_along(x), own_lines)
  code <- rest[!grepl("^\\s*#", x[rest])]
  for (w in c("car", "spde")) {
    pat <- if (w == "car") "\"car\"|car\\(|aux_car|car_type|esicar|carterms" else
      "\"spde\"|spde\\(|aux_spde|spdeterms"
    n <- sum(grepl(pat, x[code]))
    if (n) other[[w]][[f]] <- n
  }
}
cat("BASE", base, "\n")
cat("own objects, lines (code + their comments):\n"); print(tot)
for (w in names(other)) {
  cat("other R/ code lines naming", w, ":", sum(unlist(other[[w]])), "\n")
  print(unlist(other[[w]]))
}
