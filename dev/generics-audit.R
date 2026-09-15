# Which installed package OWNS each generic frmtmb defines a rival for?
# Uses parseNamespaceFile(), not grep: nlme declares its exports in a
# multi-name block and a one-name-per-line grep missed all three.
LIB <- "C:/Users/adf44/source/r/rellib-r3"
UL  <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
SYS <- file.path(R.home(), "library")
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib", UL, SYS))

frm_ns <- parseNamespaceFile("frmtmb", LIB)
frm_exports <- frm_ns$exports

# the generics frmtmb defines itself
is_gen <- function(nm) {
  f <- tryCatch(get(nm, envir = asNamespace("frmtmb")), error = function(e) NULL)
  if (!is.function(f)) return(FALSE)
  any(grepl("UseMethod", deparse(body(f)), fixed = TRUE))
}
frm_gens <- Filter(is_gen, frm_exports)

# every installed package, by library
libs <- c(SYS, UL)
inst <- do.call(rbind, lapply(libs, function(l) {
  d <- list.dirs(l, recursive = FALSE, full.names = FALSE)
  if (!length(d)) return(NULL)
  data.frame(pkg = d, lib = l, stringsAsFactors = FALSE)
}))
inst <- inst[!duplicated(inst$pkg), ]
inst <- inst[inst$pkg != "frmtmb" & !grepl("^frmtmb[.]", inst$pkg), ]

prio <- function(pkg, lib) {
  f <- file.path(lib, pkg, "DESCRIPTION")
  if (!file.exists(f)) return(NA_character_)
  d <- tryCatch(read.dcf(f, "Priority")[1, 1], error = function(e) NA_character_)
  d
}

owners <- list()
for (i in seq_len(nrow(inst))) {
  pk <- inst$pkg[i]; lb <- inst$lib[i]
  ex <- tryCatch(parseNamespaceFile(pk, lb)$exports, error = function(e) character())
  hit <- intersect(frm_gens, ex)
  for (h in hit) owners[[h]] <- c(owners[[h]], pk)
}

cat(sprintf("frmtmb defines %d exported S3 generics\n\n", length(frm_gens)))
cat(sprintf("%-22s %-10s %s\n", "generic", "priority", "also exported by"))
cat(strrep("-", 78), "\n")
for (g in sort(frm_gens)) {
  ow <- owners[[g]]
  if (is.null(ow)) { cat(sprintf("%-22s %-10s %s\n", g, "-", "(frmtmb alone)")); next }
  pr <- sapply(ow, function(p) {
    lb <- inst$lib[match(p, inst$pkg)]
    x <- prio(p, lb); if (is.na(x)) "" else x
  })
  best <- if (any(pr == "recommended")) "recommended" else if (any(pr == "base")) "base" else ""
  cat(sprintf("%-22s %-10s %s\n", g, best, paste(ow, collapse = ", ")))
}
cat("\nDONE\n")
