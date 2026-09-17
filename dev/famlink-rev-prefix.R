## Reviewer check for lane wt-famlink, priority 2: which `$name` reads
## in core, the extensions and their tests would now ERROR on some
## family object, because `name` is not a field of it but is the prefix
## of exactly one field (the case base R `$` silently partial-matched)?
##
## Field names are collected from every registry family plus
## the extension families listed below, as built (not fitted). Every
## `$name` token in R/, extensions/*/R and all tests is then checked
## against every family's field set. Output: candidate (name, family,
## would-hit, file:line) rows; each has to be read to see whether that
## read can meet that family.
ARM <- "lane"
source("dev/famlink-rev-common.R")
for (p in c("frmtmb.eam", "frmtmb.learn", "frmtmb.latent", "frmtmb.spline",
            "frmtmb.ode", "frmtmb.coupling")) {
  suppressMessages(library(p, character.only = TRUE))
}
fams <- list()
reg <- frmtmb:::family_registry
for (nm in names(reg)) {
  f <- tryCatch(reg[[nm]](), error = function(e) NULL)
  if (!is.null(f)) fams[[nm]] <- f
}
try_add <- function(label, expr) {
  f <- tryCatch(expr, error = function(e) {
    cat("could not build", label, ":", conditionMessage(e), "\n"); NULL })
  if (inherits(f, "frmtmb_family")) fams[[label]] <<- f
}
try_add("mixture", mixture(gaussian(), gaussian()))
try_add("cumulative_probit", cumulative("probit"))
exports <- c(ls("package:frmtmb.eam"), ls("package:frmtmb.learn"),
             ls("package:frmtmb.latent"), ls("package:frmtmb.spline"),
             ls("package:frmtmb.ode"), ls("package:frmtmb.coupling"))
for (e in exports) {
  fn <- get(e)
  if (!is.function(fn)) next
  f <- tryCatch(suppressWarnings(suppressMessages(fn())), error = function(err) NULL)
  if (inherits(f, "frmtmb_family")) fams[[e]] <- f
}
cat("families built:", length(fams), "\n")
cat(paste(names(fams), collapse = " "), "\n")
files <- c(list.files("R", "[.]R$", full.names = TRUE),
           Sys.glob("extensions/*/R/*.R"),
           list.files("tests/testthat", "[.]R$", full.names = TRUE),
           Sys.glob("extensions/*/tests/testthat/*.R"))
rows <- list()
for (f in files) {
  lines <- readLines(f, warn = FALSE)
  m <- gregexpr("\\$[A-Za-z_][A-Za-z0-9_.]*", lines)
  for (i in seq_along(lines)) {
    toks <- regmatches(lines[i], m[i])[[1]]
    if (!length(toks)) next
    for (t in unique(substring(toks, 2))) {
      for (fn in names(fams)) {
        nms <- names(fams[[fn]])
        if (t %in% nms) next
        hit <- nms[startsWith(nms, t)]
        if (length(hit) == 1L) {
          rows[[length(rows) + 1L]] <- data.frame(
            name = t, family = fn, hit = hit, where = paste0(f, ":", i),
            code = trimws(substr(lines[i], 1, 90)))
        }
      }
    }
  }
}
tab <- do.call(rbind, rows)
agg <- unique(tab[, c("name", "hit")])
cat("distinct (name, would-hit) pairs:", nrow(agg), "\n")
print(agg, row.names = FALSE)
utils::write.table(tab, "dev/famlink-rev-prefix.tsv", sep = "\t",
                   row.names = FALSE, quote = FALSE)
