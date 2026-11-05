# Generates the tables of dev/refuse-findings.md from the logs, into
# dev/refuse-log/gen-tables.md, for pasting.
one <- function(f, pat) {
  x <- readLines(f, warn = FALSE)
  m <- grep(pat, x, value = TRUE)
  if (length(m)) m[1L] else NA_character_
}
out <- character(0)
add <- function(...) out[length(out) + 1L] <<- paste0(...)
add("### Sweep: frm() and brms 2.23.0 per spec (dev/refuse-sweep.R, seed 4040)")
add("")
fs <- list.files("dev/refuse-sweep-out", "tsv$", full.names = TRUE)
o <- do.call(rbind, lapply(fs, function(f) suppressWarnings(read.delim(f))))
o$f <- ifelse(grepl("^REFUSED", o$frm), "refused", "fits")
o$b <- ifelse(grepl("^builds and compiles", o$brms), "builds, compiles",
       ifelse(grepl("COMPILE FAILS", o$brms), "builds, compile fails",
              "refuses"))
add("| spec | rule | frm() | brms |")
add("|---|---|---|---|")
for (i in order(o$b != "builds, compiles", o$f, o$id)) {
  add("| ", o$id[i], " | ", gsub("|", "\\|", o$rule[i], fixed = TRUE), " | ",
      o$f[i], " | ", o$b[i], " |")
}
add("")
add("Counts: ", nrow(o), " specs; brms builds and compiles ",
    sum(o$b == "builds, compiles"), "; of those frm() refuses ",
    sum(o$b == "builds, compiles" & o$f == "refused"), ".")
add("")
add("### D8 per family (dev/refuse-mi-l1.R, seed 169 + family index)")
add("")
add("| family | L1 response const | inner grad | L1 read by mi() const | inner grad | logLik - exact, per missing | max coef shift, SE |")
add("|---|---|---|---|---|---|---|")
for (f in sub("[.]txt$", "", list.files("dev/refuse-mi-out", "txt$"))) {
  p <- file.path("dev/refuse-mi-out", paste0(f, ".txt"))
  l1 <- one(p, "^L1 uni")
  l2 <- one(p, "^L1 pred")
  la <- one(p, "^laplace error")
  co <- one(p, "^mu coef diff")
  sp <- function(s) {
    if (is.na(s)) return(c("-", "-"))
    if (grepl("ERROR", s)) return(c(substr(sub(".*ERROR ", "", s), 1, 40), "-"))
    v <- strsplit(trimws(sub("^[^:]*:", "", s)), " +")[[1L]]
    c(v[1L], v[2L])
  }
  a <- sp(l1)
  b <- sp(l2)
  per <- if (is.na(la)) "-" else sub(".*per missing ", "", la)
  sh <- if (is.na(co)) "-" else {
    format(max(abs(as.numeric(strsplit(trimws(sub(".*units:", "", co)),
                                       " +")[[1L]]))), digits = 2)
  }
  add("| ", f, " | ", a[1L], " | ", a[2L], " | ", b[1L], " | ", b[2L], " | ",
      per, " | ", sh, " |")
}
add("")
add("### brms families, mi() response (dev/refuse-brms-mi.R, seed 169)")
add("")
for (nm in c("mi-families.tsv", "mi-families-pred.tsv")) {
  p <- file.path("dev/refuse-brms-mi-out", nm)
  if (!file.exists(p)) next
  r <- read.delim(p)
  add(nm, ": builds and compiles: ",
      paste(r$family[r$built %in% TRUE & r$compiled %in% TRUE], collapse = ", "))
  add("")
  add(nm, ": refused: ",
      paste(r$family[r$built %in% FALSE], collapse = ", "))
  add("")
}
add("### Variants (dev/refuse-variants*.sh, dev/equiv-variant.R)")
add("")
add("| variant | verdict | l1_const | extra points | l1_grad | l1_error |")
add("|---|---|---|---|---|---|")
for (f in sub("[.]txt$", "", list.files("dev/refuse-variant-out", "txt$"))) {
  x <- readLines(file.path("dev/refuse-variant-out", paste0(f, ".txt")),
                 warn = FALSE)
  g <- function(k) {
    m <- grep(paste0("^ [$] ", k, " *:"), x, value = TRUE)
    if (!length(m)) return("")
    trimws(gsub('"', "", sub(paste0("^ [$] ", k, " *: *(num|chr)? *"), "", m[1L])))
  }
  add("| ", f, " | ", g("verdict"), " | ", g("l1_const"), " | ",
      g("l1_extra_const"), " | ", g("l1_grad"), " | ",
      substr(g("l1_error"), 1, 50), " |")
}
writeLines(out, "dev/refuse-log/gen-tables.md")
cat(out, sep = "\n")
