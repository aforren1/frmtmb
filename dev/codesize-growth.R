# Q5 of the codesize lane: growth of core R/ from 0.60.0 to 0.83.0.
#   Rscript dev/codesize-growth.R
# Reads dev/codesize-out/growth.tsv and rel/<v>/R (dev/codesize-growth.sh)
# and parsed.rds (for the subsystem map). Prints the blocks pasted into
# dev/codesize-findings.md.
out <- "C:/Users/adf44/source/r/frmtmb-wt-codesize/dev/codesize-out"
P <- readRDS(file.path(out, "parsed.rds"))
g <- read.delim(file.path(out, "growth.tsv"), stringsAsFactors = FALSE)
g$file <- sub("[.]R$", "", g$file)
vers <- unique(g$version)
cat("BEGIN codesize-growth\n")
tot <- aggregate(cbind(bytes, lines) ~ version + date, g, sum)
tot$files <- as.vector(table(g$version)[tot$version])
tot <- tot[match(vers, tot$version), ]
tot$d_lines <- c(NA, diff(tot$lines)); tot$d_kb <- c(NA, diff(tot$bytes)) / 1024
cat("| version | date | files | lines | KB | lines added | KB added |\n|---|---|---|---|---|---|---|\n")
cat(sprintf("| %s | %s | %d | %d | %.0f | %s | %s |\n", tot$version, tot$date,
            tot$files, tot$lines, tot$bytes / 1024,
            ifelse(is.na(tot$d_lines), "", tot$d_lines),
            ifelse(is.na(tot$d_kb), "", sprintf("%.0f", tot$d_kb))), sep = "")

a <- g[g$version == "0.60.0", ]; b <- g[g$version == "0.83.0", ]
m <- merge(a[c("file", "lines")], b[c("file", "lines")], by = "file",
           all = TRUE, suffixes = c("_060", "_083"))
m[is.na(m)] <- 0
m$delta <- m$lines_083 - m$lines_060
m$subsystem <- ifelse(m$file %in% names(P$sub_of), P$sub_of[m$file],
                      "(removed or renamed)")
m <- m[order(-m$delta), ]
cat("\nFiles that grew most, 0.60.0 to 0.83.0 (lines):\n")
x <- head(m, 25)
cat(sprintf("  %s.R: %d -> %d (+%d)%s\n", x$file, x$lines_060, x$lines_083,
            x$delta, ifelse(x$lines_060 == 0, " new", "")), sep = "")
cat(sprintf("New files since 0.60.0: %d, %d lines; removed: %d (%s)\n",
            sum(m$lines_060 == 0), sum(m$lines_083[m$lines_060 == 0]),
            sum(m$lines_083 == 0),
            paste(m$file[m$lines_083 == 0], collapse = " ")))
s <- aggregate(cbind(lines_060, lines_083, delta) ~ subsystem, m, sum)
s <- s[order(-s$delta), ]
cat("\nPer subsystem (0.83.0 map; files removed by 0.83.0 grouped apart):\n")
cat("| subsystem | 0.60.0 | 0.83.0 | added | x |\n|---|---|---|---|---|\n")
cat(sprintf("| %s | %d | %d | %d | %.2f |\n", s$subsystem, s$lines_060,
            s$lines_083, s$delta, s$lines_083 / pmax(s$lines_060, 1)),
    sep = "")

# per release, the files that grew most
cat("\nPer release, the three files that grew most:\n")
for (i in 2:length(vers)) {
  p <- g[g$version == vers[i - 1], c("file", "lines")]
  q <- g[g$version == vers[i], c("file", "lines")]
  d <- merge(p, q, by = "file", all = TRUE); d[is.na(d)] <- 0
  d$dl <- d$lines.y - d$lines.x
  d <- d[order(-d$dl), ][1:3, ]
  cat(sprintf("  %s: %s\n", vers[i],
              paste0(d$file, " +", d$dl, collapse = ", ")))
}

# the line-kind split at both ends
kinds <- function(dir) {
  fs <- list.files(dir, "[.]R$", full.names = TRUE)
  r <- c(code = 0, roxygen = 0, comment = 0, blank = 0)
  for (f in fs) {
    p <- parse(f, keep.source = TRUE)
    pd <- utils::getParseData(p); pd <- pd[pd$terminal, ]
    n <- length(readLines(f, warn = FALSE))
    k <- rep("blank", n)
    cm <- pd[pd$token == "COMMENT", ]
    k[cm$line1] <- ifelse(grepl("^#'", cm$text), "roxygen", "comment")
    nc <- pd[pd$token != "COMMENT", ]
    k[unique(unlist(Map(seq, nc$line1, nc$line2)))] <- "code"
    r <- r + table(factor(k, names(r)))
  }
  r
}
k0 <- kinds(file.path(out, "rel/0.60.0/R")); k1 <- kinds(file.path(out, "rel/0.83.0/R"))
cat("\nLine kinds 0.60.0 -> 0.83.0:\n")
cat(sprintf("  %s: %d -> %d (+%d, x%.2f)\n", names(k0), k0, k1, k1 - k0,
            k1 / k0), sep = "")
# NEWS volume per release section, as a proxy for change count
news <- readLines(file.path(out, "rel/0.83.0/NEWS.md"), warn = FALSE)
h <- grep("^# frmtmb ", news)
sec <- data.frame(version = sub("^# frmtmb ([0-9.]+).*$", "\\1", news[h]),
                  start = h, end = c(h[-1] - 1, length(news)))
sec$bullets <- vapply(seq_len(nrow(sec)), function(i)
  sum(grepl("^\\s*[*-] ", news[sec$start[i]:sec$end[i]])), 1L)
sec$news_lines <- sec$end - sec$start + 1
sec <- sec[sec$version %in% vers, ]
mm <- merge(sec, tot[c("version", "d_lines")], by = "version")
mm <- mm[order(numeric_version(mm$version)), ]
cat("\nNEWS bullets per release against lines added to R/:\n")
cat(sprintf("  %s: %d bullets, %d NEWS lines, %s R/ lines added (%s per bullet)\n",
            mm$version, mm$bullets, mm$news_lines, mm$d_lines,
            ifelse(is.na(mm$d_lines) | mm$bullets == 0, "NA",
                   round(mm$d_lines / mm$bullets))), sep = "")

# A crude keyword split of the NEWS bullets of 0.61.0 to 0.83.0: a bullet
# is "guard" when its first 400 characters name a check, warning,
# refusal, limit, boundary, standard error or start, and "feature" when
# they name brms or something new that frmtmb now takes. A bullet may be
# both. This is a proxy for what drove the growth, not a measurement of
# lines.
bul <- list()
for (i in which(sec$version %in% vers[-1])) {
  x <- news[sec$start[i]:sec$end[i]]
  st <- grep("^\\s*[*-] ", x)
  en <- c(st[-1] - 1, length(x))
  for (j in seq_along(st)) bul[[length(bul) + 1]] <- data.frame(
    version = sec$version[i],
    text = substr(paste(x[st[j]:en[j]], collapse = " "), 1, 400))
}
bul <- do.call(rbind, bul)
gp <- paste0("warn|refus|message|stops? |verdict|limit|boundary|standard ",
             "error|\\bSEs?\\b|NaN|converge|start|gradient|Hessian|check|",
             "separat|escape|silent")
fp <- "brms|\\bnew\\b|now (takes|accepts|supports|fits|gives)|adds?\\b|export"
bul$guard <- grepl(gp, bul$text, ignore.case = TRUE)
bul$feature <- grepl(fp, bul$text, ignore.case = TRUE)
bul$breaking <- grepl("BREAKING", bul$text)
bul$silent <- grepl("silent|wrong", bul$text, ignore.case = TRUE)
cat(sprintf("\nNEWS bullets 0.61.0 to 0.83.0: %d; guard keywords %d (%.0f%%), feature keywords %d (%.0f%%), both %d, BREAKING %d, 'silent'/'wrong' %d\n",
            nrow(bul), sum(bul$guard), 100 * mean(bul$guard),
            sum(bul$feature), 100 * mean(bul$feature),
            sum(bul$guard & bul$feature), sum(bul$breaking),
            sum(bul$silent)))
for (v in c("0.79.0", "0.80.0", "0.82.0", "0.83.0", "0.71.0", "0.68.0")) {
  z <- bul[bul$version == v, ]
  cat(sprintf("  %s: %d bullets, guard %d, feature %d, BREAKING %d\n", v,
              nrow(z), sum(z$guard), sum(z$feature), sum(z$breaking)))
}
cat("END codesize-growth\n")
