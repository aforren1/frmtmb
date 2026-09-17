# Lane wt-conditions: paste the generated blocks into
# dev/conditions-findings.md from the logs, so no count is typed.
#   Rscript dev/conditions-fill.R
md <- "dev/conditions-findings.md"
x <- readLines(md, warn = FALSE, encoding = "UTF-8")
lg <- function(f) readLines(file.path("dev/conditions-log", f), warn = FALSE)

fence <- function(lines) c("```", lines, "```")
put <- function(x, name, body) {
  b <- grep(paste0("^<!-- BEGIN GENERATED: ", name, " -->$"), x)
  e <- grep(paste0("^<!-- END GENERATED: ", name, " -->$"), x)
  stopifnot(length(b) == 1L, length(e) == 1L, e > b)
  c(x[seq_len(b)], body, x[e:length(x)])
}

rw <- lg("rewrite.txt")
i <- grep("^== counts by package ==$", rw)
ri <- lg("rewrite-inst.txt")
counts <- c(rw[(i + 1L):length(rw)],
            ri[(grep("^== counts by package ==$", ri) + 1L):length(ri)])
counts <- counts[!grepl("^OVER80|^files compared", counts) & nzchar(counts)]
parsed <- do.call(rbind, lapply(strsplit(trimws(counts), " +"), function(p) {
  data.frame(pkg = p[1L], what = paste(p[2:(length(p) - 1L)], collapse = " "),
             n = as.integer(p[length(p)]))
}))
conv <- parsed[!grepl("^exempt", parsed$what), ]
tot <- aggregate(n ~ pkg, conv, sum)
body <- c(counts, "",
          "converted per package (every non-exempt row above summed):",
          sprintf("%-16s %5d", tot$pkg, tot$n),
          sprintf("%-16s %5d", "ALL", sum(tot$n)),
          sprintf("%-16s %5d", "exempt",
                  sum(parsed$n[grepl("^exempt", parsed$what)])),
          "",
          "proof over every file, after the hand reflow and the second",
          "pass (rewrite-verify.txt):",
          lg("rewrite-verify.txt"),
          "proof after the second pass over inst/ (rewrite-inst.txt):",
          grep("^files compared", ri, value = TRUE))
x <- put(x, "rewrite counts", fence(body))

cen <- c("lane tree (census-lane.txt):", lg("census-lane.txt"), "",
         "base sources with the census added (census-absent-base.txt):",
         substr(lg("census-absent-base.txt"), 1, 76), "",
         paste0("lane plus a bare stop() in core R/utils.R and ",
                "`error = stop`"),
         "in frmtmb.eam R/ddm-shared.R (census-absent-plant.txt):",
         substr(lg("census-absent-plant.txt"), 1, 76), "",
         "lane with fit_error_context()'s stop(e) renamed, so its",
         "exemption is stale (census-absent-stale.txt):",
         substr(lg("census-absent-stale.txt"), 1, 76))
x <- put(x, "census", fence(cen))

if (file.exists("dev/conditions-log/sweep-lane.txt")) {
  sw <- lg("sweep-lane.txt")
  j <- grep("^== sweep", sw)
  x <- put(x, "sweep", fence(c("lane build (sweep-lane.txt):",
                              substr(sw[j:length(sw)], 1, 76))))
  x <- put(x, "suites", fence(c("lane build (sweep-lane.txt):",
                               substr(sw[seq_len(j - 1L)], 1, 76))))
}
if (file.exists("dev/conditions-log/sweep-base.txt")) {
  sb <- lg("sweep-base.txt")
  j <- grep("^== sweep", sb)
  k <- grep("^package +caught", sb)
  b <- grep("^<!-- END GENERATED: sweep -->$", x)
  x <- append(x, fence(c("base build, the same test files (sweep-base.txt):",
                         sb[k:(k + 9L)])), after = b - 1L)
}
con <- file(md, open = "wb")
writeLines(enc2utf8(x), con, sep = "\n", useBytes = TRUE)
close(con)
cat("filled\n")
