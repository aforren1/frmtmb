# Lane alapse: summarize an alapse-rq.R log by design.
#   Rscript dev/alapse-rqsum.R <log>
f <- commandArgs(TRUE)[1]
L <- grep(" code ", readLines(f), value = TRUE)
num <- function(tag) {
  suppressWarnings(as.numeric(sub(paste0(".*", tag, " +([-0-9.eNA+]+).*"),
                                  "\\1", L)))
}
d <- data.frame(id = sub(" .*", "", L), gap_b = num("gap_b"),
                gap_j = num("gap_j"), db = num("[|] db"),
                warn = grepl("conv_code|conv_gradient", L),
                code = as.integer(sub(".* code ([0-9]+) .*", "\\1", L)))
for (id in unique(d$id)) {
  x <- d[d$id == id, ]
  cat(sprintf(paste0("%-6s n %2d | code!=0 %2d | conv warn %2d | ",
                     "gap_b>2e-5 %2d (max %.2e) | gap_j>2e-5 %2d ",
                     "(max %.2e, median %.2e) | short warned %d of %d, ",
                     "at-rq silent %d of %d\n"),
              id, nrow(x), sum(x$code != 0), sum(x$warn),
              sum(x$gap_b > 2e-5), max(x$gap_b), sum(x$gap_j > 2e-5),
              max(x$gap_j), median(x$gap_j),
              sum(x$warn & x$gap_b > 2e-5), sum(x$gap_b > 2e-5),
              sum(!x$warn & x$gap_b <= 2e-5), sum(x$gap_b <= 2e-5)))
}
