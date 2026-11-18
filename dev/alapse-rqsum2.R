# Lane alapse, fix round 1: summarize an alapse-rq.R log by design,
# counting the local-maximum warning (ald_local) beside the convergence
# ones. gap_j is meaningful only for designs with a constant sigma and
# quantile (the quantile process); for D and E it is alt_max() from the
# fit, which finds the fit again (see dev/alapse-findings.md).
#   Rscript dev/alapse-rqsum2.R <log>
f <- commandArgs(TRUE)[1]
L <- grep(" code ", readLines(f), value = TRUE)
num <- function(tag) {
  suppressWarnings(as.numeric(sub(paste0(".*", tag, " +([-0-9.eNA+]+).*"),
                                  "\\1", L)))
}
d <- data.frame(id = sub(" .*", "", L), gap_b = num("gap_b"),
                gap_j = num("gap_j"),
                conv = grepl("conv_code|conv_gradient", L),
                local = grepl("ald_local", L),
                code = as.integer(sub(".* code ([0-9]+) .*", "\\1", L)))
for (id in unique(d$id)) {
  x <- d[d$id == id, ]
  cat(sprintf(paste0("%-6s n %2d | code 0 %2d | conv warn %2d | local warn ",
                     "%2d | gap_b max %.2e | gap_j>2e-5 %2d (max %.2e)\n"),
              id, nrow(x), sum(x$code == 0), sum(x$conv), sum(x$local),
              max(x$gap_b), sum(x$gap_j > 2e-5), max(x$gap_j)))
}
