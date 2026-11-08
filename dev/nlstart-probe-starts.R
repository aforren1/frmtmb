# Which deterministic start sets reach the reference on the scan's rows?
# For each row: every sign pattern of the body's coefficients at
# magnitudes 0.15 and 1, each fit from that start. Rows: the 16 short rows
# of the base scan and 3 others per family (seeded choice).
# Usage: NLSTART_ARM=before Rscript dev/nlstart-probe-starts.R <family>
fam <- commandArgs(trailingOnly = TRUE)[1]
source("dev/nlstart-common.R")
src <- readLines("dev/nlstart-scan.R")
# reuse gen(), F and bodies from the scan without its loop
eval(parse(text = src[grep("^llv <-", src):(grep("^cat[(]\"ROWHEAD", src) - 1L)]))
sc <- unlist(lapply(list.files("dev/nlstart-out/scan-before", full.names = TRUE),
                    function(f) grep("^ROW ", readLines(f), value = TRUE)))
x <- read.table(text = sub("^ROW ", "", sc), stringsAsFactors = FALSE,
                col.names = c("family", "seed", "slope", "body", "default",
                              "ref", "conv", "t_def", "nwarn", "wclass",
                              "nmsg"))
x <- x[x$family == fam, ]
x$best <- pmax(x$default, x$ref, na.rm = TRUE)
short <- which(x$best - x$default > 1e-4)
set.seed(4242)
other <- sample(setdiff(seq_len(nrow(x)), short), 3)
for (i in c(short, other)) {
  row <- x[i, ]
  b <- bodies[[row$body]]
  d <- gen(row$seed * 100 + 1, row$slope)
  fm <- do.call(bf, c(list(as.formula(paste("y ~", deparse1(b$f)))),
                      lapply(b$p, function(p) as.formula(paste(p, "~ 1"))),
                      list(nl = TRUE)))
  k <- length(b$p)
  signs <- as.matrix(expand.grid(rep(list(c(1, -1)), k)))
  for (m in c(0.15, 1)) {
    vals <- apply(signs, 1, function(s) {
      st <- setNames(s * m, paste0(b$p, "_(Intercept)"))
      llv(tt(frm(fm, family = F[[fam]], data = d, start = list(beta = st))))
    })
    cat(sprintf("PROBE %s %d %.1f %s %s m=%.2f best=%.4f default=%.4f %s\n",
                fam, row$seed, row$slope, row$body,
                if (i %in% short) "SHORT" else "ok", m, row$best,
                row$default,
                paste(apply(signs, 1, function(s) {
                  paste(ifelse(s > 0, "+", "-"), collapse = "")
                }), sprintf("%.4f", vals), sep = ":", collapse = " ")))
  }
}
