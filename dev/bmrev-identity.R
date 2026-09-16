## Reviewer: the pmin identity, in both directions, and what the Rd
## says about `...` on the two diagnostics.
.libPaths(c("C:/Users/adf44/source/r/bmrev-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample)); q(library(posterior))
ds <- readRDS("dev/stan-cache/bmrev-draws.rds")
s <- summary(ds)
nr <- neff_ratio(ds)
keep <- rownames(s)
N <- posterior::ndraws(ds)
pm <- unname(pmin(s[, "Bulk_ESS"], s[, "Tail_ESS"]))
lhs <- unname(nr[keep])

cat("ndraws(ds): ", N, "\n\n")
cat("-- DIVISION direction, which is what NEWS and the new test say --\n")
cat("neff_ratio(ds)[keep]  ==  pmin(Bulk_ESS, Tail_ESS) / ndraws(ds)\n")
cat("  identical(): ", identical(lhs, pm / N), "\n")
cat("  max |diff|:  ", format(max(abs(lhs - pm / N)), digits = 17), "\n")

cat("\n-- MULTIPLICATION direction --\n")
cat("pmin(Bulk_ESS, Tail_ESS)  ==  neff_ratio(ds) * ndraws(ds)\n")
cat("  identical(): ", identical(pm, lhs * N), "\n")
cat("  max |diff|:  ", format(max(abs(pm - lhs * N)), digits = 17), "\n")
u <- max(abs(pm - lhs * N) / (.Machine$double.eps * abs(pm)))
cat("  max gap in ulp: ", format(u, digits = 4), "\n")
cat("  per-row: \n")
print(data.frame(var = keep, pmin = pm, neff_x_N = lhs * N,
                 diff = pm - lhs * N), digits = 17, row.names = FALSE)

cat("\n-- what `...` is documented to do on rhat()/neff_ratio() --\n")
rd <- file.path(system.file(package = "frmtmb.sample"), "help",
                "frmtmb.sample")
db <- tools::Rd_db("frmtmb.sample",
                   lib.loc = "C:/Users/adf44/source/r/bmrev-lib")
one <- db[["draws-diagnostics.Rd"]]
txt <- paste(utils::capture.output(print(one)), collapse = "\n")
lines <- strsplit(txt, "\n")[[1L]]
i <- grep("^\\\\usage|item\\{\\.\\.\\.\\}|\\\\arguments", lines)
cat(paste(grep("\\.\\.\\.", lines, value = TRUE), collapse = "\n"), "\n")
cat("\n-- and the rendered Two different pars rules section --\n")
out <- tempfile()
tools::Rd2txt(one, out = out)
tt <- readLines(out, warn = FALSE)
j <- grep("pars", tt)
cat(paste(tt[seq(max(1, min(j) - 2), min(length(tt), max(j) + 2))],
          collapse = "\n"), "\n")
cat("DONE\n")
