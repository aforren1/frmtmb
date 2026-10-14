# Lane dblbar: (1 + f || gr(g, dist = "student")) on the base and the
# lane. brms 2.23.0 gives the block ONE latent scale per level, shared
# by its three columns (dev/dblbar-log/brms-id.txt: udf_1, dfm_1 times
# every r_1_k). Data seed 20261013. Usage: Rscript ... base|lane
arm <- commandArgs(TRUE)[1]
lib <- if (identical(arm, "base")) character(0) else
  "C:/Users/adf44/source/r/wt-dblbar-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages(library(frmtmb))
dd <- dblbar_data()
fr <- frm(bf(y ~ f + (1 + f || gr(g, dist = "student"))), data = dd,
          family = gaussian(), dry_run = "frame")
for (b in fr$re_blocks) cat(arm, ": block", b$covstruct, "dim", b$dim,
                            paste(b$cnms, collapse = ","), "nu", b$dist_nu,
                            "\n")
fit <- suppressWarnings(frm(bf(y ~ f + (1 + f || gr(g, dist = "student"))),
                            data = dd, family = gaussian()))
cat(arm, ": names", grep("^sd_|^cor_", variables(fit), value = TRUE), "\n")
cat(arm, ": logLik", format(as.numeric(logLik(fit)), digits = 12), "\n")
# the numeric term (fix round 1, review item 1)
frn <- frm(bf(y ~ x + (1 + x || gr(g, dist = "student"))), data = dd,
           family = gaussian(), dry_run = "frame")
for (b in frn$re_blocks) cat(arm, ": numeric block", b$covstruct, "dim",
                             b$dim, paste(b$cnms, collapse = ","), "\n")
fitn <- suppressWarnings(frm(bf(y ~ x + (1 + x || gr(g, dist = "student"))),
                             data = dd, family = gaussian()))
cat(arm, ": numeric names", grep("^sd_|^cor_", variables(fitn), value = TRUE),
    "\n")
cat(arm, ": numeric logLik", format(as.numeric(logLik(fitn)), digits = 12),
    "\n")
