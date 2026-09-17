## Roxygenise and install one package of this worktree into the lane's
## private library, and nowhere else.
##   Rscript dev/brmsnames-install.R core|sample|learn|eam|latent [noroxy]
av <- commandArgs(trailingOnly = TRUE)
source("dev/brmsnames-libs.R")
brmsnames_libs("lane")
LIB <- "C:/Users/adf44/source/r/brmsnames-lib"
path <- switch(av[1L], core = ".", sample = "extensions/frmtmb.sample",
               learn = "extensions/frmtmb.learn",
               eam = "extensions/frmtmb.eam",
               latent = "extensions/frmtmb.latent", stop("which package"))
if (!"noroxy" %in% av) roxygen2::roxygenise(path)
rc <- system2(file.path(R.home("bin"), "R"),
              c("CMD", "INSTALL", "--no-multiarch", "--no-test-load",
                paste0("--library=", LIB), shQuote(path)))
cat("R CMD INSTALL exit", rc, "\n")
if (rc != 0) quit(status = 1)
