# Reviewer scratch: dump frm_compat() core-only and all-loaded.
#   Rscript dev/review-compat-dump.R <libs ;-separated> <out.rds>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(args[1], ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), find.package("frmtmb"),
    "\n")
core <- as.data.frame(frm_compat())
ext <- c("frmtmb.survival", "frmtmb.spatial", "frmtmb.latent",
         "frmtmb.sample", "frmtmb.learn", "frmtmb.eam",
         "frmtmb.coupling", "frmtmb.ode")
for (p in ext) suppressPackageStartupMessages(library(p, character.only = TRUE))
all <- as.data.frame(frm_compat())
ft <- as.data.frame(frm_compat_features())
saveRDS(list(core = core, all = all, ft = ft), args[2])
cat("core", nrow(core), "all", nrow(all), "\n")
ats <- ft$name[ft$kind == "aterm"]
only <- c("thres()", "disc", "equidistant", "sum_to_zero", "cs_pred()")
for (at in ats) {
  st <- vapply(only, function(o) frm_compat(at, o)$status, "")
  cat(sprintf("%-10s %s\n", at, paste(sprintf("%s=%s", only, st),
                                      collapse = " ")))
}
