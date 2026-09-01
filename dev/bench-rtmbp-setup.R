## one-time setup for the RTMBp (kaskr experimental parallel RTMB) benchmark.
## RTMBp and RTMB are distinct namespaces and coexist.
options(repos = c(kaskr = "https://kaskr.r-universe.dev",
                  CRAN = "https://cloud.r-project.org"))
ap <- available.packages()
cat("--- RTMBp availability ---\n")
if ("RTMBp" %in% rownames(ap)) print(ap["RTMBp", c("Version", "Repository", "NeedsCompilation")])
cat("--- binary availability ---\n")
apb <- tryCatch(available.packages(type = "win.binary"), error = function(e) NULL)
if (!is.null(apb) && "RTMBp" %in% rownames(apb))
  print(apb["RTMBp", c("Version", "Repository")]) else cat("no win.binary entry\n")

if (!requireNamespace("RTMBp", quietly = TRUE)) {
  install.packages("RTMBp")
}
cat("--- installed versions ---\n")
for (p in c("RTMB", "RTMBp", "TMB", "Matrix", "lme4", "frmtmb")) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  cat(sprintf("  %-8s %s\n", p, v))
}
