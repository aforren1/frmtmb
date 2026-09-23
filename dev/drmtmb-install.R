# Install drmTMB into the lane's private library only; see
# dev/lane-rules.md. CRAN Windows binary, no dependencies = TRUE.
LIB <- "C:/Users/adf44/source/r/drmtmb-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
options(repos = c(CRAN = "https://cloud.r-project.org"))
ap <- available.packages(type = "binary")
print(ap["drmTMB", c("Version", "Depends", "Imports", "LinkingTo",
                     "Suggests")])
deps <- tools::package_dependencies(
  "drmTMB", db = ap, which = c("Depends", "Imports", "LinkingTo"),
  recursive = TRUE)[[1]]
have <- rownames(installed.packages())
base <- rownames(installed.packages(priority = "base"))
need <- setdiff(c("drmTMB", deps), c(have, base))
print(need)
if (length(need)) {
  install.packages(need, lib = LIB, type = "binary", dependencies = FALSE)
}
print(packageVersion("drmTMB", lib.loc = LIB))
