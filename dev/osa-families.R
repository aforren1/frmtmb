# Lane osa: which core families can reach residuals(type = "osa"), with
# their type and which CDF slots they carry.
#   Rscript dev/osa-families.R before|after
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "after"
.libPaths(c(if (arm == "after") "C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
reg <- frmtmb:::family_registry
pm <- frmtmb:::osa_point_mass_families
out <- lapply(names(reg), function(nm) {
  f <- tryCatch(reg[[nm]](), error = function(e) NULL)
  if (is.null(f)) {
    return(data.frame(name = nm, family = NA, type = "deferred",
                      cdf_log = NA, lcdf = NA, point_mass = NA,
                      structure = NA))
  }
  data.frame(name = nm, family = f[["family"]], type = f[["type"]] %||% "",
             cdf_log = !is.null(f[["cdf_log"]]),
             lcdf = !is.null(f[["lcdf"]]),
             point_mass = f[["family"]] %in% pm,
             structure = !is.null(frmtmb:::fam_structure(f)))
})
`%||%` <- function(a, b) if (is.null(a)) b else a
tab <- do.call(rbind, out)
cat("BEGIN osa-families", arm, "\n")
print(tab, row.names = FALSE)
cat("END osa-families\n")
