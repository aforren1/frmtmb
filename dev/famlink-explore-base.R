.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
             "C:/Users/adf44/source/r/pinlib",
             "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
reg <- frmtmb:::family_registry
for (nm in names(reg)) {
  f <- tryCatch(if (nm == "multinomial") reg[[nm]](K = 3) else reg[[nm]](),
                error = function(e) conditionMessage(e))
  if (is.character(f)) { cat(nm, ": ERROR ", f, "\n"); next }
  cat(sprintf("%-28s family=%-28s type=%-12s dpars=%s links=%s\n", nm,
              f[["family"]], f[["type"]], paste(f[["dpars"]], collapse = ","),
              paste(vapply(f[["links"]], `[[`, "", "name"), collapse = ",")))
}
print(names(unclass(frmtmb::student())))
str(frmtmb::student()$link, max.level = 1)
