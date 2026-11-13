# Lane arspat, backlog 4320: how marginaleffects 1.0.0 builds the newdata
# it hands get_predict(), read off its namespace.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ns <- asNamespace("marginaleffects")
print(grep("sanit|pad", ls(ns), value = TRUE))
for (nm in c("get_predictions", "pad")) {
  if (exists(nm, ns)) {
    cat("\n=====", nm, "\n")
    print(get(nm, ns))
  }
}
