## Which model_info() entries differ between arms.
b <- readRDS("dev/famlink-rev-interop-base.rds"); l <- readRDS("dev/famlink-rev-interop-lane.rds")
for (nm in names(b)) { mb <- b[[nm]]$model_info; ml <- l[[nm]]$model_info
  if (!is.list(mb)) { cat(nm, "base errored; lane link_function", ml$link_function, "family", ml$family, "is_linear", ml$is_linear, "\n"); next }
  for (k in names(ml)) if (!identical(mb[[k]], ml[[k]])) cat(nm, k, "base:", format(mb[[k]])[1], " lane:", format(ml[[k]])[1], "\n") }
