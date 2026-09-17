## Compare dev/famlink-rev-interop-{base,lane}.rds element by element.
b <- readRDS("dev/famlink-rev-interop-base.rds")
l <- readRDS("dev/famlink-rev-interop-lane.rds")
for (nm in names(b)) for (k in names(b[[nm]])) {
  same <- identical(b[[nm]][[k]], l[[nm]][[k]])
  if (!same) {
    cat(sprintf("%-10s %-20s DIFFERS\n", nm, k))
    if (k %in% c("family_link", "model_info")) {
      cat("  base:", paste(capture.output(str(b[[nm]][[k]], max.level = 1))[1:6], collapse = "\n  "), "\n")
      cat("  lane:", paste(capture.output(str(l[[nm]][[k]], max.level = 1))[1:6], collapse = "\n  "), "\n")
    } else {
      print(all.equal(b[[nm]][[k]], l[[nm]][[k]]))
    }
  } else cat(sprintf("%-10s %-20s identical\n", nm, k))
}
