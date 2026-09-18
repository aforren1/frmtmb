# Does R's serializer preserve sharing of an ordinary list element?
addr <- function(x) sub("^.*(@[0-9a-fx]+).*$", "\1",
                        utils::capture.output(.Internal(inspect(x)))[1L])
set.seed(1)
mf <- data.frame(a = rnorm(20000), b = rnorm(20000))
mf2 <- data.frame(a = mf$a + 0, b = mf$b + 0)   # equal, separate SEXP
l0 <- list(x = mf)
l1 <- list(x = mf, y = mf)                      # SHARED
l2 <- list(x = mf, y = mf2)                     # equal copy
td <- tempdir()
for (nm in c("l0", "l1", "l2")) {
  p <- file.path(td, paste0(nm, ".rds"))
  saveRDS(get(nm), p, compress = FALSE)
  cat(sprintf("%s raw bytes = %d  objsize = %d\n", nm, file.size(p),
              as.numeric(utils::object.size(get(nm)))))
}
cat("in memory: addr(l1$x) == addr(l1$y):",
    identical(addr(l1$x), addr(l1$y)), "\n")
r <- readRDS(file.path(td, "l1.rds"))
cat("after readRDS: addr equal:", identical(addr(r$x), addr(r$y)),
    " identical():", identical(r$x, r$y), "\n")
cat("one frame raw size:", file.size(file.path(td, "l0.rds")) -
      file.size(file.path(td, "l0.rds")) + 0, "\n")
cat("delta l1 - l0 =", file.size(file.path(td, "l1.rds")) -
      file.size(file.path(td, "l0.rds")),
    "  delta l2 - l0 =", file.size(file.path(td, "l2.rds")) -
      file.size(file.path(td, "l0.rds")), "\n")
