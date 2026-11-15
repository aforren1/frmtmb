# Review of lane predraw2: bitwise compare of predraw2-review-c lane/base fingerprints.
a <- readRDS("dev/reviews/predraw2-review-c-lane.rds")
b <- readRDS("dev/reviews/predraw2-review-c-base.rds")
for (id in names(a)) {
  x <- a[[id]]$fp; y <- b[[id]]$fp
  cat(sprintf("%-10s identical %-5s max|diff| row means %.3g first row %.3g\n", id,
              identical(x, y), max(abs(x$rm - y$rm)), max(abs(x$c1 - y$c1))))
}
