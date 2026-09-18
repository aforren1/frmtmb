# The inverse case for the D4 pin: build the object the OTHER choice
# would produce (one copy, `data` served by a `$` method) and check that
# the three assertions fail on it. A pin is only a pin if it fails on
# the alternative.
.libPaths(c("C:/Users/adf44/source/r/adefects-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
local({
  set.seed(20260917)
  A <- diag(6); A[A == 0] <- 0.3; dimnames(A) <- list(1:6, 1:6)
  dd <- data.frame(g = factor(rep(1:6, each = 5)), x = rnorm(30))
  dd$y <- dd$x + rnorm(6)[dd$g] + rnorm(30)
  fa <- frm(y ~ x + (1 | gr(g, cov = A)), dd, data2 = list(A = A))
  # the alternative design
  alt <- fa
  alt$data <- NULL
  `$.adefects_alt` <- function(x, name) {
    if (name == "data") return(.subset2(x, "frame")[["data_frame"]])
    .subset2(x, name, exact = FALSE)
  }
  registerS3method("$", "adefects_alt", `$.adefects_alt`)
  class(alt) <- c("adefects_alt", class(alt))
  cat("alt$data is the frame          :", identical(alt$data, fa$data), "\n")
  cat("PIN 1 name in names(alt)      :",
      "data" %in% names(alt), "(want FALSE)
")
  cat("PIN 2 [[exact = TRUE]] non-NULL:",
      !is.null(alt[["data", exact = TRUE]]), "(want FALSE)\n")
  cat("PIN 3 no $ method              :",
      is.null(getS3method("$", "adefects_alt", optional = TRUE)),
      "(want FALSE)\n")
  nd <- alt; nd$data <- NULL
  cols <- length(serialize(lapply(alt$data, identity), NULL))
  delta <- length(serialize(alt, NULL)) - length(serialize(nd, NULL))
  cat("PIN 4 delta > 0.9 * cols       :", delta > 0.9 * cols,
      "(want FALSE; delta =", delta, ")\n")
})
