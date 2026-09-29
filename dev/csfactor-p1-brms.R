# Does brms 2.23.0 build both blocks for mo(m) + cs(m), poly(x, 2) +
# cs(x) and s(x) + cs(x)? Parsing only. The new refusals are departures
# only where brms accepts the model, and the vignette has to say which.
#   Rscript dev/csfactor-p1-brms.R > dev/csfactor-log/p1-brms.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
cat("brms", as.character(packageVersion("brms")), "\n")
set.seed(1907)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$m <- factor(sample(1:4, n, TRUE), ordered = TRUE)
d$yo <- sample(1:3, n, TRUE)
show <- function(lab, f) {
  cat("\n--", lab, "--\n")
  r <- tryCatch({
    sd <- standata(f, family = sratio(), data = d)
    cat("K =", sd$K, "; X:", paste(colnames(sd$X), collapse = ", "),
        "\n")
    cat("Kcs =", sd$Kcs, "; Xcs:", paste(colnames(sd$Xcs), collapse = ", "),
        "\n")
    cat("Imo =", sd$Imo %||% NA, "; Ksp =", sd$Ksp %||% NA, "\n")
    "built"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("standata:", r, "\n")
}
`%||%` <- function(a, b) if (is.null(a)) b else a
show("mo(m) + cs(m)", bf(yo ~ mo(m) + cs(m)))
show("poly(x, 2) + cs(x)", bf(yo ~ poly(x, 2) + cs(x)))
show("s(x) + cs(x)", bf(yo ~ s(x) + cs(x)))
show("cs(m) alone", bf(yo ~ cs(m)))
cat("\ndone\n")
