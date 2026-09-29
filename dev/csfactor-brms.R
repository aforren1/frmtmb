# brms 2.23.0's parameterization for the two defects, parsing only.
#   A: cs(factor), cs(character), cs(numeric factor) -> Xcs columns
#   B: yo ~ x + cs(x) -> does brms build both, or drop one?
# Seed 405, the data of dev/csfactor-repro.R.
#   Rscript dev/csfactor-brms.R > dev/csfactor-log/brms.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
cat("brms", as.character(packageVersion("brms")), "\n")
set.seed(405)
n <- 500
x <- rnorm(n)
fc <- factor(sample(c("a", "b", "c"), n, TRUE))
eff <- c(a = -1, b = 0, c = 1.5)[as.character(fc)]
p1 <- plogis(-0.3 + eff); p2 <- (1 - p1) * plogis(0.5 - eff)
u <- runif(n)
d <- data.frame(x, fc, yo = ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L)))
d$fch <- as.character(fc)
d$fnum <- factor(as.integer(fc))

show <- function(lab, f, fam = sratio()) {
  cat("\n---", lab, "---\n")
  r <- tryCatch({
    sd <- standata(f, family = fam, data = d)
    cat("K =", sd$K, "; X columns:",
        paste(colnames(sd$X), collapse = ", "), "\n")
    cat("Kcs =", sd$Kcs, "; Xcs columns:",
        paste(colnames(sd$Xcs), collapse = ", "), "\n")
    "built"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("standata:", r, "\n")
  sc <- tryCatch(stancode(f, family = fam, data = d),
                 error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("stancode parameters / model lines mentioning b or bcs:\n")
  if (inherits(sc, "character") && length(sc) == 1L &&
        grepl("^ERROR", sc)) {
    cat(sc, "\n")
  } else {
    ln <- strsplit(as.character(sc), "\n")[[1L]]
    keep <- grep("bcs|vector\\[Kc?\\] b|mu \\+=|Xcs|target \\+=", ln)
    cat(paste(ln[keep], collapse = "\n"), "\n")
  }
  cat("get_prior:\n")
  print(tryCatch(get_prior(f, family = fam, data = d),
                 error = function(e) conditionMessage(e)))
}

show("A1 cs(fc), factor", bf(yo ~ cs(fc)))
show("A2 cs(fch), character", bf(yo ~ cs(fch)))
show("A3 cs(fnum), factor with numeric-looking levels", bf(yo ~ cs(fnum)))
show("B1 yo ~ x + cs(x)", bf(yo ~ x + cs(x)))
show("B2 yo ~ fc + cs(fc)", bf(yo ~ fc + cs(fc)))
show("B3 yo ~ cs(x) alone (control)", bf(yo ~ cs(x)))
cat("\ndone\n")
