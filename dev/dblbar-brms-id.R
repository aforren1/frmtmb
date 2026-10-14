# brms 2.23.0 on |ID| spellings next to `||` and gr(cor = FALSE), and on
# gr(cov = ) under `||` (lane dblbar). Output: dev/dblbar-log/brms-id.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages(library(brms))
dd <- dblbar_data()
A <- diag(24) * 0.5 + 0.5
dimnames(A) <- list(levels(dd$g), levels(dd$g))
sink(file.path(wt, "dev/dblbar-log/brms-id.txt"), split = TRUE)
show <- function(txt) {
  cat("\n==== ", txt, "\n", sep = "")
  form <- eval(str2lang(txt))
  r <- tryCatch({
    pr <- as.data.frame(default_prior(form, data = dd, data2 = list(A = A)))
    print(pr[pr$class %in% c("sd", "cor"), c("class", "coef", "group",
                                              "dpar")], row.names = FALSE)
    sc <- strsplit(stancode(form, data = dd, data2 = list(A = A)),
                   "\n")[[1]]
    cat(trimws(grep("r_[0-9]+_[0-9]+ = |r_[0-9]+ = |L_[0-9]+;|udf",
                    sc, value = TRUE)), sep = "\n")
  }, error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
}
show("bf(y ~ x + (1 + x | p || g))")
show("bf(y ~ x + (1 + x | p | gr(g, cor = FALSE)))")
show("bf(y ~ x + (1 + x | gr(g, cor = FALSE)), sigma ~ (1 | gr(g, id = 'q')))")
show(paste0("bf(y ~ x + (1 + x | p | gr(g, cor = FALSE)), ",
            "sigma ~ (1 | p | g))"))
show(paste0("bf(y ~ x + (1 + x || g), sigma ~ (1 | g))"))
show("bf(y ~ x + (1 + x || gr(g, cov = A)))")
show("bf(y ~ f + (1 + f || gr(g, cov = A)))")
show("bf(y ~ f + (1 + f || gr(g, dist = 'student')))")
show("bf(y ~ f + (1 + f | gr(g, cor = FALSE, dist = 'student')))")
show("bf(y ~ x + (1 + x | gr(g, cor = FALSE)))")
sink()
