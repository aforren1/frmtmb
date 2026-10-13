# Reviewer: brms's Stan code for mo() interactions and group-level mo()
# terms, with no compile. Usage: Rscript dev/moparity-rev-brms-defects.R
# [brmslib]; brmslib first on the path picks the brms version.
args <- commandArgs(TRUE)
blib <- if (length(args)) args[[1]] else character(0)
.libPaths(c(blib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "from", find.package("brms"), "\n")
set.seed(1)
n <- 90
d <- data.frame(x = sample(0:3, n, TRUE), f = factor(rep(c("a", "b", "c"), 30)),
                f2 = factor(rep(c("a", "b"), 45)), z = rnorm(n), w = rnorm(n),
                g = factor(rep(1:9, 10)), xl = rnorm(n) > 0)
d$y <- rnorm(n); d$x2 <- sample(0:2, n, TRUE)
d$fa <- as.numeric(d$f == "a"); d$fb <- as.numeric(d$f == "b")
d$fc <- as.numeric(d$f == "c")
show <- function(lab, fo) {
  cat("\n==", lab, ":", deparse1(fo), "\n")
  r <- tryCatch({
    sc <- stancode(fo, data = d)
    sd <- standata(fo, data = d)
    lines <- strsplit(sc, "\n")[[1]]
    keep <- grep("mo\\(simo|Csp|simplex|r_1_2|bsp", lines, value = TRUE)
    cat(paste0("  ", trimws(keep)), sep = "\n")
    cat("  standata Csp:", paste(grep("^Csp", names(sd), value = TRUE),
                                 collapse = " "), "\n")
    for (cn in grep("^Csp", names(sd), value = TRUE)) {
      cat("   ", cn, "first 6:", head(sd[[cn]], 6), "\n")
    }
    pr <- default_prior(fo, data = d)
    cat("  bsp coefs:", paste(pr$coef[pr$class == "b" & nzchar(pr$coef)],
                              collapse = " "), "\n")
    cat("  simo coefs:", paste(pr$coef[pr$class == "simo"], collapse = " "),
        "\n")
  }, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))
  invisible(r)
}
show("factor 3 levels, main * ", y ~ mo(x) * f)
show("factor 3 levels, colon only", y ~ mo(x):f)
show("factor 3 levels, f + colon", y ~ f + mo(x):f)
show("factor 2 levels, main *", y ~ mo(x) * f2)
show("indicator spelling", y ~ fb + fc + mo(x) + mo(x):fb + mo(x):fc)
show("group-level with population", y ~ mo(x) + (mo(x) | g))
show("group-level alone", y ~ (mo(x) | g))
show("group-level alone plus another special", y ~ mo(x2) + (mo(x) | g))
show("z + mo(x):z", y ~ z + mo(x):z)
show("z * mo(x)", y ~ z * mo(x))
show("mo(x) * z", y ~ mo(x) * z)
show("mo(x) * z * w", y ~ mo(x) * z * w)
show("logical multiplier", y ~ mo(x) * xl)
show("f2 + mo(x):f2", y ~ f2 + mo(x):f2)
show("factor group-level", y ~ mo(x) * f2 + (mo(x) * f2 | g))
# cox baseline names
dc <- data.frame(t = rexp(60), x = rnorm(60))
sc <- stancode(t ~ x, data = dc, family = cox())
cat("\n== cox baseline:\n")
cat(paste0("  ", trimws(grep("sbhaz", strsplit(sc, "\n")[[1]], value = TRUE))),
    sep = "\n")
pr <- default_prior(t ~ x, data = dc, family = cox())
print(pr[pr$class == "sbhaz", ])
