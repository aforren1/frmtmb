## REVIEW: what brms 2.23.0 does when an ordinal response does not use
## every category. No Stan compile: make_standata() resolves the
## threshold count and every response check, which is the question.
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(brms))
cat("brms", as.character(utils::packageVersion("brms")), "\n\n")

try_sd <- function(tag, form, fam, data) {
  cat("==", tag, "\n")
  out <- tryCatch(brms::make_standata(form, data = data, family = fam),
                  error = function(e) e, warning = function(w) w)
  if (inherits(out, "condition")) {
    cat("   ", class(out)[1L], ":", conditionMessage(out), "\n\n")
    return(invisible(NULL))
  }
  cat("    nthres =", paste(out$nthres, collapse = "/"),
      "  ncat =", paste(out$ncat, collapse = "/"), "\n")
  if (!is.null(out$Jthres)) {
    cat("    Jthres rows:", nrow(out$Jthres), " unique:\n")
    print(unique(out$Jthres))
  }
  cat("    Y range:", range(out$Y), "  table:",
      paste(table(out$Y), collapse = "/"), "\n\n")
  invisible(NULL)
}

n <- 40
set.seed(1)
x <- rnorm(n)

## 1. TOP category absent, integer response
d1 <- data.frame(x = x, y = rep(1:3, length.out = n))
try_sd("integer 1..3 (top of a 4-category scale absent), cumulative",
       brms::bf(y ~ x), cumulative(), d1)

## 2. TOP category absent, ordered factor DECLARING 4 levels
d2 <- d1
d2$y <- factor(d1$y, levels = 1:4, ordered = TRUE)
try_sd("ordered factor with 4 declared levels, only 3 used, cumulative",
       brms::bf(y ~ x), cumulative(), d2)
try_sd("same, sratio", brms::bf(y ~ x), sratio(), d2)
try_sd("same, acat", brms::bf(y ~ x), acat(), d2)

## 3. INTERIOR category absent, integer 1,3,4
d3 <- data.frame(x = x, y = rep(c(1L, 3L, 4L), length.out = n))
try_sd("integer 1,3,4 (interior category 2 absent), cumulative",
       brms::bf(y ~ x), cumulative(), d3)
try_sd("integer 1,3,4, sratio", brms::bf(y ~ x), sratio(), d3)

## 4. INTERIOR absent, ordered factor declaring 4 levels
d4 <- d3
d4$y <- factor(d3$y, levels = 1:4, ordered = TRUE)
try_sd("ordered factor 4 levels, interior level 2 unused, cumulative",
       brms::bf(y ~ x), cumulative(), d4)

## 5. thres() pinned counts
try_sd("integer 1..3 with thres(3) pinned, cumulative",
       brms::bf(y | thres(3) ~ x), cumulative(), d1)
try_sd("integer 1,3,4 with thres(3) pinned, cumulative",
       brms::bf(y | thres(3) ~ x), cumulative(), d3)

## 6. grouped thresholds with one level losing its top category
g <- factor(rep(c("a", "b"), length.out = n))
d6 <- data.frame(x = x, g = g,
                 y = ifelse(g == "a", rep(1:4, length.out = n),
                            rep(1:3, length.out = n)))
try_sd("thres(gr = g), level b tops out lower, cumulative",
       brms::bf(y | thres(gr = g) ~ x), cumulative(), d6)
d7 <- d6
d7$y <- factor(d6$y, levels = 1:4, ordered = TRUE)
try_sd("thres(gr = g) with an ordered factor of 4 declared levels",
       brms::bf(y | thres(gr = g) ~ x), cumulative(), d7)
cat("DONE rev-04\n")
