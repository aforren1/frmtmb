# Lane prilabel: more brms 2.23.0 group labels (no compile). Data:
# dev/dblbar-data.R, seed 20261013, plus member weights w1, w2.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages(library(brms))
dd <- dblbar_data()
set.seed(7)
dd$w1 <- runif(nrow(dd))
dd$w2 <- 1 - dd$w1
dd$g3 <- factor(sample(c("r1", "r2"), nrow(dd), TRUE))
sink(file.path(wt, "dev/prilabel-log/brms-ref2.txt"), split = TRUE)
cat("brms", format(packageVersion("brms")), "\n")
show_rows <- function(f) {
  cat("\n== default_prior:", deparse1(f$formula), "\n")
  p <- tryCatch(as.data.frame(default_prior(f, data = dd)),
                error = function(e) {
                  cat("ERROR:", conditionMessage(e), "\n"); NULL
                })
  if (!is.null(p)) {
    print(p[p$class %in% c("sd", "cor"),
            c("class", "coef", "group", "dpar")], row.names = FALSE)
  }
  sd <- tryCatch(standata(f, data = dd), error = function(e) NULL)
  if (!is.null(sd)) cat("N_ groups:", grep("^N_[0-9]", names(sd), value = TRUE),
                        "\n")
}
cases <- list(
  bf(y ~ x + (1 | mm(m1, m2, weights = cbind(w1, w2)))),
  bf(y ~ x + (1 | mm(m1, m2, scale = FALSE))),
  bf(y ~ x + (1 | g1/g2/g3)),
  bf(y ~ x + (1 + x || g1/g2)),
  bf(y ~ x + (1 | p | g1/g2)),
  bf(y ~ x + (1 | g2:g1)),
  bf(y ~ x + (1 | gr(g1/g2))),
  bf(y ~ x + diag(1 + x | g1/g2)),
  bf(y ~ x + (1 | g) + (0 + x | g)),
  bf(y ~ x + (1 | g1 / g2) + (0 + x | g1:g2))
)
for (f in cases) show_rows(f)
cat("\n== brms ranef level labels and Stan names via frame_re\n")
for (f in list(bf(y ~ x + (1 | g1/g2/g3)),
               bf(y ~ x + (1 | mm(m1, m2, weights = cbind(w1, w2)))))) {
  re <- brms:::frame_re(brms:::brmsterms(f), data = dd)
  print(as.data.frame(re)[, c("id", "group", "gn", "coef", "cor")])
}
cat("\n== per-coef sd prior on Intercept of g1:g2:g3\n")
sc <- stancode(bf(y ~ x + (1 | g1/g2/g3)), data = dd,
               prior = set_prior("normal(0, 0.5)", class = "sd",
                                 group = "g1:g2:g3", coef = "Intercept"))
cat(grep("sd_[0-9]", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
cat("\n== sd prior by group, explicit mm with weights label\n")
r <- tryCatch({
  sc <- stancode(bf(y ~ x + (1 | mm(m1, m2, weights = cbind(w1, w2)))),
                 data = dd,
                 prior = set_prior("normal(0, 0.5)", class = "sd",
                                   group = "mmm1m2"))
  grep("sd_[0-9]", strsplit(sc, "\n")[[1]], value = TRUE)
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat(r, sep = "\n")
cat("\n== sd prior on group g2:g1 for a g1/g2 model\n")
r <- tryCatch({
  sc <- stancode(bf(y ~ x + (1 | g1/g2)), data = dd,
                 prior = set_prior("normal(0, 0.5)", class = "sd",
                                   group = "g2:g1"))
  grep("sd_[0-9]", strsplit(sc, "\n")[[1]], value = TRUE)
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat(r, sep = "\n")
cat("\n== sd prior group mm(m1, m2)\n")
r <- tryCatch({
  sc <- stancode(bf(y ~ x + (1 | mm(m1, m2))), data = dd,
                 prior = set_prior("normal(0, 0.5)", class = "sd",
                                   group = "mm(m1, m2)"))
  grep("sd_[0-9]", strsplit(sc, "\n")[[1]], value = TRUE)
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat(r, sep = "\n")
cat("\n== sd coef prior with no group on a one-group model\n")
r <- tryCatch({
  sc <- stancode(bf(y ~ x + (1 + x | g)), data = dd,
                 prior = set_prior("normal(0, 0.5)", class = "sd",
                                   coef = "x"))
  grep("sd_[0-9]", strsplit(sc, "\n")[[1]], value = TRUE)
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat(r, sep = "\n")
cat("\n== sd group prior and coef prior together, which wins\n")
sc <- stancode(bf(y ~ x + (1 + x | g)), data = dd,
               prior = c(set_prior("normal(0, 2)", class = "sd", group = "g"),
                         set_prior("normal(0, 0.5)", class = "sd",
                                   group = "g", coef = "x")))
cat(grep("sd_[0-9]", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
cat("\n== Student-t on two terms with an |ID| key\n")
sc <- tryCatch(stancode(bf(y ~ x + (1 | p | gr(g, dist = "student")) +
                    (0 + x | p | gr(g, dist = "student"))), data = dd),
               error = function(e) paste("ERROR:", conditionMessage(e)))
cat(grep("df|dfm|udf|ERROR|r_1", strsplit(sc, "\n")[[1]], value = TRUE),
    sep = "\n")
cat("\n== Student-t: (1 + x | g) and (0 + z | g)\n")
sc <- stancode(bf(y ~ x + (1 + x | gr(g, dist = "student")) +
                    (0 + z | gr(g, dist = "student"))), data = dd)
cat(grep("df|dfm|udf|r_[0-9] |r_[0-9]_", strsplit(sc, "\n")[[1]],
         value = TRUE), sep = "\n")
cat("\n== Student-t: df prior rows\n")
p <- as.data.frame(default_prior(bf(y ~ x + (1 | gr(g, dist = "student")) +
                    (0 + x | gr(g, dist = "student"))), data = dd))
print(p[p$class %in% c("df"), c("prior", "class", "coef", "group")])
sink()
