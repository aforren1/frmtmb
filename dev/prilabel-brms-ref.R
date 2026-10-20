# Lane prilabel: brms 2.23.0's prior rows, group labels and Stan code for
# the cases this lane matches. No compile: default_prior(), stancode()
# and standata() only. Data: dev/dblbar-data.R, seed 20261013.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
dd <- dblbar_data()
dir.create(file.path(wt, "dev/prilabel-log"), showWarnings = FALSE)
sink(file.path(wt, "dev/prilabel-log/brms-ref.txt"), split = TRUE)

show_rows <- function(f, family = gaussian()) {
  cat("\n== default_prior:", deparse1(f$formula %||% f), "\n")
  p <- as.data.frame(default_prior(f, data = dd, family = family))
  print(p[p$class %in% c("sd", "cor"),
          c("prior", "class", "coef", "group", "resp", "dpar", "nlpar")],
        row.names = FALSE)
}
cases <- list(
  bf(y ~ x + (1 + x | g)),
  bf(y ~ x + (1 + x || g)),
  bf(y ~ x + (1 + x | gr(g, cor = FALSE))),
  bf(y ~ x + (1 + x | g1/g2)),
  bf(y ~ x + (1 + x | g1:g2)),
  bf(y ~ x + (1 | g1) + (1 | g1:g2)),
  bf(y ~ x + (1 + x | mm(m1, m2))),
  bf(y ~ x + (1 | g), sigma ~ (1 + x | g)),
  bf(y ~ x + (1 | p | g), sigma ~ (1 | p | g)),
  bf(y ~ x + (f | g)),
  bf(y ~ x + (1 | gr(g, id = "a")) + (0 + x | gr(g, id = "a")))
)
for (f in cases) show_rows(f)

cat("\n== Stan code: per-coefficient sd prior on g, coef x\n")
sc <- stancode(bf(y ~ x + (1 + x | g)), data = dd,
               prior = set_prior("normal(0, 0.5)", class = "sd", group = "g",
                                 coef = "x"))
cat(grep("sd_1", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")

cat("\n== Stan code: per-coefficient sd prior, || term\n")
sc <- stancode(bf(y ~ x + (1 + x || g)), data = dd,
               prior = set_prior("normal(0, 0.5)", class = "sd", group = "g",
                                 coef = "x"))
cat(grep("sd_1", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")

cat("\n== Stan code: per-coefficient sd prior, nested g1:g2, coef x\n")
sc <- stancode(bf(y ~ x + (1 + x | g1/g2)), data = dd,
               prior = set_prior("normal(0, 0.5)", class = "sd",
                                 group = "g1:g2", coef = "x"))
cat(grep("sd_[0-9]", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")

cat("\n== Stan code: per-coefficient sd prior, mm, group mmm1m2\n")
sc <- tryCatch(stancode(bf(y ~ x + (1 + x | mm(m1, m2))), data = dd,
               prior = set_prior("normal(0, 0.5)", class = "sd",
                                 group = "mmm1m2", coef = "x")),
               error = function(e) paste("ERROR:", conditionMessage(e)))
cat(grep("sd_[0-9]|ERROR", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")

cat("\n== coef without group on class sd\n")
r <- tryCatch({
  sc <- stancode(bf(y ~ x + (1 + x | g) + (1 + x | g1)), data = dd,
                 prior = set_prior("normal(0, 0.5)", class = "sd", coef = "x"))
  grep("sd_[0-9]", strsplit(sc, "\n")[[1]], value = TRUE)
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat(r, sep = "\n")

cat("\n== coef that is not in the group (sd, group g, coef z)\n")
r <- tryCatch({
  sc <- stancode(bf(y ~ x + (1 + x | g)), data = dd,
                 prior = set_prior("normal(0, 0.5)", class = "sd", group = "g",
                                   coef = "z"))
  "accepted"
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat(r, sep = "\n")

cat("\n== per-coef sd on a dpar block\n")
sc <- stancode(bf(y ~ x + (1 | g), sigma ~ (1 + x | g)), data = dd,
               prior = set_prior("normal(0, 0.5)", class = "sd", group = "g",
                                 coef = "x", dpar = "sigma"))
cat(grep("sd_[0-9]", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")

cat("\n== per-coef sd on a factor (f | g), coef fb\n")
sc <- stancode(bf(y ~ x + (f | g)), data = dd,
               prior = set_prior("normal(0, 0.5)", class = "sd", group = "g",
                                 coef = "fb"))
cat(grep("sd_[0-9]", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")

cat("\n== Student-t: two terms on g\n")
sc <- stancode(bf(y ~ x + (1 | gr(g, dist = "student")) +
                    (0 + x | gr(g, dist = "student"))), data = dd)
cat(grep("df|dfm|udf|student|sd_|z_|r_", strsplit(sc, "\n")[[1]],
         value = TRUE), sep = "\n")

cat("\n== Student-t: two terms on g, one gaussian\n")
sc <- stancode(bf(y ~ x + (1 | gr(g, dist = "student")) +
                    (0 + x | g)), data = dd)
cat(grep("df|dfm|udf|student", strsplit(sc, "\n")[[1]],
         value = TRUE), sep = "\n")

cat("\n== Student-t: one on g, one on g1\n")
sc <- stancode(bf(y ~ x + (1 | gr(g, dist = "student")) +
                    (1 | gr(g1, dist = "student"))), data = dd)
cat(grep("df|dfm|udf|student", strsplit(sc, "\n")[[1]],
         value = TRUE), sep = "\n")

cat("\n== Student-t: on g in mu and in sigma\n")
sc <- tryCatch(stancode(bf(y ~ x + (1 | gr(g, dist = "student")),
                     sigma ~ (1 | gr(g, dist = "student"))), data = dd),
               error = function(e) paste("ERROR:", conditionMessage(e)))
cat(grep("df|dfm|udf|student|ERROR", strsplit(sc, "\n")[[1]],
         value = TRUE), sep = "\n")
sink()
