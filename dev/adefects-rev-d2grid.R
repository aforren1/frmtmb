source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")
suppressPackageStartupMessages(library(frmtmb))

# The design of dev/adefects-repro-base.R, so logLik values compare.
set.seed(20260917)
d <- data.frame(g1 = rep(c(1, 1, 2, 2), each = 6),
                g2 = rep(c(1, 2, 1, 2), each = 6),
                t = c(1:6, 1:6, 1:6, 7:12))
d$x <- rnorm(24)
d$g <- interaction(d$g1, d$g2)
d$y <- d$x + as.numeric(stats::filter(rnorm(24), 0.6, "recursive"))
d$gf1 <- factor(d$g1)
d$gf2 <- factor(d$g2)
d$g3 <- rep(c(1, 2), 12)
d$g.dot <- d$g1
d[[".lead"]] <- d$g1
d[["my var"]] <- d$g1
d$t.dot <- d$t
d[["my t"]] <- d$t
gv <- d$g1              # a variable holding the codes, not a column
gfac <- factor(d$g1)    # a variable holding a factor

verdict <- function(expr) {
  r <- tryCatch({
    v <- suppressWarnings(suppressMessages(eval(expr)))
    paste0("OK ", v)
  }, error = function(e) paste0("REFUSE: ",
                                gsub("[\r\n]+", " ", conditionMessage(e))))
  substr(r, 1, 150)
}

fm <- function(txt) {
  f <- stats::as.formula(paste("y ~ x +", txt))
  verdict(bquote(round(as.numeric(stats::logLik(frm(.(f), d))), 5)))
}

cases <- c(
  # --- gr grammar
  "ar(t, g1, cov = TRUE)",
  "ar(t, factor(g1), cov = TRUE)",
  "ar(t, interaction(g1, g2), cov = TRUE)",
  "ar(t, g1:g2, cov = TRUE)",
  "ar(t, gf1:gf2, cov = TRUE)",
  "ar(t, g1:g2:g3, cov = TRUE)",
  "ar(t, g.dot, cov = TRUE)",
  "ar(t, .lead, cov = TRUE)",
  "ar(t, `my var`, cov = TRUE)",
  "ar(t, d$g1, cov = TRUE)",
  "ar(t, gv, cov = TRUE)",
  "ar(t, gfac, cov = TRUE)",
  "ar(t, get(\"g1\"), cov = TRUE)",
  "ar(t, g1/g2, cov = TRUE)",
  "ar(t, g1 + g2, cov = TRUE)",
  "ar(t, g1 * g2, cov = TRUE)",
  "ar(t, (g1), cov = TRUE)",
  "ar(t, 1, cov = TRUE)",
  "ar(t, \"g1\", cov = TRUE)",
  "ar(t, NA, cov = TRUE)",
  # --- time grammar
  "ar(x + t, g, cov = TRUE)",
  "ar(t.dot, g, cov = TRUE)",
  "ar(`my t`, g, cov = TRUE)",
  "ar(d$t, g, cov = TRUE)",
  "ar(as.integer(t), g, cov = TRUE)",
  "ar(1, g, cov = TRUE)",
  "ar(\"t\", g, cov = TRUE)",
  "ar(NA, g, cov = TRUE)",
  "ar(gr = g, cov = TRUE)",
  # --- other constructors, same grammar
  "cosy(gr = g1:g2)",
  "cosy(gr = factor(g1))",
  "unstr(t, g1:g2)",
  "ma(t, g1:g2, cov = TRUE)",
  "arma(t, g1:g2, cov = TRUE)"
)

cat("CASE\tFRMTMB\n")
for (cs in cases) cat(cs, "\t", fm(cs), "\n", sep = "")
